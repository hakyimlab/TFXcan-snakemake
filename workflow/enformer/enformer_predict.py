#!/usr/bin/env python3
# Author: Temi
# Date: Tue Sep 29 2026
# Description: Personalized Enformer predictions for many individuals at many loci.
#   For each locus and individual: build both haplotype sequences from a phased VCF,
#   predict with Enformer, average the bins around the locus (per haplotype), and sum
#   the two haplotypes. Writes the same files as TFXcan-snakemake's
#   enformer_predict.py -> enformer_merge.py -> enformer_process.py chain:
#       {output_basename}.metadata.tsv   (locus, individual)
#       {output_basename}.matrix.h5.gz   (pandas HDF, key 'matrix', one row per locus x individual, 5313 columns)
# Usage:
#   python3 enformer_predict.py --config config.yaml                          # predict on all visible GPUs (or CPU), then merge
#   python3 enformer_predict.py --config config.yaml --shard 0 --n-shards 10 --no-merge   # one piece of a job array
#   python3 enformer_predict.py --config config.yaml --locus chr1_100_101 --no-merge      # just one locus (one job per locus)
#   python3 enformer_predict.py --config config.yaml --merge-only             # merge after all shards finish

import argparse, os, sys, time, glob, queue, threading, subprocess
import multiprocessing as mp
import numpy as np
import yaml

SEQ_LEN = 393216
N_BINS = 896
BIN_SIZE = 128
N_TRACKS = 5313

DEFAULTS = {
    'n_individuals': -1,      # -1 = everyone in the individuals file
    'pad_bins': 1,            # bins added on each side of the two middle bins; 1 -> average of bins 446..449
    'aggregation': 'mean',    # how to combine bins within a haplotype: mean (aggByMean) or sum (aggBySum)
    'batch_size': 1,          # sequences per Enformer call; >1 is no faster on A40s and changes values slightly (~0.02%)
    'devices': 'auto',        # auto (all visible GPUs, else CPU), cpu, or a list of GPU ids e.g. [0, 1]
    'cpu_threads': None,      # threads for TensorFlow on CPU; None = all available
    'loci_dir': None,         # where per-locus results go; None = {output_dir}/loci
}

# one-hot lookup by ASCII code; A C G T order and N = 0.25 as in kipoiseq.one_hot_dna
ONE_HOT = np.full((256, 4), np.nan, dtype=np.float32)
for _i, _b in enumerate('ACGT'):
    ONE_HOT[ord(_b)] = np.eye(4, dtype=np.float32)[_i]
ONE_HOT[ord('N')] = 0.25
BASE_CODE = {b: i for i, b in enumerate('ACGT')}
EYE4 = np.eye(4, dtype=np.float32)


def log(msg, worker=None):
    prefix = f'[{time.strftime("%H:%M:%S")}]' + (f' [worker {worker}]' if worker is not None else '')
    print(f'{prefix} {msg}', flush=True)


# ---------------------------------------------------------------- inputs

def load_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    for k, v in DEFAULTS.items():
        cfg.setdefault(k, v)
    required = ['loci_file', 'individuals', 'vcf_pattern', 'fasta_file', 'model_path', 'output_dir', 'output_basename']
    missing = [k for k in required if k not in cfg]
    if missing:
        sys.exit(f'ERROR - config is missing: {", ".join(missing)}')
    if cfg['loci_dir'] is None:
        cfg['loci_dir'] = os.path.join(cfg['output_dir'], 'loci')
    if not 1 <= int(cfg['batch_size']) <= 8:
        # 16 sequences exceed a cuDNN tensor-size limit inside Enformer and abort the process
        sys.exit(f"ERROR - batch_size must be between 1 and 8, not {cfg['batch_size']}")
    if cfg['aggregation'] not in ('mean', 'sum'):
        sys.exit(f"ERROR - aggregation must be mean or sum, not {cfg['aggregation']}")
    return cfg


def read_loci(path):
    # first whitespace-separated field of each line, e.g. chr1_16984695_16984696; order kept, duplicates dropped
    loci = []
    with open(path) as f:
        for line in f:
            if line.strip():
                loci.append(line.split()[0])
    return list(dict.fromkeys(loci))


def read_individuals(path, n):
    with open(path) as f:
        ids = [line.split()[0] for line in f if line.strip()]
    return ids if n == -1 else ids[:n]


def vcf_for(cfg, chrom):
    return cfg['vcf_pattern'].format(chrom=chrom)


# ---------------------------------------------------------------- sequences

def locus_window(locus):
    # same as kipoiseq.Interval(chrom, start, end).resize(393216) on the + strand
    chrom, start, end = locus.split('_')
    start, end = int(start), int(end)
    center = (start + end) // 2 + (start + end) % 2
    win_start = center - SEQ_LEN // 2
    return chrom, win_start, win_start + SEQ_LEN


def bins_to_aggregate(pad_bins):
    # the window is centered on the locus, so the locus sits where the two middle bins (447, 448) meet;
    # with 128-bp bins it is hard to say which of the two holds the variant, so take both and pad each side.
    # pad_bins = 1 -> bins 446..449 (4 bins); returned as a [start, stop) slice
    return N_BINS // 2 - 1 - pad_bins, N_BINS // 2 + 1 + pad_bins


def reference_one_hot(fasta, chrom, win_start, win_end):
    # returns (393216, 4) float32, or a reason string if the window cannot be used
    if chrom not in fasta.keys():
        return f'{chrom} not in fasta'
    chrom_len = len(fasta[chrom])
    seq = fasta[chrom][max(win_start, 0):min(win_end, chrom_len)].seq.upper()
    seq = 'N' * max(-win_start, 0) + seq + 'N' * max(win_end - chrom_len, 0)
    x = ONE_HOT[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
    if np.isnan(x).any():
        return 'non-ACGTN (e.g. R or Y) base in reference'
    if np.all(x == 0.25):
        return 'reference window is all N'
    return x


def read_haplotypes(vcf_path, samples, chrom, win_start):
    # returns (positions within window, base codes of shape (n_variants, n_samples, 2), VCF sample order, n missing calls)
    import cyvcf2, warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # IDs absent from the VCF; reported once at startup
        vcf = cyvcf2.VCF(vcf_path, samples=samples, lazy=True, threads=2)
    positions, codes, n_missing = [], [], 0
    # same query as the old code: 1-based window start to one past the window end
    for v in vcf(f'{chrom}:{win_start + 1}-{win_start + SEQ_LEN + 1}'):
        idx = v.POS - (win_start + 1)
        if idx < 0 or idx >= SEQ_LEN:
            continue
        alleles = [v.REF] + v.ALT
        if any(a not in BASE_CODE for a in alleles):
            continue  # not a biallelic/multiallelic SNP with plain bases
        gt = v.genotype.array()[:, :2]
        missing = gt < 0
        if missing.any():
            n_missing += int(missing.sum())
            gt = np.where(missing, 0, gt)  # missing calls get the REF base
        allele_codes = np.array([BASE_CODE[a] for a in alleles], dtype=np.uint8)
        positions.append(idx)
        codes.append(allele_codes[gt])
    order = list(vcf.samples)
    vcf.close()
    if not positions:
        return np.zeros(0, dtype=np.int64), np.zeros((0, len(order), 2), dtype=np.uint8), order, n_missing
    return np.array(positions, dtype=np.int64), np.stack(codes), order, n_missing


def build_batch(ref, positions, hap_codes):
    # hap_codes: (batch, n_variants) base codes -> (batch, 393216, 4) one-hot sequences
    x = np.broadcast_to(ref, (hap_codes.shape[0],) + ref.shape).copy()
    if positions.size:
        x[np.arange(hap_codes.shape[0])[:, None], positions[None, :], :] = EYE4[hap_codes]
    return x


# ---------------------------------------------------------------- prediction

def load_model(model_path, device, cpu_threads=None):
    # device: GPU id (int) or 'cpu'; must run before tensorflow is imported in this process
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    os.environ['CUDA_VISIBLE_DEVICES'] = '' if device == 'cpu' else str(device)
    import tensorflow as tf
    for gpu in tf.config.list_physical_devices('GPU'):
        tf.config.experimental.set_memory_growth(gpu, True)
    if device == 'cpu' and cpu_threads:
        tf.config.threading.set_intra_op_parallelism_threads(int(cpu_threads))
    model = tf.saved_model.load(model_path).model
    return model, tf


def predict_locus(model, tf, ref, positions, codes, bins, aggregation, batch_size):
    # codes: (n_variants, n_samples, 2) -> (n_samples, 5313) = agg(hap1 bins) + agg(hap2 bins)
    n_samples = codes.shape[1]
    haps = codes.transpose(1, 2, 0).reshape(n_samples * 2, -1)  # row 2i = hap1 of sample i, 2i+1 = hap2
    out = np.empty((haps.shape[0], N_TRACKS), dtype=np.float32)
    for i in range(0, haps.shape[0], batch_size):
        x = build_batch(ref, positions, haps[i:i + batch_size])
        pred = model.predict_on_batch(tf.constant(x))['human']  # (batch, 896, 5313)
        window = pred[:, bins[0]:bins[1], :].numpy().astype(np.float32)
        out[i:i + x.shape[0]] = window.mean(axis=1) if aggregation == 'mean' else window.sum(axis=1)
    return out.reshape(n_samples, 2, N_TRACKS).sum(axis=1)


def locus_file(cfg, locus):
    return os.path.join(cfg['loci_dir'], f'{locus}.npz')


def invalid_file(cfg, locus):
    return os.path.join(cfg['loci_dir'], f'{locus}.invalid')


def worker(cfg, loci, device, worker_id):
    import pyfaidx
    individuals = read_individuals(cfg['individuals'], cfg['n_individuals'])
    fasta = pyfaidx.Fasta(cfg['fasta_file'])
    log(f'{len(loci)} loci on {"CPU" if device == "cpu" else f"GPU {device}"}; loading model', worker_id)
    model, tf = load_model(cfg['model_path'], device, cfg['cpu_threads'])

    # read sequences for the next locus while the current one is predicted
    todo = queue.Queue(maxsize=2)

    def producer():
        for locus in loci:
            chrom, win_start, win_end = locus_window(locus)
            ref = reference_one_hot(fasta, chrom, win_start, win_end)
            if isinstance(ref, str):
                todo.put((locus, ref))
                continue
            vcf_path = vcf_for(cfg, chrom)
            if not os.path.exists(vcf_path):
                todo.put((locus, f'no VCF at {vcf_path}'))
                continue
            positions, codes, samples, n_missing = read_haplotypes(vcf_path, individuals, chrom, win_start)
            todo.put((locus, (ref, positions, codes, samples, n_missing)))
        todo.put(None)

    threading.Thread(target=producer, daemon=True).start()
    done = 0
    while (item := todo.get()) is not None:
        locus, payload = item
        if isinstance(payload, str):
            log(f'WARNING - skipping {locus}: {payload}', worker_id)
            with open(invalid_file(cfg, locus), 'w') as f:
                f.write(payload + '\n')
            continue
        ref, positions, codes, samples, n_missing = payload
        if n_missing:
            log(f'WARNING - {locus}: {n_missing} missing genotype calls set to REF', worker_id)
        tic = time.perf_counter()
        values = predict_locus(model, tf, ref, positions, codes, bins_to_aggregate(cfg['pad_bins']),
                               cfg['aggregation'], cfg['batch_size'])
        # write then rename, so a killed job never leaves a half-written locus behind
        tmp = locus_file(cfg, locus) + '.tmp.npz'
        np.savez(tmp, values=values, individuals=np.array(samples))
        os.replace(tmp, locus_file(cfg, locus))
        done += 1
        secs = time.perf_counter() - tic
        log(f'{locus}: {len(samples)} individuals, {positions.size} variants, {secs:.0f}s '
            f'({secs / (2 * len(samples)):.3f}s per sequence) [{done}/{len(loci)}]', worker_id)
    fasta.close()


def visible_gpus():
    env = os.environ.get('CUDA_VISIBLE_DEVICES')
    if env is not None:
        return [g for g in env.split(',') if g.strip() != '']
    try:
        out = subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True, check=True).stdout
        return [str(i) for i, line in enumerate(out.splitlines()) if line.startswith('GPU')]
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []


def run_predictions(cfg, shard, n_shards, only_loci=None):
    os.makedirs(cfg['loci_dir'], exist_ok=True)
    loci = read_loci(cfg['loci_file'])
    if only_loci:
        unknown = [l for l in only_loci if l not in loci]
        if unknown:
            sys.exit(f'ERROR - not in {cfg["loci_file"]}: {", ".join(unknown)}')
        loci = [l for l in loci if l in only_loci]
    loci = loci[shard::n_shards]
    remaining = [l for l in loci if not (os.path.exists(locus_file(cfg, l)) or os.path.exists(invalid_file(cfg, l)))]
    log(f'{len(loci)} loci in this shard; {len(loci) - len(remaining)} already done; {len(remaining)} to predict')
    if not remaining:
        return

    if cfg['devices'] == 'cpu':
        devices = ['cpu']
    elif cfg['devices'] == 'auto':
        devices = visible_gpus() or ['cpu']
    else:
        devices = [str(d) for d in cfg['devices']]
    log(f'devices: {", ".join("CPU" if d == "cpu" else f"GPU {d}" for d in devices)}')

    ctx = mp.get_context('spawn')  # fresh processes so each can pin its own GPU before importing tensorflow
    procs = [ctx.Process(target=worker, args=(cfg, remaining[i::len(devices)], d, i))
             for i, d in enumerate(devices) if remaining[i::len(devices)]]
    for p in procs:
        p.start()
    for p in procs:
        p.join()
    failed = [i for i, p in enumerate(procs) if p.exitcode != 0]
    if failed:
        sys.exit(f'ERROR - worker(s) {failed} failed; rerun the same command to resume')


# ---------------------------------------------------------------- merge

def merge(cfg):
    import pandas as pd
    loci = read_loci(cfg['loci_file'])
    values, metadata, invalid, missing = [], [], [], []
    for locus in loci:
        if os.path.exists(locus_file(cfg, locus)):
            d = np.load(locus_file(cfg, locus))
            values.append(d['values'])
            metadata.append(pd.DataFrame({'locus': locus, 'individual': d['individuals']}))
        elif os.path.exists(invalid_file(cfg, locus)):
            invalid.append(locus)
        else:
            missing.append(locus)
    if missing:
        sys.exit(f'ERROR - {len(missing)} loci have no predictions yet (e.g. {missing[0]}); finish predicting before merging')
    if invalid:
        log(f'WARNING - {len(invalid)} loci could not be predicted and are left out; see {cfg["loci_dir"]}/*.invalid')

    os.makedirs(cfg['output_dir'], exist_ok=True)
    basename = os.path.join(cfg['output_dir'], cfg['output_basename'])
    df_metadata = pd.concat(metadata, ignore_index=True)
    df_metadata.to_csv(f'{basename}.metadata.tsv', sep='\t', index=False)
    pd.DataFrame(np.concatenate(values, axis=0)).to_hdf(f'{basename}.matrix.h5.gz', key='matrix', mode='w', complevel=9)
    log(f'wrote {basename}.metadata.tsv and {basename}.matrix.h5.gz: '
        f'{df_metadata.shape[0]} rows ({len(values)} loci x {df_metadata.individual.nunique()} individuals)')


def main():
    parser = argparse.ArgumentParser(description='Personalized Enformer predictions: both haplotypes, bins averaged, haplotypes summed.')
    parser.add_argument('--config', required=True, help='YAML config (see config.example.yaml)')
    parser.add_argument('--shard', type=int, default=0, help='this piece of the loci, for job arrays (0-based)')
    parser.add_argument('--n-shards', type=int, default=1, help='total pieces the loci are split into')
    parser.add_argument('--locus', nargs='+', default=None, help='only predict these loci (must be in loci_file), e.g. one locus per job')
    parser.add_argument('--no-merge', action='store_true', help='only predict; merge later with --merge-only')
    parser.add_argument('--merge-only', action='store_true', help='only merge finished predictions into the output files')
    args = parser.parse_args()

    cfg = load_config(args.config)
    individuals = read_individuals(cfg['individuals'], cfg['n_individuals'])
    loci = read_loci(cfg['loci_file'])
    log(f'{len(loci)} loci x {len(individuals)} individuals requested; output in {cfg["output_dir"]}')
    first_vcf = vcf_for(cfg, loci[0].split('_')[0])
    if os.path.exists(first_vcf):
        import cyvcf2
        in_vcf = set(cyvcf2.VCF(first_vcf).samples)
        absent = [i for i in individuals if i not in in_vcf]
        if absent:
            log(f'WARNING - {len(absent)} of {len(individuals)} individuals are not in the VCF and will be skipped (e.g. {", ".join(absent[:3])})')

    if not args.merge_only:
        run_predictions(cfg, args.shard, args.n_shards, args.locus)
    if not args.no_merge:
        merge(cfg)


if __name__ == '__main__':
    main()
