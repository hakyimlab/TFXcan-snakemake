
# input summary statistic -> split by chromosomes
checkpoint process_summary_statistics:
    output: directory(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'))
    params:
        runmeta = runmeta,
        jobname = '{phenotype}',
        diag_file = os.path.join(DATA_DIR, 'diagnostics', f'{{phenotype}}.gwas_diagnostics.summary'),
        reference_annotations = REFERENCE_ANNOTATIONS,
        input_sumstats = lambda wildcards: os.path.join(INPUT_SUMSTATS, run_list[wildcards.phenotype]),
        pthreshold = config['processing']['GWAS_pvalue_threshold']
    message: "working on {wildcards}" 
    resources:
        mem_cpu=8,
        cpu_task=8
    shell:
        """
        Rscript workflow/process/process_summary_statistics.R --summary_stats_file {params.input_sumstats} --output_folder {output} --annotation_file {params.reference_annotations} --diagnostics_file {params.diag_file} --pvalue_threshold {params.pthreshold}
        """

# select either the top n snps or the most significant snp per LD region
checkpoint select_top_snps: 
    input: lambda wildcards: checkpoints.process_summary_statistics.get(phenotype = wildcards.phenotype).output[0]
    output:
        directory(os.path.join(FILTERING_DIR, '{phenotype}'))
    params:
        runmeta = runmeta,
        jobname = '{phenotype}',
        input_sumstats = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'chr{{}}.sumstats.txt.gz'),
        ld_blocks = config['processing']['LD_blocks'],
        chroms = collect_chromosomes,
        diag_file = os.path.join(DATA_DIR, 'diagnostics', f'{{phenotype}}.chr{{}}.topSNPs_diagnostics.summary'),
        selection_method = config['processing']['selection_method'],
        select_n_snps = config['processing']['select_n_snps'],
        rank_by = config['processing'].get('rank_by', 'pval')
    message: "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=12
    shell:
        """
        module load parallel;
        printf "%s\\n" {params.chroms} | parallel -j 12 "Rscript workflow/process/select_top_snps.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {output} --phenotype {wildcards.phenotype} --diagnostics_file {params.diag_file} --selection_method {params.selection_method} --select_n_snps {params.select_n_snps} --rank_by {params.rank_by}"
        """

# collect the result above into one data
# a checkpoint: once it runs, snakemake reads the chosen loci and creates one Enformer job per locus
checkpoint collect_top_snps_results:
    input: lambda wildcards: checkpoints.select_top_snps.get(**wildcards).output[0]
    output:
        filtered_sumstats = os.path.join(COLLECTION_DIR, '{phenotype}.filteredGWAS.topSNPs.txt.gz'),
        enformer_loci = os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.txt')
    params:
        runmeta = runmeta,
        jobname = '{phenotype}',
        # optional genome-wide cap on loci, applied after per-chromosome selection
        limit_flag = f"--limit_number_of_loci {config['processing']['limit_number_of_loci']}" if config['processing'].get('limit_number_of_loci') else '',
        rank_by = config['processing'].get('rank_by', 'pval')
    message: "working on {wildcards}"
    resources:
        partition="caslake"
    shell:
        """
        Rscript workflow/process/collect_topsnps_results.R --selection_dir {input} --phenotype {wildcards.phenotype} --filtered_sumstats {output.filtered_sumstats} --enformer_loci {output.enformer_loci} --rank_by {params.rank_by} {params.limit_flag}
        """

# ==== Enformer (workflow/enformer/enformer_predict.py): one GPU job per locus, then one merge ====
# each locus: both haplotypes of every individual, middle bins averaged, haplotypes summed

# write the enformer_predict.py config for this phenotype
rule write_enformer_config:
    input:
        loci = rules.collect_top_snps_results.output.enformer_loci
    output:
        os.path.join(ENFORMER_PARAMETERS, f'enformer_config_{runname}_{{phenotype}}.yaml')
    params:
        loci_dir = os.path.join(ENFORMER_PREDICTIONS, runmeta, '{phenotype}')
    run:
        import yaml
        with open(output[0], 'w') as f:
            yaml.dump({
                'loci_file': os.path.abspath(input.loci),
                'individuals': ENFORMER_SETTINGS['individuals'],
                'n_individuals': ENFORMER_SETTINGS['n_individuals'],
                'vcf_pattern': ENFORMER_SETTINGS['vcf_pattern'],
                'fasta_file': config['genome']['fasta'],
                'model_path': config['enformer']['model'],
                'pad_bins': ENFORMER_SETTINGS['pad_bins'],
                'aggregation': ENFORMER_SETTINGS['aggregation'],
                'batch_size': 1,
                'devices': 'auto',
                'output_dir': os.path.abspath(AGGREGATED_PREDICTIONS),
                'output_basename': f'{wildcards.phenotype}.{runmeta}.processed',
                'loci_dir': os.path.abspath(params.loci_dir)
            }, f, sort_keys=False)

# predict one locus on one GPU; the .done marker is written for predicted and for unusable (.invalid) loci
rule predict_enformer_locus:
    input:
        rules.write_enformer_config.output
    output:
        touch(os.path.join(ENFORMER_PREDICTIONS, runmeta, '{phenotype}', '{locus}.done'))
    wildcard_constraints:
        locus = r'chr[0-9XY]+_[0-9]+_[0-9]+'
    params:
        jobname = '{phenotype}_{locus}',
        runmeta = runmeta,
        script = ENFORMER_SETTINGS['script'],
        conda_lib = ENFORMER_SETTINGS['conda_lib']
    resources:
        partition = ENFORMER_SETTINGS['partition'],
        account = ENFORMER_SETTINGS['account'],
        gpu = 1,
        cpu_task = 4,
        mem_cpu = 6,
        time = ENFORMER_SETTINGS['time_per_locus']
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{{locus}}.predict_enformer_locus.tsv")
    shell: "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{params.conda_lib}; nvidia-smi -L || true; python3 {params.script} --config {input} --locus {wildcards.locus} --no-merge"

def enformer_loci_done(wildcards):
    # read the loci chosen by the checkpoint and ask for one .done marker per locus
    loci_file = checkpoints.collect_top_snps_results.get(phenotype = wildcards.phenotype).output.enformer_loci
    with open(loci_file) as f:
        loci = list(dict.fromkeys(line.split()[0] for line in f if line.strip()))
    return expand(os.path.join(ENFORMER_PREDICTIONS, runmeta, wildcards.phenotype, '{locus}.done'), locus = loci)

# merge the per-locus results into the two files prepare_files_for_predictDB reads
rule merge_enformer_predictions:
    input:
        config = rules.write_enformer_config.output,
        done = enformer_loci_done
    output:
        metadata = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.metadata.tsv'),
        matrix = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.matrix.h5.gz')
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        script = ENFORMER_SETTINGS['script']
    resources:
        partition = "caslake",
        mem_cpu = 8,
        cpu_task = 2,
        time = "01:00:00"
    message: "working on {wildcards}"
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.merge_enformer_predictions.tsv")
    shell: "python3 {params.script} --config {input.config} --merge-only"

# prepare for snp-enpact training
rule prepare_files_for_predictDB:
    input: 
        # matrix = rules.process_predictions.output.matrix,
        # metadata = rules.process_predictions.output.metadata
        metadata = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.metadata.tsv'),
        matrix = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.matrix.h5.gz')
    output: 
        enpact_scores = expand(os.path.join(PREDICTDB_DATA, "{{phenotype}}", "{{phenotype}}.{model}.enpact_scores.txt"), model = enpact_models_list),
        annotations = expand(os.path.join(PREDICTDB_DATA, "{{phenotype}}", "{{phenotype}}.{model}.annotation.txt"), model = enpact_models_list)
    params:
        runmeta = runmeta,
        jobname = '{phenotype}',
        blacklist = config['predictdb']['blacklist_regions'],
        output_basename = os.path.join(PREDICTDB_DATA, '{phenotype}', '{phenotype}'),
        enpact_weights = config['enpact_weights'],
        loci_subset = os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.txt') #rules.collect_top_snps_results.output.enformer_loci
    resources:
        partition="caslake",
        mem_cpu=12,
        #mem_mb=24000,
        time = "04:00:00"
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.prepare_files_for_predictDB.tsv")
    shell: 
        """
        python3 workflow/process/enpact_predict.py --matrix {input.matrix} --weights {params.enpact_weights} --metadata {input.metadata} --split --output_basename {params.output_basename} --subset_of_loci {params.loci_subset}
        """

# linearize or train snps
rule generate_lEnpact_models:
    input:
        enpact_scores = os.path.join(PREDICTDB_DATA, "{phenotype}", "{phenotype}.{model}.enpact_scores.txt"),
        annot_file = os.path.join(PREDICTDB_DATA, "{phenotype}", "{phenotype}.{model}.annotation.txt")
    output:
        # temp: snakemake deletes it once format_covariances has written Covariances.varID.txt.gz from it
        covariances_model = temp(os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/predict_db_{phenotype}_filtered.txt.gz')),
        lEnpact_model = os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/predict_db_{phenotype}_filtered.db')
    params:
        jobname = '{phenotype}_{model}',
        runmeta = runmeta,
        output_dir = os.path.abspath(os.path.join(LENPACT_DIR, '{phenotype}', "{model}")),
        annot_file = lambda wildcards, input: os.path.abspath(input.annot_file),
        enpact_scores = lambda wildcards, input: os.path.abspath(input.enpact_scores),
        generate_sbatch = os.path.abspath("workflow/predictdb/generate_snp_predictors.sbatch"),
        reference_genotypes = os.path.abspath(config['predictdb']['reference_genotypes']),
        reference_annotations = os.path.abspath(REFERENCE_ANNOTATIONS),
        nextflow_executable = os.path.abspath(config['predictdb']['nextflow_main_executable'])
    resources:
        partition="caslake",
        time="36:00:00",
        load=5
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{{model}}.generate_lEnpact_models.tsv")
    shell: "cd {params.output_dir} && {params.generate_sbatch} {wildcards.phenotype} {params.output_dir} {params.annot_file} {params.enpact_scores} {params.reference_genotypes} {params.reference_annotations} {params.nextflow_executable}"

# process covariances
rule format_covariances:
    input:
        covariances = rules.generate_lEnpact_models.output.covariances_model
    output:
        formatted_covariances = os.path.join(LENPACT_DIR, "{phenotype}", "{model}", 'models/filtered_db/Covariances.varID.txt.gz')
    params:
        jobname = '{phenotype}_{model}',
        runmeta = runmeta
    resources:
        partition="caslake",
        time="00:30:00"
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{{model}}.format_covariances.tsv")
    shell: "workflow/src/format_covariances.sbatch {input.covariances} {output.formatted_covariances}"

# run summary TFXcan
rule summary_TFXcan:
    input:
        snp_model = rules.generate_lEnpact_models.output.lEnpact_model,
        cov = rules.format_covariances.output.formatted_covariances
    output:
        summary_tfxcan = os.path.join(SUMMARYTFXCAN_DIR, "{phenotype}", "{model}-{phenotype}.enpactScores.spredixcan.csv")
    params:
        jobname = '{phenotype}-{model}',
        runmeta = runmeta,
        gwas_folder = os.path.abspath(os.path.join(PROCESSED_SUMSTATS, '{phenotype}')),
        gwas_pattern = '.*.sumstats.txt.gz',
        executable = config['summaryTFXcan']['summaryXcan_executable'],
        environment = config['summaryTFXcan']['conda_environment']
    resources:
        partition="caslake",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{{model}}.summary_TFXcan.tsv")
    shell: "workflow/process/summary_TFXcan.sbatch {wildcards.phenotype} {input.snp_model} {output.summary_tfxcan} {params.gwas_folder} {params.gwas_pattern} {input.cov} {params.executable} {params.environment}" 

# prepare to collect summary TFXcan
rule prepare_to_collect_summaryTFXcan_results:
    input: 
        lambda wildcards: expand(os.path.join(SUMMARYTFXCAN_DIR, f"{wildcards.phenotype}", f"{{model}}-{wildcards.phenotype}.enpactScores.spredixcan.csv"), model = enpact_models_list)
    output: os.path.join(COLLECTION_DIR, '{phenotype}.summaryTFXcan.paths.txt')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta
    resources:
        partition="caslake",
        time="00:30:00",
        mem_cpu=4,
        cpu_task=4
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.prepare_to_collect_summaryTFXcan_results.tsv")
    run:
        with open(output[0], 'w') as outfile:
            for fname in input:
                outfile.write(fname + "\n")

# collect summary TFXcan into one file
rule collect_summaryTFXcan_results:
    input: rules.prepare_to_collect_summaryTFXcan_results.output[0]
    output:
        summary_tfxcan = os.path.join(SUMMARY_OUTPUT, f'{{phenotype}}.enpactScores.{rundate}.spredixcan.txt')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta
    resources:
        partition="caslake",
        time="00:30:00",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.collect_summaryTFXcan_results.tsv")
    shell: "Rscript workflow/process/collect_summaryTFXcan_results.R --input_files_pattern {input} --phenotype {wildcards.phenotype} --output_file {output.summary_tfxcan}"


rule prepare_summary_matrices:
    input: 
        gwas_input= os.path.join(COLLECTION_DIR, '{phenotype}.filteredGWAS.topSNPs.txt.gz'), #rules.collect_top_snps_results.output.filtered_sumstats,
        tfxcan_input=rules.collect_summaryTFXcan_results.output.summary_tfxcan
    output: 
        zscores_matrices = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.zscores.matrices.rds'),
        zratio_matrix = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.zratios.matrices.rds'),
        subsets_list = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.random_subsets.txt')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta,
        output_basename = f'{{phenotype}}.{runmeta}',
        output_directory=FACTORIZATION_DIR,
        temporary_directory=TEMPORARY_DIR
    resources:
        partition="caslake",
        time="02:00:00",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.prepare_summary_matrices.tsv")
    shell: "Rscript workflow/factorize/prepare_summary_matrices.R --gwas_summary_statistic {input.gwas_input} --tfxcan_summary_statistic {input.tfxcan_input} --output_basename {params.output_basename} --output_directory {params.output_directory} --temporary_directory {params.temporary_directory} --number_resamples 1000 --number_batches 10"

rule repeat_factorization_on_subsets_and_gather:
    input: rules.prepare_summary_matrices.output.zratio_matrix
    output: os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.Iters.rds.gz')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta,
        output_basename=f'{{phenotype}}.{runmeta}',
        output_directory=FACTORIZATION_DIR,
        temporary_directory=TEMPORARY_DIR,
        priorL='ebnm_point_exponential',
        priorF='ebnm_point_exponential',
        greedy_Kmax=40,
        batch_list=rules.prepare_summary_matrices.output.subsets_list,
        splits = os.path.join(TEMPORARY_DIR, f'{{phenotype}}.{runmeta}.random_subsets.{{}}.rds')
    resources:
        partition="caslake",
        time="03:00:00",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.repeat_factorization_on_subsets_and_gather.tsv")
    shell: 
        """
        module load parallel;
        parallel -a {params.batch_list} -j 20 "Rscript workflow/factorize/repeat_flash.R --data {input} --splits {params.splits} --batch {{}} --output_basename {params.temporary_directory}/{params.output_basename} --priorL {params.priorL} --priorF {params.priorF} --greedy_Kmax {params.greedy_Kmax}";
        status=$?
        if [ $status -eq 0 ]; then
            Rscript workflow/factorize/gather_repeats.R --input_basename {params.temporary_directory}/{params.output_basename} --output_basename {params.output_directory}/{params.output_basename} --batch_list {params.batch_list} --priorL {params.priorL} --priorF {params.priorF}
        fi
        """

rule cluster_TFXcan_programs:
    input: rules.repeat_factorization_on_subsets_and_gather.output
    output: 
        best_cluster = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.programs_clara_silhouette.txt.gz'),
        programs_matrix = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.programs_matrix.rds.gz'),
        programs_corrmatrix = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.programs_corrmatrix.rds.gz'),
        programs_clusters = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.program_clusters.txt.gz'),
        loci_clusters = os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.loci_assignments.txt.gz')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta,
        output_basename = f'{{phenotype}}.{runmeta}',
        output_directory=FACTORIZATION_DIR
    resources:
        partition="caslake",
        time="02:00:00",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.cluster_TFXcan_programs.tsv")
    shell: "Rscript workflow/factorize/cluster_programs.R --flash_results {input} --output_basename {params.output_directory}/{params.output_basename}"