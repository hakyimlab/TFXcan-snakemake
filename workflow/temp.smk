
TEMPORARY_DIR = os.path.join(config['scratch_dir'], 'factorize.temp') if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'factorize.temp')
FACTORIZATION_DIR = os.path.join(DATA_DIR, 'factorize')


expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.zscores.matrices.rds'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.zratios.matrices.rds'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.random_subsets.txt'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.Iters.rds.gz'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.programs_clara_silhouette.txt.gz'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.programs_matrix.rds.gz'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.programs_corrmatrix.rds.gz'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}program_clusters.txt.gz'), phenotype = run_list.keys(), runmeta = [runmeta]),
expand(os.path.join(FACTORIZATION_DIR, '{phenotype}.{runmeta}.loci_assignments.txt.gz'), phenotype = run_list.keys(), runmeta = [runmeta])


# collect summary TFXcan into one file
rule collect_summaryTFXcan_results:
    input: lambda wildcards: checkpoints.prepare_to_collect_summaryTFXcan_results.get(phenotype = wildcards.phenotype).output[0]
    output: summary_tfxcan = os.path.join(SUMMARY_OUTPUT, f'{{phenotype}}.enpactScores.{rundate}.spredixcan.txt')
    message: "working on {wildcards}"
    params:
        jobname = runmeta,
        runmeta = runmeta
    resources:
        partition="caslake",
        time="01:00:00",
        mem_cpu=4,
        cpu_task=8
    benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.collect_summaryTFXcan_results.tsv")
    shell: "Rscript workflow/process/collect_summaryTFXcan_results.R --input_files_pattern {input} --phenotype {wildcards.phenotype} --output_file {output.summary_tfxcan}"

rule prepare_summary_matrices:
    input: 
        gwas_input=rules.collect_top_snps_results.output.filtered_sumstats,
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
        parallel -a {params.batch_list} -j 20 "Rscript workflow/factorize/repeat_flash.R --data {input} --splits {params.splits} --batch {{}} --output_basename {params.output_directory}/{params.output_basename} --priorL {params.priorL} --priorF {params.priorF} --greedy_Kmax {params.greedy_Kmax}";
        status=$?
        if [ $status -eq 0 ]; then
            Rscript workflow/factorize/gather_repeats.R --input_basename {params.temporary_directory}/{params.output_basename} --output_basename {params.output_directory}/{params.output_basename} --batch_list {params.batch_list} --priorL {params.priorL} --priorF {params.priorF}
        fi
        """

rule cluster_TFXcan_programs:
    input: rules.compile_mini_factorizations.output
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
















# rule compile_mini_factorizations:
#     input: rules.prepare_summary_matrices.output.zratio_matrix
#     output: os.path.join(FACTORIZATION_DIR, f'{{phenotype}}.{runmeta}.Iters.rds.gz')
#     message: "working on {wildcards}"
#     params:
#         jobname = runmeta,
#         runmeta = runmeta,
#         output_basename = f'{{phenotype}}.{runmeta}',
#         output_directory=FACTORIZATION_DIR,
#         temporary_directory=TEMPORARY_DIR,
#         priorL='ebnm_point_exponential',
#         priorF='ebnm_point_exponential',
#         batch_list=rules.prepare_summary_matrices.output.subsets_list
#     resources:
#         partition="caslake",
#         time="02:00:00",
#         mem_cpu=4,
#         cpu_task=8
#     benchmark: os.path.join(f"{BENCHMARK_DIR}/{{phenotype}}.{rundate}.compile_mini_factorizations.tsv")
#     shell: "Rscript workflow/factorize/gather_repeats.R --input_basename {params.temporary_directory}/{params.output_basename} --output_basename {params.output_directory}/{params.output_basename} --batch_list {params.batch_list} --priorL {params.priorL} --priorF {params.priorF}"