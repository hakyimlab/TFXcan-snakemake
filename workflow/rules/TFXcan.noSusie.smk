






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
        select_n_snps = config['processing']['select_n_snps']
    message: "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=12
    shell:
        """
        module load parallel;
        printf "%s\\n" {params.chroms} | parallel -j 12 "Rscript workflow/process/select_top_snps.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {output} --phenotype {wildcards.phenotype} --diagnostics_file {params.diag_file} --selection_method {params.selection_method} --select_n_snps {params.select_n_snps}"
        """

# collect the result above into one data
rule collect_top_snps_results:
    input: lambda wildcards: checkpoints.select_top_snps.get(**wildcards).output[0]
    output:
        filtered_sumstats = os.path.join(COLLECTION_DIR, '{phenotype}.filteredGWAS.topSNPs.txt.gz'),
        enformer_loci = os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.txt')
    params:
        runmeta = runmeta,
        jobname = '{phenotype}'
    message: "working on {wildcards}"
    resources:
        partition="caslake"
    shell:
        """
        Rscript workflow/process/collect_topsnps_results.R --selection_dir {input} --phenotype {wildcards.phenotype} --filtered_sumstats {output.filtered_sumstats} --enformer_loci {output.enformer_loci}
        """

# create configuration file for enformer to run
rule create_enformer_configuration:
    input: rules.collect_top_snps_results.output.enformer_loci
    output: os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.yaml')
    message: "working on {wildcards}"
    resources:
        partition="caslake"
    params:
        runmeta = runmeta,
        bdirectives = config['enformer']['base_directives'],
        dset = runname,
        model = config['enformer']['model'],
        fasta_file = config['genome']['fasta'],
        pdir = ENFORMER_PREDICTIONS,
        ddate = rundate,
        jobname = '{phenotype}',
        personalized_predictions = config['personalized_predictions'],
        personalized_directives = config['enformer']['personalized_directives'],
        cp_aggregation = os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.yaml')
    run:  
        if params.personalized_predictions == True:
            shell("Rscript workflow/process/create_enformer_config.R --runname {params.dset} --phenotype {wildcards.phenotype} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate} --personalized_parameters_file {params.personalized_directives} --copy_aggregation_config {params.cp_aggregation}")
        elif params.personalized_predictions == False:
            shell("Rscript workflow/process/create_enformer_config.R --runname {params.dset} --phenotype {wildcards.phenotype} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate} --copy_aggregation_config {params.cp_aggregation}")

# predict with enformer
rule predict_with_enformer:
    input:
        rules.create_enformer_configuration.output
    output:
        os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.yaml')
    resources:
        partition="caslake",
        time = "36:00:00",
        mem_cpu=4,
        cpu_task=4
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        enformer_predict_script = config['enformer']['predict']
    message: 
        "working on {params.jobname}"
    shell:
        """
        sbatch workflow/enformer/enformer_predict.sbatch {params.enformer_predict_script} {input}
        """

# collect the predictions into one database file
rule aggregate_predictions:
    input:
        rules.predict_with_enformer.output
    output:
        os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.h5')
    message: 
        "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=12,
        #mem_mb=24000,
        time = "04:00:00"
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        delete_enformer_outputs = config["delete_enformer_outputs"],
        aggregation_config = rules.predict_with_enformer.output,
        output_dir = AGGREGATED_PREDICTIONS,
        output_filename = os.path.join(f'{{phenotype}}.{runmeta}.h5')
    shell:
        """
        python3 workflow/enformer/enformer_merge.py --config {input} --output_directory {params.output_dir} --output_filename {params.output_filename}
        """

# process the predictions (sum haplotypes)
rule process_predictions:
    input:
        rules.aggregate_predictions.output
    output:
        metadata = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.metadata.tsv'),
        matrix = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.matrix.h5.gz')
    message: 
        "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=12,
        #mem_mb=24000,
        time = "04:00:00"
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        basename = os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed')
    shell:
        """
        python3 workflow/enformer/enformer_process.py --merged_h5_file {input} --process_by_haplotype --process_function 'sum' --output_basename {params.basename}
        """

# prepare for snp-enpact training
checkpoint prepare_files_for_predictDB:
    input: 
        matrix = rules.process_predictions.output.matrix,
        metadata = rules.process_predictions.output.metadata
    output: 
        enpact_scores = expand(os.path.join(PREDICTDB_DATA, "{{phenotype}}", "{{phenotype}}.{model}.enpact_scores.txt"), model = enpact_models_list),
        annotations = expand(os.path.join(PREDICTDB_DATA, "{{phenotype}}", "{{phenotype}}.{model}.annotation.txt"), model = enpact_models_list)
    params:
        runmeta = runmeta,
        jobname = '{phenotype}',
        blacklist = config['predictdb']['blacklist_regions'],
        output_basename = os.path.join(PREDICTDB_DATA, '{phenotype}', '{phenotype}'),
        enpact_weights = config['enpact_weights'],
        loci_subset = rules.collect_top_snps_results.output.enformer_loci
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
        covariances_model = os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/predict_db_{phenotype}_filtered.txt.gz'),
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
checkpoint summary_TFXcan:
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
checkpoint prepare_to_collect_summaryTFXcan_results:
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
    input: lambda wildcards: checkpoints.prepare_to_collect_summaryTFXcan_results.get(phenotype = wildcards.phenotype).output[0]
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