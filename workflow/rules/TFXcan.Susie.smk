rule process_summary_statistics:
    input: os.path.join(INPUT_SUMSTATS, '{phenotype}.liftover.logistic.assoc.tsv.gz')
    output: directory(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'))
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        jobname = '{phenotype}',
        diag_file = os.path.join(DATA_DIR, 'diagnostics', f'{{phenotype}}.gwas_diagnostics.summary')
    message: "working on {wildcards}" 
    resources:
        mem_mb = 10000
    shell:
        """
        {params.rscript} workflow/src/process_summary_statistics.R --summary_stats_file {input} --output_folder {output} --phenotype {wildcards.phenotype} --diagnostics_file {params.diag_file}
        """

rule run_susie_on_summary_statistics: 
    input: 
        rules.process_summary_statistics.output,
        collect_processed_summary_statistics
    output:
        directory(os.path.join(FINEMAPPING_DIR, '{phenotype}'))
    params:
        runmeta = runmeta,
        rscript = config['rscript'],
        jobname = '{phenotype}',
        input_sumstats = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'chr{{}}.sumstats.txt.gz'),
        ld_blocks = config['finemapping']['LD_blocks'],
        patterns = config['finemapping']['genotypes_dosages_pattern'],
        chroms = collect_chromosomes,
        diag_file = os.path.join(DATA_DIR, 'diagnostics', f'{{phenotype}}.chr{{}}.susie_diagnostics.summary')
    message: "working on {wildcards}"
    resources:
        partition="beagle3",
        mem_cpu=12
    shell:
        """
        module load parallel;
        printf "%s\\n" {params.chroms} | parallel -j 12 "{params.rscript} workflow/src/run_susie_on_summary_statistics.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {output} --phenotype {wildcards.phenotype} --pip_threshold 0.5 --diagnostics_file {params.diag_file} --genotypes_dosages_pattern {params.patterns}"
        """

rule collect_finemapping_results:
    input: rules.run_susie_on_summary_statistics.output
    output:
        finemapped_sumstats = os.path.join(COLLECTION_DIR, '{phenotype}.filteredGWAS.txt.gz'),
        enformer_loci = os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.txt')
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        jobname = '{phenotype}'
    message: "working on {wildcards}"
    resources:
        partition="beagle3"
    shell:
        """
        {params.rscript} workflow/src/collect_finemapping_results.R --finemapping_dir {input} --phenotype {wildcards.phenotype} --filtered_sumstats {output.finemapped_sumstats} --enformer_loci {output.enformer_loci}
        """

rule create_enformer_configuration:
    input: rules.collect_finemapping_results.output.enformer_loci
    output: os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.json')
    message: "working on {wildcards}"
    resources:
        partition="beagle3"
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        bdirectives = config['enformer']['base_directives'],
        dset = runname,
        model = config['enformer']['model'],
        fasta_file = config['genome']['fasta'],
        pdir = ENFORMER_PREDICTIONS,
        ddate = rundate,
        jobname = '{phenotype}',
        personalized_predictions = config['personalized_predictions'],
        personalized_directives = config['enformer']['personalized_directives']
    run:  
        if params.personalized_predictions == True:
            shell("{params.rscript} workflow/src/create_enformer_config.R --runname {params.dset} --phenotype {wildcards.phenotype} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate} --personalized_parameters_file {params.personalized_directives}")
        elif params.personalized_predictions == False:
            shell("{params.rscript} workflow/src/create_enformer_config.R --runname {params.dset} --phenotype {wildcards.phenotype} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate}")

rule predict_with_enformer:
    input:
        rules.create_enformer_configuration.output
    output:
        os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.json')
    resources:
        partition="beagle3",
        time = "12:00:00",
        gpu=4,
        mem_cpu=8,
        cpu_task=8
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        enformer_predict_script = config['enformer']['predict']
    message: 
        "working on {params.jobname}"
    shell:
        """
        python3 {params.enformer_predict_script} --parameters {input}
        """

rule aggregate_predictions:
    input:
        rules.predict_with_enformer.output
    output:
        #touch(os.path.join(CHECKPOINTS_DIR, '{phenotype}.checkpoint')),
        directory(os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}'))
    message: 
        "working on {wildcards}"
    resources:
        partition="beagle3",
        mem_cpu=8,
        cpu_task=8,
        mem_mb=24000
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        aggregation_script = config['enformer']['aggregate'],
        aggtype = config['enformer']['aggtype'],
        hpc = "caslake",
        parsl_executor = "local",
        delete_enformer_outputs = config["delete_enformer_outputs"],
        aggregation_config = rules.predict_with_enformer.output,
        output_dir = os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}')
    run:
        if params.delete_enformer_outputs == True:
            shell("mkdir -p {params.output_dir} && python3 {params.aggregation_script} --metadata_file {params.aggregation_config} --agg_types {params.aggtype} --output_directory {params.output_dir} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
        elif params.delete_enformer_outputs == False: # don't delete the outputs
            shell("mkdir -p {params.output_dir} && python3 {params.aggregation_script} --metadata_file {params.aggregation_config} --agg_types {params.aggtype} --output_directory {params.output_dir} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")


rule calculate_enpact_scores:
    input: os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}') #rules.aggregate_predictions.output 
    output: directory(os.path.join(ENPACT_PREDICTIONS, '{phenotype}'))
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        individuals = lambda wildcards: collect_aggregated_individuals(),
        input_file = lambda wildcards, input: os.path.join(AGGREGATED_PREDICTIONS, wildcards.phenotype, f'{{}}_{config["enformer"]["aggtype"]}_{wildcards.phenotype}.csv.gz'),
        output_file = lambda wildcards, output: os.path.join(ENPACT_PREDICTIONS, wildcards.phenotype, f'{{}}.{wildcards.phenotype}.{config["enformer"]["aggtype"]}.{rundate}.csv.gz'),
        models_directory = config['enpact_models']['directory'],
        models_metadata = config['enpact_models']['metadata'],
        jobname = '{phenotype}'
    message: 
        "working on {params.jobname}"
    resources:
        partition="beagle3",
        mem_cpu=12
    shell:
        """
        module load parallel;
        mkdir -p {output};
        printf "%s\\n" {params.individuals} | parallel -j 10 '{params.rscript} workflow/src/calculate_enpact_scores.R --input_file {params.input_file} --output_file {params.output_file} --enpact_models_directory {params.models_directory} --enpact_models_metadata {params.models_metadata}'
        """

rule create_enpact_scores_database:
    input: rules.calculate_enpact_scores.output
    output: 
        f1 = os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.array.rds.gz"),
        f2 = os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.txt.gz")
    message: "working on {wildcards.phenotype}"
    resources:
        partition="beagle3"
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        jobname = '{phenotype}',
        input_pattern = os.path.join(ENPACT_PREDICTIONS, '{phenotype}', f'{{}}.{{phenotype}}.{config["enformer"]["aggtype"]}.{rundate}.csv.gz'),
        individuals = lambda wildcards: ','.join(collect_enpact_individuals())
    shell:
        """
            {params.rscript} workflow/src/create_enpact_scores_database.R --input_files {params.input_pattern} --output_file {output.f2} --output_db {output.f1} --individuals {params.individuals}
        """

rule prepare_files_for_predictDB:
    input: rules.create_enpact_scores_database.output.f2
    output: 
        enpact_scores = os.path.join(PREDICTDB_DATA, '{phenotype}.enpact_scores.txt'),
        annot_file = os.path.join(PREDICTDB_DATA, '{phenotype}.tf_tissue_annot.txt')
    params:
        rscript = config['rscript'],
        runmeta = runmeta,
        jobname = '{phenotype}',
        output_dir = os.path.join(PREDICTDB_DATA, '{phenotype}'),
        filtered_gwas = rules.collect_finemapping_results.output.finemapped_sumstats,
        blacklist = config['predictdb']['blacklist_regions']
    resources:
        partition="caslake"
    shell: "{params.rscript} workflow/src/prepare_files_for_predictDB.R --enpact_scores_file {input} --formatted_escores_file {output.enpact_scores} --formatted_annot_file {output.annot_file} --filtered_GWAS_file {params.filtered_gwas} --blacklist {params.blacklist}"

rule generate_lEnpact_models:
    input:
        enpact_scores = rules.prepare_files_for_predictDB.output.enpact_scores,
        annot_file = rules.prepare_files_for_predictDB.output.annot_file
    output:
        lEnpact_model = os.path.join('output', 'lEnpact', '{phenotype}', 'models/filtered_db/predict_db_{phenotype}_filtered.db'),
        covariances_model = os.path.join('output', 'lEnpact', '{phenotype}', 'models/database/predict_db_{phenotype}.txt.gz')
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        output_dir = os.path.abspath(os.path.join('output', 'lEnpact')),
        annot_file = os.path.abspath(rules.prepare_files_for_predictDB.output.annot_file),
        enpact_scores = os.path.abspath(rules.prepare_files_for_predictDB.output.enpact_scores)
    resources:
        partition="caslake",
        time="36:00:00"
    shell: "workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch {wildcards.phenotype} {params.output_dir} {params.annot_file} {params.enpact_scores}"

            # gzip -d -k {input.covariances};
            # sed -e 's/:/_/g' {output.f1} > {output.f2}
rule format_covariances:
    input:
        covariances = rules.generate_lEnpact_models.output.covariances_model
    output:
        f1 = os.path.join('output', 'lEnpact', '{phenotype}/models/database/predict_db_{phenotype}.txt'),
        f2 = os.path.join('output', 'lEnpact', '{phenotype}/models/database/Covariances.varID.txt')
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        # c1 = os.path.abspath(os.path.join('output', 'lEnpact', '{phenotype}/models/database/predict_db_{phenotype}.txt.gz')),
        # f1 = os.path.abspath(os.path.join('output', 'lEnpact', '{phenotype}/models/database/predict_db_{phenotype}.txt')),
        # f2 = os.path.abspath(os.path.join('output', 'lEnpact', '{phenotype}/models/database/Covariances.varID.txt'))
    resources:
        partition="caslake",
        time="01:00:00"
    shell:
        """
            touch {input.covariances};
            touch {output.f1};
            touch {output.f2};
        """

rule summary_TFXcan:
    input:
        model=rules.generate_lEnpact_models.output.lEnpact_model,
        cov=rules.format_covariances.output.f2
    output:
        summary_tfxcan = os.path.join('output', 'summary', '{phenotype}/{phenotype}.enpactScores.spredixcan.csv')
    params:
        jobname = '{phenotype}',
        runmeta = runmeta,
        covariances = os.path.abspath(os.path.join('output', 'lEnpact', '{phenotype}/models/database/Covariances.varID.txt')),
        gwas_folder = os.path.abspath(os.path.join(PROCESSED_SUMSTATS, '{phenotype}')),
        gwas_pattern = '.*.sumstats.txt.gz',
        inp = os.path.abspath(rules.generate_lEnpact_models.output.lEnpact_model),
        outp = os.path.abspath(os.path.join('output', 'summary', '{phenotype}/{phenotype}.enpactScores.spredixcan.csv'))
    resources:
        partition="bigmem",
        time="04:00:00"
    shell: "workflow/src/summary_TFXcan.sbatch {wildcards.phenotype} {input.model} {output.summary_tfxcan} {params.gwas_folder} {params.gwas_pattern} {input.cov}" #"workflow/src/summary_TFXcan.sbatch {wildcards.phenotype} {params.inp} {params.outp} {params.gwas_folder} {params.gwas_pattern} {params.covariances}"