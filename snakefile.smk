# Description: Given a TF and tissue (or context) , this pipeline trains logistic elastic net models of that TF binding activity in that tissue
# Author: Temi
# Date: Wed Mar 29 2023
# Usage: --> see README

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards, Wildcards, expand
import numpy as np
import collections

sys.path.append('workflow/src')

#import module
print_progress = False
# configfile: '/project2/haky/temi/projects/TFXcan-snakemake/config/parameters.yaml'
# from snakemake.utils import validate
# validate(config, "config.schema.yaml")

# directories
DATA_DIR = 'data' 
INPUT_SUMSTATS = os.path.join(DATA_DIR, 'sumstats')
PROCESSED_SUMSTATS = os.path.join(DATA_DIR, 'processed_sumstats')
FINEMAPPING_DIR = os.path.join(DATA_DIR, 'finemapping')
COLLECTION_DIR = os.path.join(DATA_DIR, 'collection')
ENFORMER_PARAMETERS = os.path.join(DATA_DIR, 'enformer_parameters')
ENFORMER_PREDICTIONS = os.path.join(config['scratch_dir'], f"{config['runname']}_{config['date']}", 'predictions_folder') if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'enformer_predictions')
AGGREGATED_PREDICTIONS = os.path.join(DATA_DIR, 'aggregated_predictions')
ENPACT_PREDICTIONS = os.path.join(DATA_DIR, 'enpact_predictions')

def read_metadata(mtdt_file):
    dd = pd.read_csv(mtdt_file)
    #return({'phenotypes': dd.phenotype.tolist(), 'sumstats': dd.sumstat.tolist()})
    return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))
    #return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))

run_list = read_metadata(config["metadata"])
print(run_list)


# printf '%s\\n' {{1..22}} X | parallel -j 23 "plink --r --ld-window-r2 0 --bfile {params.bfile_pattern} --ld-snp-list {params.snplist} --out {params.outputLD}"
#phe, chrom_ = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, '{phenotype}', 'chr{chrom}.sumstats.txt.gz'))
# print(phe)
# print(chrom_)

def collect_processed_summary_statistics(wildcards):
    #checkpoint_output = checkpoints.process_summary_statistics.get(**wildcards).output[0]
    chromosomes = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, '{phenotype}', 'chr{chrom}.sumstats.txt.gz')).chrom #glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom
    output = expand(os.path.join(PROCESSED_SUMSTATS, "{phenotype}", 'chr{chrom}.sumstats.txt.gz'), chrom=chromosomes, phenotype = run_list.keys())
    # print(output)
    return output

def collect_chromosomes(wildcards):
    #checkpoint_output = checkpoints.process_summary_statistics.get(**wildcards).output[0]
    chromosomes = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, "{phenotype}", 'chr{chrom}.sumstats.txt.gz')).chrom
    return(chromosomes)

# print(glob_wildcards(os.path.join(PROCESSED_SUMSTATS, '{phenotype}', 'chr{chrom}.sumstats.txt.gz')).chrom)

rule all:
    input:
        #expand(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(FINEMAPPING_DIR, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(COLLECTION_DIR, '{phenotype}', f'{{phenotype}}.filteredGWAS.txt.gz'), phenotype = run_list.keys()),
        expand(os.path.join(COLLECTION_DIR, '{phenotype}', f'{{phenotype}}.EnformerLoci.txt'), phenotype = run_list.keys()),
        expand(os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{config["runname"]}_{{phenotype}}.json'), phenotype = run_list.keys()),
        expand(os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{config["runname"]}_{{phenotype}}.json'), phenotype = run_list.keys()),
        expand(os.path.join(AGGREGATED_PREDICTIONS, f'{config["runname"]}_{config["enformer"]["aggtype"]}_{{phenotype}}.csv.gz'), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_PREDICTIONS, f'{config["runname"]}_{{phenotype}}_{config["date"]}.csv.gz'), phenotype = run_list.keys())
        #expand('{phenotype}_report.txt', phenotype = run_list.keys())
        #input_file=expand(os.path.join(INPUT_SUMSTATS, '{sumstat}'), zip, phenotype = run_list.keys(), sumstat = run_list.values()),
        
checkpoint process_summary_statistics:
    input: lambda wildcards: os.path.join(INPUT_SUMSTATS, '{phenotype}')
    output: directory(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'))
    params:
        rscript = config['rscript'],
        jobname = '{phenotype}',
        #input_file = os.path.join(INPUT_SUMSTATS, run_list[phenotype]),
    message: "working on {wildcards}" 
    resources:
        mem_mb = 10000
    shell:
        """
        {params.rscript} workflow/src/process_summary_statistics.R --summary_stats_file {input} --output_folder {output} --phenotype {wildcards.phenotype}
        """

checkpoint run_susie_on_summary_statistics: 
    input: 
        collect_processed_summary_statistics
    output:
        directory(os.path.join(FINEMAPPING_DIR, '{phenotype}'))
    params:
        rscript = config['rscript'],
        jobname = '{phenotype}',
        input_sumstats = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'chr{{}}.sumstats.txt.gz'),
        ld_blocks = config['finemapping']['LD_blocks'],
        #input_files = aggregate_for_susie,
        chroms = collect_chromosomes
    message: "working on {wildcards}"
    resources:
        partition="caslake"
    shell:
        """
        module load parallel;
        printf "%s\\n" {params.chroms} | parallel -j 5 "{params.rscript} workflow/src/run_susie_on_summary_statistics.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {output} --phenotype {wildcards.phenotype}"
        """

rule collect_finemapping_results:
    input:
        #lambda wildcards: checkpoints.run_susie_on_summary_statistics.get(**wildcards).output
        os.path.join(FINEMAPPING_DIR, '{phenotype}')
        #os.path.join(FINEMAPPING_DIR, '{phenotype}', f'{{phenotype}}.chr{{chrom}}.filteredGWAS.txt.gz')
    output:
        finemapped_sumstats = os.path.join(COLLECTION_DIR, '{phenotype}', f'{{phenotype}}.filteredGWAS.txt.gz'),
        enformer_loci = os.path.join(COLLECTION_DIR, '{phenotype}', f'{{phenotype}}.EnformerLoci.txt')
    params:
        rscript = config['rscript'],
        jobname = '{phenotype}',
        #input_dir = os.path.join(FINEMAPPING_DIR, '{phenotype}')
    message: "working on {wildcards}"
    resources:
        partition="caslake"
    shell:
        """
        {params.rscript} workflow/src/collect_finemapping_results.R --finemapping_dir {input} --phenotype {wildcards.phenotype} --filtered_sumstats {output.finemapped_sumstats} --enformer_loci {output.enformer_loci}
        """


rule create_enformer_configuration:
    input: rules.collect_finemapping_results.output.enformer_loci
    output: os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{config["runname"]}_{{phenotype}}.json')
    message: "working on {wildcards}"
    resources:
        partition="beagle3"
    params:
        rscript = config['rscript'],
        bdirectives = config['enformer']['base_directives'],
        dset = config['runname'],
        model = config['enformer']['model'],
        fasta_file = config['genome']['fasta'],
        pdir = DATA_DIR,
        ddate = config['date'],
        jobname = '{phenotype}'
    shell:
        """
            {params.rscript} workflow/src/create_enformer_config.R --runname {params.dset} --phenotype {wildcards.phenotype} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate}
        """

rule predict_with_enformer:
    input:
        rules.create_enformer_configuration.output
    output:
        os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{config["runname"]}_{{phenotype}}.json')
    resources:
        partition="beagle3",
        time="04:00:00",
        gpu=4,
        mem_cpu=8,
        cpu_task=8
    params:
        jobname = '{phenotype}',
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
        os.path.join(AGGREGATED_PREDICTIONS, f'{config["runname"]}_{config["enformer"]["aggtype"]}_{{phenotype}}.csv.gz')
    message: 
        "working on {wildcards}"
    resources:
        partition="beagle3",
        mem_cpu=8,
        cpu_task=8,
        mem_mb=24000
    params:
        jobname = '{phenotype}',
        aggregation_script = config['enformer']['aggregate'],
        aggtype = config['enformer']['aggtype'],
        output_folder = AGGREGATED_PREDICTIONS,
        hpc = "beagle3",
        parsl_executor = "local",
        delete_enformer_outputs = config["delete_enformer_outputs"]
    run:
        if params.delete_enformer_outputs == True:
            shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
        elif params.delete_enformer_outputs == False: # don't delete the outputs
            shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")

rule calculate_enpact_score:
    input:
        rules.aggregate_predictions.output
    output:
        os.path.join(ENPACT_PREDICTIONS, f'{config["runname"]}_{{phenotype}}_{config["date"]}.csv.gz')
    params:
        models_directory = config['enpact_models']['directory'],
        models_filters_date = config['enpact_models']['filters']['date'],
        models_filters_type = config['enpact_models']['filters']['type'],
        jobname = '{phenotype}',
        rscript = config['rscript']
    message: 
        "working on {params.jobname}"
    shell:
        """
            {params.rscript} workflow/src/calculate_enpact_score.R --input_file {input} --output_file {output} --enpact_models_directory {params.models_directory} --filters_date {params.models_filters_date} --filters_type {params.models_filters_type}
        """






#print(checkpoints.process_summary_statistics.get(**wildcards).output)





# checkpoint report:
#     input:
#         collect_processed_summary_statistics
#     output:
#         '{phenotype}_report.txt'
#     params:
#         jobname = '{phenotype}',
#         out = lambda wildcards: f'{wildcards.phenotype}_report.txt',
#         chroms = collect_chromosomes
#     message: "working on {wildcards}"
#     shell:
#         "echo {input} {params.chroms} > {output}"








# rule run_susie_on_summary_statistics: 
#     input: 
#         aggregate_for_susie
#     output:
#         directory(os.path.join(FINEMAPPING_DIR, '{phenotype}'))
#     params:
#         rscript = config['rscript'],
#         jobname = '{phenotype}',
#         input_sumstats = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'chr{{}}.sumstats.txt.gz'),
#         ld_blocks = config['finemapping']['LD_blocks'],
#         input_files = aggregate_for_susie,
#         chroms = get_chromosomes
#     message: "working on {wildcards}"
#     resources:
#         partition="caslake"
#     shell:
#         """
#         module load parallel;
#         printf "%s\\n" {params.chroms} | parallel -j 5 "{params.rscript} workflow/src/run_susie_on_summary_statistics.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {output} --phenotype {wildcards.phenotype}"
#         """





# # module load parallel;
# #         printf "%s\\n" {{1..22}} X | parallel -j 1 "{params.rscript} workflow/src/run_susie_on_summary_statistics.R --chromosome {{}} --sumstats {params.input_sumstats} --LDBlocks_info {params.ld_blocks} --output_folder {params.output_folder} --phenotype {wildcards.phenotype}"









# rule calculate_LD: # uses plink
#     input:
#         rules.process_summary_statistics.output.LD_SNPs_folder
#     output:
#         directory(os.path.join(LD_DIR, '{phenotype}'))
#     params:
#         rscript = config['rscript'],
#         jobname = '{phenotype}',
#         bfile_pattern = os.path.join(config['plink']['folder'], config['plink']['basename']),
#         snplist = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'{wildcards.phenotype}.chr{{}}.SNPsForLD'),
#         outputLD = lambda wildcards: os.path.join(LD_DIR, wildcards.phenotype, f'{wildcards.phenotype}.chr{{}}')
#     message: "working on {wildcards}"
#     resources:
#         partition="caslake"
#     shell:
#         """
#         module load parallel;
#         dname=`dirname {params.outputLD}`;
#         mkdir -p ${{dname}};
#         {params.rscript} workflow/src/calculate_LD.R --bfile_pattern {params.bfile_pattern} --snplist_pattern {params.snplist} --output_pattern {params.outputLD}
#         """

# rule prepare_LD_matrix:
#     input: 
#         LD_folder = rules.calculate_LD.output,
#         SNPList_folder = rules.process_summary_statistics.output.LD_SNPs_folder
#     output: rules.calculate_LD.output
#     params:
#         rscript = config['rscript'],
#         jobname = '{phenotype}',
#         reference_genome_txt_file = config['finemapping']['genotypes'],
#         bfile_pattern = os.path.join(config['plink']['folder'], config['plink']['basename']),
#         snplist = lambda wildcards: os.path.join(PROCESSED_SUMSTATS, wildcards.phenotype, f'{wildcards.phenotype}.chr{{}}.SNPsForLD'),
#         outputLD = lambda wildcards: os.path.join(LD_DIR, wildcards.phenotype, f'{wildcards.phenotype}.chr{{}}')

# rule count:
#     input: rules.calculate_LD.output
#     output: os.path.join(DATA_DIR, '{phenotype}.count.txt')
#     shell:
#         """
#         ls {input} | wc -l > '{output}'
#         """



#print(checkpoints.process_summary_statistics.get('phenotype'))

# lambda wildcards: print(wildcards)

# rule calculate_LD: # uses plink
#     input:
#         checkpoints.process_summary_statistics.get()
#     output:
#         os.path.join(LD_DIR, '{phenotype}', f'{{phenotype}}.chr{{}}')
#         #touch(os.path.join(CORRELATION_MATRIX, 'correlation_matrix.rds'))
#     params:
#         rscript = config['rscript'],
#         jobname = 'correlation_matrix',
#         reference_genome_txt_file = config['finemapping']['genotypes'],
#         bfile_pattern = os.path.join(config['plink']['folder'], config['plink']['basename'])
#     message: "working on {wildcards}"
#     resources:
#         partition="caslake"
#     shell:
#         """
#         printf "%s\n" {1..22} X | parallel -j 23 "plink --r square gz --ld-window-r2 0 --bfile {params.bfile_pattern} --ld-snp-list os.path.join({{input}}, {{wildcard.phenotype}}.chr{}.SNPsForLD --out {output}
#         """





# rule check_genome_build:
#     output: 
#         touch(os.path.join(DATA_DIR, 'genome_build.txt'))
#     params:
#         rscript = config['rscript'],
#         jobname = 'genome_build'
#     message: "working on {wildcards}"
#     resources:
#         mem_mb = 10000
#     shell:
#         """
#         {params.rscript} workflow/src/check_genome_build.R --output_file {output} --reference_genome_txt_file {config['finemapping']['genotypes']}
#         """

# rule run_susie_on_summary_statistics:
#     input:
#         rules.process_summary_statistics.output.processed_summary_stats
#     output:
#         os.path.join(DATA_DIR, 'susie_output', '{phenotype}.susie.rds')
#     params:
#         rscript = config['rscript'],
#         jobname = '{phenotype}'
#     message: "working on {wildcards}"
#     resources:
#         mem_mb = 10000
#     shell:
#         """
#         {params.rscript} workflow/src/run_susie_on_summary_statistics.R --processed_sumstats_file {input} --output_file {output} --correlation_matrix {os.path.join(DATA_DIR, 'correlation_matrix.rds')} --zscores_column 'zscore' 
#         """




    



# checks consistency of summary statistics
# rule process_summary_statistics: 
#     input: os.path.join(PROCESSED_SUMSTATS, '{phenotype}')
#     output:
#         lambda wildcards: dynamic(os.path.join(PROCESSED_SUMSTATS, '{wildcards.phenotype}', f'{{wildcards.phenotype}}.chr{{chrom}}.sumstats.txt.gz'))
#     params:
#         rscript = config['rscript'],
#         jobname = '{phenotype}',
#         outdir = os.path.join(PROCESSED_SUMSTATS, '{phenotype}')
#     message: "working on {wildcards}" 
#     resources:
#         mem_mb = 10000
#     shell:
#         """
#         {params.rscript} workflow/src/process_summary_statistics.R --summary_stats_file {input} --output_folder {params.outdir} --phenotype {wildcards.phenotype}
#         """