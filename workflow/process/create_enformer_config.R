# Author: Temi
# Date: Thursday July 27 2023
# Description: script used to create enformer predict parameters file
# Usage: Rscript create_training_sets.R [OPTIONS]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--phenotype", help='A GWAS phenotype e.g. alzheimers'),
    make_option("--runname", help="The name of the dataset"),
    make_option("--base_directives", help='A yaml file containing directives for enformer; will be used to create a json file for the enformer predict pipeline'),
    make_option("--project_directory", help='A project directory for enformer predict'),
	make_option("--predictors_file", help='predictor file containing the intervals to predict on'),
    make_option("--model", help='enformer model'),
	make_option("--fasta_file", help='fasta file, typically hg38'),
    make_option("--date", help='fasta file, typically hg38'),
    make_option("--parameters_file", help='the json file that will be created'),
    make_option("--personalized_parameters_file", default=NULL, help='the json file that will be created'),
    make_option("--copy_aggregation_config", default=NULL, help='the yaml file that will be created')
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# opt <- list()
# opt$runname <- 'Asthma_GWAS' 
# opt$phenotype <- 'asthma_children' 
# opt$base_directives <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/config/enformer_base.yaml' 
# opt$project_directory <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data'
# opt$parameters_file <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/misc/enformer_parameters_Asthma_GWAS_asthma_children.yaml'
# # --project_directory data --predictors_file data/collection/asthma_children/asthma_children.EnformerLoci.txt --model /project2/haky/Data/enformer/raw --fasta_file /project2/haky/Data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta --parameters_file data/enformer_parameters/enformer_parameters_Asthma_GWAS_asthma_children.json --date 2024-01-31

# opt$personalized_parameters_file <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/config/personalized_base.yaml'

library(yaml)

# read and write the enformer config file
directives <- yaml::yaml.load_file(opt$base_directives)

enformer_parameters_json <- directives$enformer$prediction_directives
# you may change these as appropriate
enformer_parameters_json[['project_dir']] <- normalizePath(opt$project_directory)
enformer_parameters_json[["interval_list_file"]] <- opt$predictors_file
enformer_parameters_json[['prediction_data_name']] <- paste(opt$runname, opt$phenotype, sep = '_')
# enformer_parameters_json[['prediction_id']] <- opt$phenotype
enformer_parameters_json[['date']] <- opt$date
enformer_parameters_json[['model_path']] <- opt$model
enformer_parameters_json[['fasta_file']] <- opt$fasta_file
enformer_parameters_json[['output_dir']] <- opt$project_directory

# chANGE the metadata dir
enformer_parameters_json[['metadata_dir']] <- dirname(opt$parameters_file)
enformer_parameters_json[['copy_aggregation_yaml']] <- opt$copy_aggregation_config

# this ensures that personalized parameters are used
if(!is.null(opt$personalized_parameters_file)){
    enformer_parameters_json[['sequence_source']] <- 'personalized'

    # create vcf files patterns and write to a list
    personalized_parameters <- yaml::yaml.load_file(opt$personalized_parameters_file)
    vcf_path <- file.path(personalized_parameters[['vcf_files']][['folder']], personalized_parameters[['vcf_files']][['files_pattern']])

    chrom_filter <- c(1:22)
    valid_chroms <- c()
    chr_vcfs <- sapply(chrom_filter, function(cc){
        pp <- gsub('\\{\\}', cc, vcf_path)
        if(file.exists(pp)){
            valid_chroms <<- append(valid_chroms, cc)
            return(basename(pp))
        }
    }) 

    names(chr_vcfs) <- paste0('chr', valid_chroms)
    chr_vcfs <- as.list(chr_vcfs)
    #chr_vcfs <- Filter(file.exists, chr_vcfs) |> as.list()
    # nn <- paste0('chr', names(chr_vcfs))
    # names(chr_vcfs) <- nn

    enformer_parameters_json[['vcf_files']][['folder']] <- personalized_parameters[['vcf_files']][['folder']]
    enformer_parameters_json[['vcf_files']][['files']] <- chr_vcfs

    enformer_parameters_json[['individuals']] <- personalized_parameters[['individuals']]
    enformer_parameters_json[['n_individuals']] <- personalized_parameters[['n_individuals']]
    enformer_parameters_json[['batch_individuals']] <- personalized_parameters[['batch_individuals']]

}

yaml::write_yaml(enformer_parameters_json, opt$parameters_file, fileEncoding = "UTF-8", handlers = list(
  logical = function(x) {
    result <- ifelse(x, "true", "false")
    class(result) <- "verbatim"
    return(result)
  })
)

# write(
#     jsonlite::toJSON(enformer_parameters_json, na='null', pretty=TRUE, auto_unbox=T),
#     file=opt$parameters_file
# )


# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript prepare/workflow/scripts/create_enformer_config.R --transcription_factor AR --tissue Breast --base_directives config/enformer_base.yaml --project_directory data/predictions_folder --predictors_file data/predictor_files/AR_Breast.predictors.txt --model "/project2/haky/Data/enformer/raw" --fasta_file "/project2/haky/Data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta" --parameters_file data/prediction_parameters/enformer_parameters_cistrome_AR_Breast.json