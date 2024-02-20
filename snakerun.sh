snakemake -s snakefile.smk --configfile config/pipeline.yaml -np


screen
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/project2/haky/temi/software/homer/bin
snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ >> run.txt 2>&1  


snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ >> run.txt 2>&1

snakemake -s snakefile.smk --configfile config/pipeline_2023-12-01_enformer_vs_tfpred.yaml --profile profiles/simple/ --resources load=50 >> run_2023-12-01_enformer_vs_tfpred.txt 2>&1