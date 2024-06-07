snakemake -s snakefile.smk --configfile config/pipeline.yaml -np


screen
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/project2/haky/temi/software/homer/bin
snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ >> run.txt 2>&1  


snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ >> run.txt 2>&1

snakemake -s snakefile.smk --configfile config/pipeline_2023-12-01_enformer_vs_tfpred.yaml --profile profiles/simple/ --resources load=50 >> run_2023-12-01_enformer_vs_tfpred.txt 2>&1

snakemake -s snakefile.smk --configfile config/pipeline_pcr.yaml --profile profiles/simple/ -np

snakemake -s snakefile.smk --configfile config/pipeline_pcr.yaml --profile profiles/simple/ -F

pp = '/project/haky/users/temi/projects/TFXcan-snakemake/data/BC_GWAS/enpact_predictions'

pp = '/scratch/midway3/temi/enpact_predictions'

snakemake.io.glob_wildcards(os.path.join(pp, '{phenotype}', '{idi}.*.csv.gz')).idi

snakemake.io.glob_wildcards(os.path.join(pp, '{phenotype}', '{idi}.{phenop}.aggByCollect.2024-06-04.csv.gz'))

snakemake -s snakefile.smk --configfile config/pipeline_bc2.yaml --profile profiles/simple/ -np


snakemake -s snakefile.smk --configfile config/pipeline_bc.yaml --profile profiles/simple/ -R create_enformer_configuration -np

snakemake -s snakefile.smk --configfile config/pipeline_bc.yaml --profile profiles/simple/ -R create_enformer_configuration -np > txt.out

snakemake -s snakefile.smk --configfile config/pipeline_bc.yaml --profile profiles/simple/ --allowed-rules create_enformer_configuration -np > txt.out

snakemake -s snakefile.smk --configfile config/pipeline_bc.yaml --profile profiles/simple/ -np