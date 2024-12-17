

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))

mdir <- '/project/haky/users/temi/projects/TFXcan-snakemake/output/lEnpact/pc_risk/models'

filtered_db <- file.path(mdir, 'filtered_db', 'predict_db_pc_risk_filtered.db')
unfiltered_cov <- file.path(mdir, 'database', 'predict_db_pc_risk.txt.gz')
filtered_cov <- args[3]


# load db
driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from extra')
dbDisconnect(in_conn)

# load covariances
all_covs <- read.table(unfiltered_cov, header = TRUE, stringsAsFactors = FALSE)

# Filter out models with low performance
all_covs <- all_covs %>% 
    filter(GENE %in% model_summaries$gene)

# Write out the covariance file filtered
write.table(all_covs, file = filtered_cov, sep = "\t",
            row.names = FALSE, quote = FALSE)



#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=summary_TFXcan
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/summary_TFXcan.out
#SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/summary_TFXcan.err
#SBATCH --time=06:00:00



# SBATCH --partition=bigmem

# module load openmpi
# module load parallel




 sbatch stfxcan.sbatch pc_risk output/lEnpact/pc_risk/models/filtered_db/predict_db_pc_risk_filtered.db output/summary/pc_risk/pc_risk.enpactScores.spredixcan.csv /project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/processed_sumstats/pc_risk .*.sumstats.txt.gz output/lEnpact/pc_risk/models/database/Covariances.varID.txt


 squeue -u temi | awk '{ print $1 }'
 squeue -u temi | grep "beagle3" | grep "smk-pred" | awk '{ print $1 }' | tail -n+2 | xargs scancel

 sinteractive --account=pi-haky --partition=beagle3 --gres=gpu:1

 ['N', 'R', 'W', 'Y', 'S', 'Q', 'K', 'M']




 Loading Loci: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 48696/48696 [06:45<00:00, 120.01it/s]
Loading Loci: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3063/3063 [00:06<00:00, 441.21it/s]
Training Set Size:  48696
Validation Set Size:  3063
Epoch   Iteration       Training Time   Validation Time Training MNLL   Training Count MSE      Validation MNLL Validation Profile Pearson      Validation Count Pearson        Validation Count MSE    Saved?
Traceback (most recent call last):
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/bin/bpnet", line 432, in <module>
    model.fit(training_data, optimizer, X_valid=valid_sequences,
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/bpnetlite/bpnet.py", line 407, in fit
    y_profile, y_counts = self(X, X_ctl)
                          ^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/module.py", line 1532, in _wrapped_call_impl
    return self._call_impl(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/module.py", line 1541, in _call_impl
    return forward_call(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/bpnetlite/bpnet.py", line 301, in forward
    y_profile = self.fconv(X_w_ctl)[:, :, start:end]
                ^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/module.py", line 1532, in _wrapped_call_impl
    return self._call_impl(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/module.py", line 1541, in _call_impl
    return forward_call(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/conv.py", line 310, in forward
    return self._conv_forward(input, self.weight, self.bias)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/project/haky/users/temi/software/conda_envs/bpnet-lite/lib/python3.12/site-packages/torch/nn/modules/conv.py", line 306, in _conv_forward
    return F.conv1d(input, weight, bias, self.stride,
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Given groups=1, weight of size [1, 66, 75], expected input[64, 65, 2114] to have 66 channels, but got 65 channels instead