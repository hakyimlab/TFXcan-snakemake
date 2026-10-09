# Installs the R packages of the tfxcan environment that are not on conda-forge (used by factorization).
# Run inside the activated environment: Rscript container/envs/install_r_packages.R

repos <- "https://cloud.r-project.org"
pkgs <- c(ebnm = "1.1-42", flashier = "1.0.7")   # the versions the Midway runs used

for (p in names(pkgs)) {
    if (!requireNamespace(p, quietly = TRUE) || as.character(packageVersion(p)) != pkgs[[p]]) {
        remotes::install_version(p, version = pkgs[[p]], repos = repos, upgrade = "never")
    }
}
for (p in names(pkgs)) {
    stopifnot(requireNamespace(p, quietly = TRUE))
    cat(p, as.character(packageVersion(p)), "\n")
}
