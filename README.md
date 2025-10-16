# Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis

This is a repository accompanying this [note](https://www.biorxiv.org/content/10.1101/2025.04.15.648932v1), with the code and data to reproduce the analyses.

## Installing the modified version of muscat

The bbhw-related functions are now in the devel bioconductor version of muscat (and will soon be in the release version). You may also install the latest from github:

```{r}
BiocManager::install("HelenaLC/muscat")
```

If you wish to use the exact version of muscat with which the analyses were run, use this:

```{r}
BiocManager::install("HelenaLC/muscat", ref="ff6d5a020c2171deaae67c9686a7974e296ca78d", build_vignettes=FALSE)
```

### Docker image

We also provide a docker image (`plger/bbhw` on dockerhub) to run the scripts:

```
# navigate to the directory you want to work at
# clone the repository:
git clone https://github.com/plger/bbhw.git
# launch the docker image, mount the pwd, and run all scripts:
docker run -v $PWD:/home/build/package plger/bbhw bash run_all.sh
```

## Structure of the repository

- `01_toySim` contains the scripts generating the toy simulations.
- `02_MSdata` contains the scripts doing the bulk aggregation, downsampling, DEA and bbhw for the real (human MS) data.
- `03_simulation1` contains the scripts to generate the first set of simulations, and perform DEA and bbhw on them.
- `04_simulation2` contains a snakemake to generate the second set of simulations, and perform DEA and bbhw on them. Note that this is sitting in a separate repository, consult the readme there.
- `05_figures` contains the scripts to generate the figures, based on the outputs of the previous scripts. Make sure to run the scripts in the previous folders (in alphabetical order) first.
