# Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis

This is a repository accompanying this [note](https://www.biorxiv.org/content/10.1101/2025.04.15.648932v1), with the code and data to reproduce the analyses.

## Installing the modified version of muscat

To use the methods (until the update is ported to the bioconductor release), install the bulkIhw development branch of muscat, i.e.:

```{r}
BiocManager::install("HelenaLC/muscat", ref="bulkIhw", build_vignettes=TRUE)
```

Then see the related documentation with:

```{r}
vignette("bbhw", package="muscat")
```

(The code in the present repository assumes that this version of the muscat package is installed)

## Structure of the repository

- `01_toySim` contains the scripts generating the toy simulations.
- `02_MSdata` contains the scripts doing the bulk aggregation, downsampling, DEA and bbhw for the real (human MS) data.
- `03_simulation1` contains the scripts to generate the first set of simulations, and perform DEA and bbhw on them.
- `04_simulation2` contains a snakemake to generate the second set of simulations, and perform DEA and bbhw on them.
- `05_figures` contains the scripts to generate the figures, based on the outputs of the previous scripts. Make sure to run the scripts in the previous folders (in alphabetical order) first.
