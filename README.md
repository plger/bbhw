# Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis

This is a repository accompanying this [note](https://www.biorxiv.org/content/10.1101/2025.04.15.648932v1), with the code and data to reproduce the analyses.

To use the methods (until the update is ported to the bioconductor release), install the bulkIhw development branch of muscat, i.e.:

```{r}
BiocManager::install("HelenaLC/muscat", ref="bulkIhw", build_vignettes=TRUE)
```

Then see the related documentation with:

```{r}
vignette("bbhw", package="muscat")
```

(The code in the present repository assumed that this version of the muscat package is installed)
