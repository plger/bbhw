# Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis

This is a repository accompanying the above-titled note.

To use the methods (until the update is ported to the bioconductor release), install the bulkIhw development branch of muscat, i.e.:

```{r}
BiocManager::install("HelenaLC/muscat", ref="bulkIhw", build_vignettes=TRUE)
```

Then see the related documentation with:

```{r}
vignette("bbhw", package="muscat")
```
