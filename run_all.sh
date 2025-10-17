# main actual analyses:
Rscript -e 'rmarkdown::render("01_toySim/grouped_fdr_sim.Rmd")'
Rscript -e 'rmarkdown::render("02_MSdata/01_ms_ihw.Rmd")'
Rscript -e 'rmarkdown::render("02_MSdata/02_ms_nbins.Rmd")'
Rscript -e 'rmarkdown::render("03_simulation1/01_simulation.Rmd")'
Rscript -e 'rmarkdown::render("03_simulation1/02_sim1_nbins.Rmd")'
# figures:
Rscript -e 'rmarkdown::render("05_figures/Figure1_example.Rmd")'
Rscript -e 'rmarkdown::render("05_figures/Figure3_toySims.Rmd")'
Rscript -e 'rmarkdown::render("05_figures/Figure4_mainResults.Rmd")'
Rscript -e 'rmarkdown::render("05_figures/Figure5_combinations.Rmd")'
Rscript -e 'rmarkdown::render("05_figures/Figure6_logFC.Rmd")'
# for the last figure (added sims in revision), first 
# adjust your environment in 04_simulation2 (see readme there), then
# run the snakemake
# snakemake 04_simulation2/Snakefile
# Rscript -e 'rmarkdown::render("05_figures/Figure7_sampleSize.Rmd")'
