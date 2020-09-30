# TRAP datasets

We acquired the deposited raw sequencing fastq files from the online repositories [GSE100579](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100579) (10 sequencing runs from [Marrocco J. et al. 2017](https://doi.org/10.1038/s41467-017-01014-4) and [Gray J. D. et al. 2018](https://dx.doi.org/10.1038%2Fmp.2016.219)) and [GSE131972](www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131972) (10 sequencing runs from [Marrocco J. et al. 2019](https://doi.org/10.3389/fnbeh.2019.00157)) and used [kallisto](https://pachterlab.github.io/kallisto/about) for the pseudoalignment of reads on the GENCODE M17 transcriptome, with an estimated fragment length of 200 ±20.

The uniformly-processed data for the three publications is provided as a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html): [AllData.kallisto.SE.rds](AllData.kallisto.SE.rds).

In addition, since in the lastest publications the authors used [salmon](https://combine-lab.github.io/salmon/) for quantification, we also provide a salmon-based quantification for this data: [GSE131972.salmonv2.SE.rds](GSE131972.salmonv2.SE.rds).

