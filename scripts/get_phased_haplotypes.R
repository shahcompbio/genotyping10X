x <- readRDS(snakemake@input[[1]])
data.table::fwrite(x$hscn$haplotype_phasing,snakemake@output[[1]])