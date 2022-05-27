## Genotyping SNPs in 10X bam files

Snakemake pipeline to genotype SNPs in 10X bam files. Pipeline takes output of inferhaps, converts it to vcf, then uses (cellSNP)[https://github.com/single-cell-genetics/cellsnp-lite] to count reads in single cell bam files. Output is a merged file with read counts per SNP in all cells.