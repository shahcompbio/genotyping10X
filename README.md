# Genotyping SNPs in 10X bam files

## What is this?

Snakemake pipeline to genotype SNPs in 10X bam files. Pipeline takes output of inferhaps, converts it to vcf, then uses [cellSNP](https://github.com/single-cell-genetics/cellsnp-lite) to count reads in single cell bam files. Output is a merged file with read counts per SNP in all cells for each patient of interest for either scRNA or scATAC or both.

## Running pipeline

1. Set up a snakemake profile for lsf job submissions, see [here](https://github.com/Snakemake-Profiles/lsf).
2. Build conda env `envs/isablapi.yml`
3. Change the `metadata/atac.csv` and `metadata/rna.csv` to include the samples of interest. These should be isabl_patient_id's
4. Submit the pipeline: `bsub < launch_job_lsf.sh`

Outputs:
* `results/scrna/counthaps/{patient}/{patient}_allele_counts.csv.gz` contains the read counts for each SNP in every cell, alleles are phased according to INFER-HAPS output (ie population based phasing using SHAPEIT).
* `results/scrna/counthaps/{patient}/{patient}_allele_counts_phased.csv.gz` contains the read counts for each SNP in every cell, alleles are phased according to the refined phasing provided by signals. Note, that this file will contain fewer SNPs as signals drops any SNPs in noisy bins.

## Notes

* Pipeline assumes the DLP results are in hg19 and does some liftovers for compatibility with scRNA. This could be changed if using hg38 results directly.
* Currently this is not working for scATACseq, not sure why. Attempt at running for SA1049 can be found here: `/juno/work/shah/users/william1/projects/smallprojects/genotyping10X`. The allele_counts files are empty.
* Pipeline relies on the output from INFER-HAPS in the DLP pipeline, therefore only samples with normal bulk WGS will work.
* Pipeline also relies on results from signals (`SCDNA-SCHNAPPS` in isabl), this is only used to generate the phased allele table so could easily be removed if it's not possible to run signals.
