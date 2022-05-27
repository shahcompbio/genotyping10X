configfile: "config.yaml"

import pandas as pd

df = pd.read_csv(config["metadata"])

rule all:
    input:
        expand("results/vcf/{sample}/haplotypes_hg38.vcf.gz", sample = df["sample"])

include: "rules/haps_to_vcf.smk"