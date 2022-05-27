configfile: "config.yaml"

import pandas as pd
import shahlabdata.isabl
import os
import time
from random import random
os.environ['ISABL_API_URL']='https://isabl.shahlab.mskcc.org/api/v1/'

samples_rna = pd.read_csv(config["samples_rna"])
samples_atac = pd.read_csv(config["samples_atac"])

#get metadata from isabl
scrna = shahlabdata.isabl.get_results("CELLRANGER", most_recent=True)
scrna = scrna[scrna["isabl_patient_id"].isin(samples_rna["sample"])]
scrna = scrna[scrna["isabl_aliquot_id"] != "Unspecified Aliquot"]
scrna_md = scrna[["isabl_patient_id", "isabl_sample_id", "isabl_aliquot_id"]].drop_duplicates()
scrna_md["filename"] = "results/scrna/counthaps/" + scrna_md["isabl_patient_id"] + "/" + scrna_md["isabl_sample_id"] + "/allele_counts.csv.gz"


scatac = shahlabdata.isabl.get_results("CELLRANGER_ATAC", most_recent=True)
scatac = scatac[scatac["isabl_patient_id"].isin(samples_atac["sample"])]
scatac = scatac[scatac["isabl_aliquot_id"] != "Unspecified Aliquot"]
scatac_md = scatac[["isabl_patient_id", "isabl_sample_id", "isabl_aliquot_id"]].drop_duplicates()
scatac_md["filename"] = "results/scatac/counthaps/" + scatac_md["isabl_patient_id"] + "/" + scatac_md["isabl_sample_id"] + "/allele_counts.csv.gz"
print(scatac_md.head())

rule all:
    input:
        expand("results/scrna/counthaps/{patient}/{patient}_allele_counts.csv.gz", patient = samples_rna["sample"]),
        expand("results/scatac/counthaps/{patient}/{patient}_allele_counts.csv.gz", patient = samples_atac["sample"])

include: "rules/haps_to_vcf.smk"
include: "rules/scRNA.smk"
include: "rules/scATAC.smk"