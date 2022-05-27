def _get_bam_file(wildcards):
    df = scrna[(scrna["isabl_patient_id"] == wildcards.patient) & (scrna["isabl_sample_id"] == wildcards.sample_id)]
    df = df[df["result_type"] == "bam"]
    path = df["result_filepath"].to_list()[0] #this should be checked to make sure there is only one file
    return path

def _get_barcodes_file(wildcards):
    df = scrna[(scrna["isabl_patient_id"] == wildcards.patient) & (scrna["isabl_sample_id"] == wildcards.sample_id)]
    df = df[df["result_type"] == "filtered_matrices"]
    path = df["result_filepath"].to_list()[0] + "/barcodes.tsv.gz"
    return path

rule cellsnp:
    output:
        results_DP_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.OTH.mtx",
        results_sample = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.samples.tsv",
        vcf = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.base.vcf.gz",
    input: 
        vcf =  "results/vcf/{patient}/haplotypes_hg38.vcf.gz",
        bam = _get_bam_file,
    params:
        barcodes = _get_barcodes_file,
        mydir = "results/scrna/counthaps/{patient}/{sample_id}/",
    threads: 40
    resources: mem_mb=200 * 1
    conda: "../envs/cellsnp.yml"
    shell: 
        """
        cellsnp-lite -s {input.bam} \
            -b {params.barcodes} \
            -O {params.mydir} \
            -R {input.vcf} \
            -p {threads} \
            --minMAF 0.0 \
            --minCOUNT 0 \
            --gzip
        """

rule format_per_cell_counts:
    input: 
        results_DP_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.tag.OTH.mtx",
        results_sample = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.samples.tsv",
        vcf = "results/scrna/counthaps/{patient}/{sample_id}/cellSNP.base.vcf.gz",
        vcf_lift = "results/vcf/{patient}/haplotypes_hg38.vcf.gz",
        haplotypes = "results/vcf/{patient}/haplotypes_positions.csv.gz",
        phasedhaplotypes = "results/vcf/{patient}/haplotypes_phased.csv.gz"
    output: 
        alldata = "results/scrna/counthaps/{patient}/{sample_id}/allele_counts.csv.gz",
        alldata_phased_dna = "results/scrna/counthaps/{patient}/{sample_id}/allele_counts_phased.csv.gz",
    threads: 2
    resources: mem_mb=1024 * 50
    singularity: config["singularityR"]
    script: "../scripts/format_vcf_haps.R"

def _get_allele_count_files(wildcards):
    return scrna_md[scrna_md["isabl_patient_id"] == wildcards.patient]["filename"].to_list()

rule combinefiles:
    input: _get_allele_count_files
    output: "results/scrna/counthaps/{patient}/{patient}_allele_counts.csv.gz"
    resources: mem_mb=1024*50
    run:
        dist_list = []
        for f in input:
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output[0], sep=',')