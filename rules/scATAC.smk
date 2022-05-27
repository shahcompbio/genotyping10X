def _get_bam_file_atac(wildcards):
    df = scatac[(scatac["isabl_patient_id"] == wildcards.patient) & (scatac["isabl_sample_id"] == wildcards.sample_id)]
    df = df[df["result_type"] == "bam"]
    path = df["result_filepath"].to_list()[0] #this should be checked to make sure there is only one file
    return path

def _get_barcodes_file_atac(wildcards):
    df = scatac[(scatac["isabl_patient_id"] == wildcards.patient) & (scatac["isabl_sample_id"] == wildcards.sample_id)]
    df = df[df["result_type"] == "filtered_peak_bc_matrix"]
    path = df["result_filepath"].to_list()[0].replace(".h5", "") + "/barcodes.tsv"
    print(path)
    return path

rule cellsnp_atac:
    output:
        results_DP_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.OTH.mtx",
        results_sample = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.samples.tsv",
        vcf = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.base.vcf.gz",
    input: 
        vcf =  "results/vcf/{patient}/haplotypes_hg38.vcf.gz",
        bam = _get_bam_file_atac,
    params:
        barcodes = _get_barcodes_file_atac,
        mydir = "results/scatac/counthaps/{patient}/{sample_id}/",
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

rule format_per_cell_counts_atac:
    input: 
        results_DP_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.tag.OTH.mtx",
        results_sample = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.samples.tsv",
        vcf = "results/scatac/counthaps/{patient}/{sample_id}/cellSNP.base.vcf.gz",
        vcf_lift = "results/vcf/{patient}/haplotypes_hg38.vcf.gz",
        haplotypes = "results/vcf/{patient}/haplotypes_positions.csv.gz",
        phasedhaplotypes = "results/vcf/{patient}/haplotypes_phased.csv.gz"
    output: 
        alldata = "results/scatac/counthaps/{patient}/{sample_id}/allele_counts.csv.gz",
        alldata_phased_dna = "results/scatac/counthaps/{patient}/{sample_id}/allele_counts_phased.csv.gz",
    threads: 2
    resources: mem_mb=1024 * 50
    singularity: config["singularityR"]
    script: "../scripts/format_vcf_haps.R"

def _get_allele_count_files_atac(wildcards):
    return scatac_md[scatac_md["isabl_patient_id"] == wildcards.patient]["filename"].to_list()

rule combinefilesatac:
    input: _get_allele_count_files_atac
    output: "results/scatac/counthaps/{patient}/{patient}_allele_counts.csv.gz"
    resources: mem_mb=1024*50
    run:
        dist_list = []
        for f in input:
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output[0], sep=',')