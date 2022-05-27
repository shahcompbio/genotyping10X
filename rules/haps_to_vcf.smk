rule downloadhaplotypes:
    output: "results/vcf/{sample}/haplotypes_positions.csv.gz"
    run:
        import shahlabdata.isabl
        import os
        import time
        from random import random
        #time.sleep(random() * 10 * 60) #seems to crash sometimes if there's a lot of requests at the same time
        os.environ['ISABL_API_URL']='https://isabl.shahlab.mskcc.org/api/v1/'
        haps_res = shahlabdata.isabl.get_results("SCDNA-INFERHAPS", patient_id=wildcards.sample, most_recent=True)
        haps_path = haps_res[haps_res["result_type"] == "haplotypes"].reset_index()["result_filepath"][0]
        os.symlink(haps_path, output[0])

def get_signals_object(wildcards):
    import shahlabdata.isabl
    import os
    import time
    from random import random
    #time.sleep(random() * 10 * 60) #seems to crash sometimes if there's a lot of requests at the same time
    os.environ['ISABL_API_URL']='https://isabl.shahlab.mskcc.org/api/v1/'
    paths = shahlabdata.isabl.get_results("SCDNA-SCHNAPPS", patient_id=wildcards.sample, most_recent=True)
    print(paths)
    rdata_path = paths[paths["result_type"] == "rdatafile"].reset_index()["result_filepath"][0]
    return rdata_path


rule getsignalsphasedhaplotypes:
    output: "results/vcf/{sample}/haplotypes_phased.csv.gz"
    input: get_signals_object
    threads: 1
    resources: mem_mb=1024 * 50
    singularity: config["singularityR"]
    script: "../scripts/get_phased_haplotypes.R"

rule createvcf:
    input: haplotypes = "results/vcf/{sample}/haplotypes_positions.csv.gz",
    output: haplotypesvcf = "results/vcf/{sample}/haplotypes_hg19.vcf",
    threads: 1
    resources: mem_mb=1024 * 50
    singularity: config["singularityR"]
    script: "../scripts/create_vcf_file.R"

rule liftoverinputvcf:
    input: "results/vcf/{sample}/haplotypes_hg19.vcf"
    output: "results/vcf/{sample}/haplotypes_hg38.vcf.gz"
    params:
        chain = "metadata/liftOver/hg19ToHg38.over.chain.gz"
    conda: "../envs/cellsnp.yml"
    shell:
        """
        python scripts/lift_over_vcf.py -c {params.chain} -i {input} -o {output}
        """