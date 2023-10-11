#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/prepare"

rule all:
    input:
        # expand(outdir + "/sra/{srr}.sra", srr=srr_list),
        # expand(outdir + "/fastqs/{srr}.fastq.gz", srr=srr_list),
        expand(outdir + "/split/{srr}", srr=srr_list),

rule prefetch:
    output:
        sra = outdir + "/sra/{srr}.sra"
    shell:
        """
        prefetch --max-size 200000000 -o {output.sra} {wildcards.srr}
        """
    
rule fasterq_dump:
    input:
        sra = rules.prefetch.output.sra
    output:
        fq1 = temp(outdir + "/fastqs/{srr}.fastq"),
        fq2 = outdir + "/fastqs/{srr}.fastq.gz"
    threads:
        6
    shell:
        """
        fasterq-dump -e {threads} -o {output.fq1} {input.sra}
        pigz -p {threads} -c {output.fq1} > {output.fq2}
        """

rule split:
    input:
        fq = rules.fasterq_dump.output.fq2
    output:
        directory(outdir + "/split/{srr}")
    threads:
        8
    shell:
        """
        mkdir -p {output}
        zcat {input.fq} | split --additional-suffix=.fastq -l 1000000 - {output}/
        pigz -p {threads} {output}/*.fastq
        """