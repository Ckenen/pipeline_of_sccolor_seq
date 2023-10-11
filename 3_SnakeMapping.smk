#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/mapping"

GENOME_SPLICE_MMI = "/home/chenzonggui/species/mix_hg38_mm10/GRCh38.p13_GRCm38.p6.mm2.splice.mmi"
TRANSCRIPT_BED = "/home/chenzonggui/species/mix_hg38_mm10/gencode.v39_vM25.transcripts.bed"

rule all:
    input:
        expand(outdir + "/minimap2/{srr}.bam", srr=srr_list),
        expand(outdir + "/minimap2/{srr}.flagstat", srr=srr_list),
        expand(outdir + "/filtered/{srr}.bam", srr=srr_list),
        expand(outdir + "/filtered/{srr}.flagstat", srr=srr_list),
        expand(outdir + "/stat_clip/{srr}.bam", srr=srr_list),
        expand(outdir + "/stat_clip/{srr}.flagstat", srr=srr_list),


rule minimap2:
    input:
        fq = "results/tallynn/remove_polya/{srr}.fastq.gz",
        mmi = GENOME_SPLICE_MMI,
        bed = TRANSCRIPT_BED
    output:
        bam = outdir + "/minimap2/{srr}.bam"
    log:
        outdir + "/minimap2/{srr}.log"
    threads:
        24
    shell: 
        """(
        minimap2 -ax splice -u f -Y --MD --junc-bed {input.bed} -t {threads} {input.mmi} {input.fq} \
            | samtools view -q 60 -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = outdir + "/filtered/{srr}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^(hg|mm)_chr([0-9]+|[XY])$"' -q 30 -m 200 -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/stat_clip/{srr}.bam",
        txt = outdir + "/stat_clip/{srr}.tsv"
    log:
        outdir + "/stat_clip/{srr}.log"
    threads:
        4
    shell:
        """
        nasctools StatClip --max-clip 5 -s {output.txt} -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """