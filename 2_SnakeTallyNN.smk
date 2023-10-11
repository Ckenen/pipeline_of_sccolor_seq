#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/tallynn"

rule all:
    input:
        # expand(outdir + "/complement_polyA/{srr_batch}.fastq", srr_batch=srr_batch_list),
        # expand(outdir + "/identify_perfect/{srr_batch}_unambiguous_barcode_R1.fastq", srr_batch=srr_batch_list),
        # expand(outdir + "/identify_perfect_merged/{srr}_unambiguous_barcode_R1.fastq", srr=srr_list),
        expand(outdir + "/umitools/{srr}.whitelist.txt", srr=srr_list),
        # expand(outdir + "/correct_barcode/{srr_batch}_unambiguous_fixed_barcode_R1.fastq", srr_batch=srr_batch_list),
        # expand(outdir + "/correct_barcode_merged/{srr}_unambiguous_fixed_barcode_R1.fastq", srr=srr_list),
        expand(outdir + "/merge_full/{srr}_R1.fastq.gz", srr=srr_list),
        # expand(outdir + "/extract_umibc_readname/{srr}.fastq.1.gz", srr=srr_list),
        #expand(outdir + "/umitools_2/{srr}.whitelist.txt", srr=srr_list),
        expand(outdir + "/umitools_extract/{srr}_R1.fastq.gz", srr=srr_list),
        expand(outdir + "/remove_adapter/{srr}.fastq.gz", srr=srr_list),
        expand(outdir + "/reverse_complement/{srr}.fastq.gz", srr=srr_list),
        expand(outdir + "/remove_polya/{srr}.fastq.gz", srr=srr_list),


rule complement_polyA:
    input:
        fq = "results/prepare/split/{srr}/{batch}.fastq.gz"
    output:
        fq = outdir + "/complement_polyA/{srr}/{batch}.fastq"
    log:
        outdir + "/complement_polyA/{srr}/{batch}.log"
    shell:
        """
        ./scripts/complement_polyA.py --infile {input} --outname {output} &> {log}
        """

rule identify_perfect:
    input:
        fq = rules.complement_polyA.output.fq
    output:
        fq1 = outdir + "/identify_perfect/{srr}/{batch}_unambiguous_barcode_R1.fastq",
        fq2 = outdir + "/identify_perfect/{srr}/{batch}_unambiguous_barcode_R2.fastq",
        fq3 = outdir + "/identify_perfect/{srr}/{batch}_ambiguous_barcode_R1.fastq",
        fq4 = outdir + "/identify_perfect/{srr}/{batch}_ambiguous_barcode_R2.fastq",
        txt = outdir + "/identify_perfect/{srr}/{batch}.whitelist.txt"
    log:
        outdir + "/identify_perfect/{srr}/{batch}.log"
    params:
        prefix = outdir + "/identify_perfect/{srr}/{batch}"
    shell:
        """
        set +u; source activate tallynn
        python software/TallyNN-master/tallynn/python/identify_perfect_nano.py \
            --outname={params.prefix} --infile={input.fq} --whitelist={output.txt} &> {log}
        """

rule merge_unambigous:
    input:
        fqs1 = lambda wildcards: [outdir + "/identify_perfect/{srr}/%s_unambiguous_barcode_R1.fastq" % x for x in srr_to_batch[wildcards.srr]],
        fqs2 = lambda wildcards: [outdir + "/identify_perfect/{srr}/%s_unambiguous_barcode_R2.fastq" % x for x in srr_to_batch[wildcards.srr]]
    output:
        fq1 = outdir + "/identify_perfect_merged/{srr}_unambiguous_barcode_R1.fastq",
        fq2 = outdir + "/identify_perfect_merged/{srr}_unambiguous_barcode_R2.fastq"
    shell:
        """
        cat {input.fqs1} > {output.fq1}
        cat {input.fqs2} > {output.fq2}
        """

rule umitools_whitelist:
    input:
        fq = rules.merge_unambigous.output.fq1
    output:
        txt = outdir + "/umitools/{srr}.whitelist.txt"
    log:
        outdir + "/umitools/{srr}.whitelist.log"
    params:
        cell_number = lambda wildcards: dat[dat["Run"] == wildcards.srr]["Cell_number"].values[0]
    shell:
        """
        set +u; source activate tallynn
        umi_tools whitelist --stdin={input.fq} --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN \
            --set-cell-number={params.cell_number} -L {log} > {output.txt}
        """

rule correct_barcode:
    input:
        fq1 = rules.identify_perfect.output.fq3,
        fq2 = rules.identify_perfect.output.fq4,
        txt = rules.umitools_whitelist.output.txt
    output:
        fq1 = outdir + "/correct_barcode/{srr}/{batch}_unambiguous_fixed_barcode_R1.fastq",
        fq2 = outdir + "/correct_barcode/{srr}/{batch}_unambiguous_fixed_barcode_R2.fastq"
    params:
        prefix = outdir + "/correct_barcode/{srr}/{batch}"
    shell:
        """
        set +u; source activate tallynn
        python ./scripts/correct_barcode_nano.py \
            --whitelist={input.txt} --read1={input.fq1} --read2={input.fq2} \
            --outname={params.prefix} --distance=4
        """

rule merge_corrected:
    input:
        fqs1 = lambda wildcards: [outdir + "/correct_barcode/{srr}/%s_unambiguous_fixed_barcode_R1.fastq" % x for x in srr_to_batch[wildcards.srr]],
        fqs2 = lambda wildcards: [outdir + "/correct_barcode/{srr}/%s_unambiguous_fixed_barcode_R2.fastq" % x for x in srr_to_batch[wildcards.srr]]
    output:
        fq1 = outdir + "/correct_barcode_merged/{srr}_unambiguous_fixed_barcode_R1.fastq",
        fq2 = outdir + "/correct_barcode_merged/{srr}_unambiguous_fixed_barcode_R2.fastq"
    shell:
        """
        cat {input.fqs1} > {output.fq1}
        cat {input.fqs2} > {output.fq2}
        """

rule merge_full:
    input:
        fq1 = rules.merge_unambigous.output.fq1,
        fq2 = rules.merge_unambigous.output.fq2,
        fq3 = rules.merge_corrected.output.fq1,
        fq4 = rules.merge_corrected.output.fq2
    output:
        fq1 = outdir + "/merge_full/{srr}_R1.fastq.gz",
        fq2 = outdir + "/merge_full/{srr}_R2.fastq.gz"
    threads:
        8
    shell:
        """
        cat {input.fq1} {input.fq3} | pigz -p {threads} -c > {output.fq1}
        cat {input.fq2} {input.fq4} | pigz -p {threads} -c > {output.fq2}
        """

rule extract_umibc_readname:
    input:
        fq1 = rules.merge_full.output.fq1,
        fq2 = rules.merge_full.output.fq2
    output:
        fq1 = outdir + "/extract_umibc_readname/{srr}.fastq.1.gz",
        fq2 = outdir + "/extract_umibc_readname/{srr}.fastq.2.gz"
    params:
        prefix = outdir + "/extract_umibc_readname/{srr}"
    shell:
        """
        set +u; source activate tallynn
        python software/TallyNN-master/tallynn/python/extract_umibc_readname.py \
            --read1 {input.fq1} --read2 {input.fq2} --outname {params.prefix}
        """

rule umitools_whitelist_2:
    input:
        fq = rules.merge_full.output.fq1
    output:
        txt = outdir + "/umitools_2/{srr}.whitelist.txt"
    log:
        outdir + "/umitools_2/{srr}.whitelist.log"
    params:
        cell_number = lambda wildcards: dat[dat["Run"] == wildcards.srr]["Cell_number"].values[0]
    shell:
        """
        set +u; source activate tallynn
        umi_tools whitelist --stdin={input.fq} --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN \
            --set-cell-number={params.cell_number} -L {log} > {output.txt}
        """

rule umitools_extract_bc:
    input:
        fq1 = rules.merge_full.output.fq1,
        fq2 = rules.merge_full.output.fq2,
        txt = rules.umitools_whitelist_2.output.txt
    output:
        fq1 = outdir + "/umitools_extract/{srr}_R1.fastq.gz",
        fq2 = outdir + "/umitools_extract/{srr}_R2.fastq.gz"
    log:
        outdir + "/umitools_extract/{srr}.log"
    shell:
        """
        set +u; source activate tallynn
        umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN --stdin {input.fq1} \
            --stdout={output.fq1} --read2-in {input.fq2} --read2-out={output.fq2} --whitelist={input.txt} &> {log}
        """

rule remove_adapter:
    input:
        fq = rules.umitools_extract_bc.output.fq2
    output:
        fq = outdir + "/remove_adapter/{srr}.fastq.gz",
    log:
        outdir + "/remove_adapter/{srr}.log"
    threads:
        12
    shell:
        """
        cutadapt -j {threads} -e 0.2 -m 200 -a ACTCTGCGTTGATACCACTGCTT -o {output.fq} {input.fq} &> {log}
        """

rule reverse_complement:
    input:
        fq = rules.remove_adapter.output.fq
    output:
        fq = outdir + "/reverse_complement/{srr}.fastq.gz"
    threads:
        8
    shell:
        """
        zcat {input.fq} | ./scripts/reverse_fastq.py | pigz -p {threads} -c > {output.fq}
        """

rule remove_polya:
    input:
        fq = rules.reverse_complement.output.fq
    output:
        fq = outdir + "/remove_polya/{srr}.fastq.gz"
    log:
        outdir + "/remove_polya/{srr}.log"
    threads:
        12
    shell:
        """
        cutadapt -j {threads} -e 0.2 -m 200 -a AAAAAAAAAAAAAAA -o {output.fq} {input.fq} &> {log}
        """
