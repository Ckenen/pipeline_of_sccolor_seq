#!/usr/bin/env runsnakemake
import glob
import pandas as pd
dat = pd.read_csv("data/SraRunTable.scCOLORseq.csv")
srr_list = dat["Run"]
srr_list = [
    "SRR13120810",
    "SRR13397411",
    "SRR13397412",
    "SRR13397413",
    "SRR13397414",
    "SRR13509039"
    ]


srr_batch_list = []
srr_to_batch = dict()
for srr in srr_list:
    srr_to_batch[srr] = list()
    d = "results/prepare/split/%s" % srr
    if os.path.exists(d):
        for path in sorted(glob.glob(d + "/*.fastq.gz")):
            batch = path.split("/")[-1][:-9]
            srr_batch_list.append("%s/%s" % (srr, batch))
            srr_to_batch[srr].append(batch)
print("SRR batch:", len(srr_batch_list))

samples = srr_list