#!/usr/bin/env python

import pysam
import sys
import csv
import imp
import os
import re
import argparse
import warnings
import subprocess
import tempfile
import itertools
from Bio import SeqIO
import mafUtils
import bamUtils

NUCLEOTIDES = ["A", "T", "C", "G"]
NA_STRINGS = ["NA", "None", ""]

def read_maf_file(path):
    """Read lines from MAF file. 

    Return a dictionary of lines, keyed by (chromosome, start, ref, alt).
    """
    rows = {}
    with open(path) as f:
        lines = filter(lambda x: not x.startswith("#"), f)
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            key = (row["Chromosome"], 
                   row["Start_Position"], 
                   row["Reference_Allele"], 
                   mafUtils.get_nref_allele(row))
            rows[key] = row
    return rows

def dict_append(d, key, value):
    try:
        d[key].append(value)
    except KeyError:
        d[key] = [value]
    return d

if __name__ == "__main__":

    desc = "Collect VAF from multiple BAM files."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("metadata")
    parser.add_argument("reference")
    args = parser.parse_args()

    reffile = pysam.Fastafile(args.reference)
    reader = csv.DictReader(open(args.metadata), delimiter="\t")

    fields = ["patient", "sample", "chr", "start", "end", "ref", "alt",
              "ref.count", "alt.count", "depth", "gene", "class", "rs",
              "cosmic", "esp", "prot.change"]
    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    # collect files
    maf_files = {}
    bam_files = {}
    for row in reader:
        patient = row["patient"]
        sample = row["sample"]
        snv_maf = row["snv.maf"]
        indel_maf = row["indel.maf"]
        bam_file = row["bam.file"]

        if os.path.exists(snv_maf):
            maf_files = dict_append(maf_files, patient, snv_maf)
        elif snv_maf not in NA_STRINGS:
            warnings.warn("MAF file {} does not exist".format(snv_maf))

        if os.path.exists(indel_maf):
            maf_files = dict_append(maf_files, patient, indel_maf)
        elif indel_maf not in NA_STRINGS:
            warnings.warn("MAF file {} does not exist".format(indel_maf))

        if os.path.exists(bam_file):
            try:
                bam_files[patient][sample] = bam_file
            except KeyError:
                bam_files[patient] = {sample: bam_file}
        elif bam_file not in NA_STRINGS:
            warnings.warn("BAM file {} does not exist".format(bam_file))

    # collect VAF for each patient
    for patient in maf_files.keys():

        # get list of positions and alleles for each of the patient's samples
        all_rows = {}
        for maf_file in maf_files[patient]:
            all_rows.update(read_maf_file(maf_file))
        
        # get the VAF at each of the listed positions
        for sample, bam_file in bam_files[patient].items():
            samfile = pysam.Samfile(bam_file)
            for key, in_row in all_rows.items():

                # fill out an output row
                row = dict.fromkeys(fields)

                # basics
                row["patient"] = patient
                row["sample"] = sample
                row["chr"], row["start"], row["ref"], row["alt"] = key
                row["end"] = int(row["start"]) + len(row["alt"]) - 1

                # get allele counts
                if mafUtils.is_snv(in_row):
                    counts = bamUtils.count_bases(samfile, reffile, row["chr"], int(row["start"]))
                else: # indel
                    counts = bamUtils.count_indels(samfile, reffile, row["chr"], int(row["start"]), row["ref"], row["alt"])
                row["ref.count"] = counts[row["ref"]]
                row["alt.count"] = counts[row["alt"]]
                row["depth"] = sum(counts.values())

                # info from MAF file
                row["gene"] = in_row["Hugo_Symbol"]
                row["class"] = in_row["Variant_Classification"]
                row["cosmic"] = ",".join(re.findall("COSM\d+", in_row["dbSNP_RS"]))
                row["rs"] = ",".join(re.findall("rs\d+", in_row["dbSNP_RS"]))
                row["esp"] = ",".join(re.findall("TMP_ESP_[\dXY]{1,2}_\d+", in_row["dbSNP_RS"]))
                row["prot.change"] = mafUtils.get_protein_change(in_row)

                writer.writerow(row)
