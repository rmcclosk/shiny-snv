#!/usr/bin/env python

import pysam
import sys
import csv
import imp
import os
import re
import argparse
import warnings

NUCLEOTIDES = ["A", "T", "C", "G"]
NA_STRINGS = ["NA", "None", ""]

def count_bases_pileup(pileup, position, min_baseq=15, min_mapq=20):
    """
    Count the number of times each nucleotide occurs at a given position in a
    pileup object, subject to minimum base quality and map quality constraints.
    Return a dictionary keyed by nucleotides.
    """
    counts = dict.fromkeys(NUCLEOTIDES, 0)
    for x in pileup:
        if position == x.pos + 1:
            for read in x.pileups:
                dup = read.alignment.is_duplicate
                qc_fail = read.alignment.is_qcfail
                low_mapq = read.alignment.mapq < min_mapq
                if not (dup or qc_fail or low_mapq):
                    base_qual = ord(read.alignment.qual[read.qpos])-33
                    if base_qual >= min_baseq:
                        base = read.alignment.seq[read.qpos]
                        counts[base] += 1

    return counts

def count_bases(samfile, reffile, chrom, pos):
    """
    Count the number of times each nucleotide occurs at a given position in a
    bam file. Return a dictionary keyed by nucleotides.
    """
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
    pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    return count_bases_pileup(pileup, pos)

def read_maf_file(path):
    """Read lines from MAF file. 

    Return a dictionary of lines, keyed by (chromosome, start, ref, alt).
    """
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = (row["Chromosome"], 
                   row["Start_Position"], 
                   row["Reference_Allele"], 
                   row["Tumor_Seq_Allele1"])
            rows[key] = row
    return rows

if __name__ == "__main__":

    desc = "Collect VAF from multiple BAM files."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("metadata", metavar="M")
    parser.add_argument("reference", metavar="R")
    args = parser.parse_args()

    reffile = pysam.Fastafile(args.reference)
    reader = csv.DictReader(open(args.metadata), delimiter="\t")

    fields = ["patient", "sample", "chrom", "pos", "ref", "alt", "ref.count", 
              "alt.count", "depth", "gene", "class", "rs", "cosmic", "esp",
              "prot.change"]
    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    # collect files
    maf_files = {}
    bam_files = {}
    for row in reader:
        patient = row["patient"]
        sample = row["sample"]
        maf_file = row["maf.file"]
        bam_file = row["bam.file"]

        if os.path.exists(maf_file):
            try:
                maf_files[patient].append(maf_file)
            except KeyError:
                maf_files[patient] = [maf_file]
        elif maf_file not in NA_STRINGS:
            warnings.warn("MAF file {} does not exist".format(maf_file))

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
                row["chrom"], row["pos"], row["ref"], row["alt"] = key

                # get allele counts
                counts = count_bases(samfile, reffile, row["chrom"], int(row["pos"]))
                row["ref.count"] = counts[in_row["Reference_Allele"]]
                row["alt.count"] = counts[in_row["Tumor_Seq_Allele1"]]
                row["depth"] = sum(counts.values())

                # info from MAF file
                row["gene"] = in_row["Hugo_Symbol"]
                row["class"] = in_row["Variant_Classification"]
                row["cosmic"] = ",".join(re.findall("COSM\d+", in_row["dbSNP_RS"]))
                row["rs"] = ",".join(re.findall("rs\d+", in_row["dbSNP_RS"]))
                row["esp"] = ",".join(re.findall("TMP_ESP_[\dXY]{1,2}_\d+", in_row["dbSNP_RS"]))
                row["prot.change"] = in_row["Protein_Change"].replace("p.", "")

                writer.writerow(row)
