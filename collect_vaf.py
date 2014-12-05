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

NUCLEOTIDES = ["A", "T", "C", "G"]
NA_STRINGS = ["NA", "None", ""]

def num_mismatches(ref, query):
    """Number of matches between query and reference."""
    mismatches = 0
    for i in range(len(ref)):
        if query[i] != ref[i]:
            mismatches += 1
    return mismatches

def multi_align(ref1, ref2, queries):
    """
    Do MSA of several query sequences to two references. Return the reference
    index with the fewest mismatches for each sequence.
    """
    tfile = tempfile.NamedTemporaryFile()
    tfile.write(">ref1\n{}\n".format(ref1))
    tfile.write(">ref2\n{}\n".format(ref2))
    n_queries = 0
    for query in queries:
        tfile.write(">query{}\n{}\n".format(n_queries, query))
        n_queries += 1
    tfile.flush()
    
    outfile = tempfile.TemporaryFile()
    cmd = ["mafft", "--auto", tfile.name]
    subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE).communicate()
    outfile.seek(0)
    align = SeqIO.parse(outfile, "fasta")

    ref1 = next(align)
    ref2 = next(align)
    scores = [0, 0]
    ties = 0
    for query in align:
        mismatch1 = num_mismatches(ref1, query)
        mismatch2 = num_mismatches(ref2, query)
        if mismatch1 < mismatch2:
            scores[0] += 1
        elif mismatch1 > mismatch2:
            scores[1] += 1
        else:
            ties += 1
    if ties > 0:
        msg = "Discarding {} reads aligning equally to ref and alt".format(ties)
        warnings.warn(msg)

    return scores

def count_indels(samfile, reffile, chrom, pos, ref, alt, min_mapq=20):
    """
    Count occurences of the reference and alternate indel allele at a given
    position, by alignment score.
    """
    start = max(pos-100, 0)
    end = pos+100
    ref_seq = reffile.fetch(reference=chrom, start=start, end=end).decode("utf-8")
    alt_seq = ref_seq[:pos-start-1] + alt + ref_seq[pos-start-1+len(ref):]
    
    reads = samfile.fetch(chrom, pos, pos+len(ref))
    counts = [0, 0]
    
    reads = [r.seq for r in reads if r.mapq >= min_mapq]
    scores = multi_align(ref_seq, alt_seq, reads)
    
    return scores

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
                row["end"] = row["start"] + len(row["alt"]) - 1

                # get allele counts
                if len(row["ref"]) == len(row["alt"]) == 1: # snv
                    counts = count_bases(samfile, reffile, row["chr"], int(row["pos"]))
                    row["ref.count"] = counts[in_row["Reference_Allele"]]
                    row["alt.count"] = counts[in_row["Tumor_Seq_Allele1"]]
                    row["depth"] = sum(counts.values())
                else: # indel
                    counts = count_indels(samfile, reffile, row["chr"], int(row["pos"]), row["ref"], row["alt"])
                    row["ref.count"] = counts[0]
                    row["alt.count"] = counts[1]
                    row["depth"] = sum(counts)

                # info from MAF file
                row["gene"] = in_row["Hugo_Symbol"]
                row["class"] = in_row["Variant_Classification"]
                row["cosmic"] = ",".join(re.findall("COSM\d+", in_row["dbSNP_RS"]))
                row["rs"] = ",".join(re.findall("rs\d+", in_row["dbSNP_RS"]))
                row["esp"] = ",".join(re.findall("TMP_ESP_[\dXY]{1,2}_\d+", in_row["dbSNP_RS"]))
                if in_row["Protein_Change"] is None:
                    row["prot.change"] = ""
                else:
                    row["prot.change"] = in_row["Protein_Change"].replace("p.", "")

                writer.writerow(row)
