Dependencies
============

This visualization requires several pieces of software.

- R version 3.1.1 or higher
- Python version 2.7.8 or higher
- the Pysam python library

In addition, several R libraries are used.

- shiny
- ggvis

Running the Visualization
=========================

Copy "metadata_sample.tsv" to "metadata.tsv", and edit it for your data. The
"time.point" column is used to differentiate samples taken contemporaneously
vs. from different time points.  These should be numbered in chronological
order, but the actual numbers used are not important. The special time point 0
indicates a normal sample.

Run "collect_vaf.py metadata.tsv /path/to/reference.fa > vaf.tsv". You will
also need to produce a file called "segments.tsv", which has the columns
(sample, chrom, start, end, copy.number, prevalence). This can be done with the
segmentation tool of your choice.
