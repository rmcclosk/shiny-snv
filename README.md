Dependencies
============

This visualization requires several pieces of software to be installed and
available.

 - R version 3.1.1 or higher
 - Python version 2.7.8 or higher
 - the Pysam python library
 - bedtools

In addition, the following R libraries are used.

 - shiny
 - ggvis
 - ggplot2

Running the Visualization
=========================

Copy "metadata_sample.tsv" to "metadata.tsv", and edit it for your data. The
"time.point" column is used to differentiate samples taken contemporaneously
vs. from different time points.  These should be numbered in chronological
order, but the actual numbers used are not important. The special time point 0
indicates a normal sample.

You also need to produce a file of copy number segments, which looks like
"segments_sample.tsv", called "segments.tsv". These can be made with the
segmentation tool of your choice (eg. DNAcopy, ExomeCNV, TITAN).

Finally, run "collect_vaf.py metadata.tsv /path/to/reference.fa > variants.tsv"
to produce a TSV file with variant allelic fractions.

The first time the visualization is run, it will take a long time, and may even
time out. This is because it is precomputing statistics about the variants and
copy number segments. Subsequent accesses will be much faster.

