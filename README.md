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

Edit "metadata.tsv" for your data. The "time.point" column is used to
differentiate samples taken contemporaneously vs. from different time points.
These should be numbered in chronological order, but the actual numbers used
are not important. The special time point 0 indicates a normal sample.
