#!/usr/bin/env Rscript
library(ggplot2)

sample.data <- read.table("metadata.csv", header=T)
have.data <- file.exists(file.path("data", sample.data$tumor.sample))
sample.data <- sample.data[have.data,]
sample.data$date <- as.Date(sample.data$date, "%m/%d/%Y")

seg.files <- Sys.glob(file.path("data", sample.data$tumor.sample, "*segs.txt"))
all.segs <- do.call(rbind, lapply(seg.files, read.table, header=T))
all.segs <- merge(all.segs, sample.data, by.x=c("Sample"), 
                  by.y=c("tumor.sample"))

#by(all.segs, all.segs$patient, function (d) {
#    usamp <- unique(d$Sample)
#    sample.num <- sapply(d$Sample, function (s) {
#        which(usamp==s) - as.integer(length(usamp)/2)
#    })
#    d$Copy_Number <- 0.1*sample.num + d$Copy_Number
#
#    ggplot(d, aes(x=Start_Position.bp., y=Copy_Number)) +
#        geom_segment(aes(xend=End_Position.bp., yend=Copy_Number, colour=Sample)) +
#        facet_wrap(~Chromosome, scales="free_x") +
#        scale_x_discrete(breaks=NULL, name="") +
#        scale_y_discrete(breaks=0:5, name="copy number") 
#    outfile <- file.path("plots", paste0(d[1,"patient"], ".pdf"))
#    ggsave(file=outfile, width=10)
#})
