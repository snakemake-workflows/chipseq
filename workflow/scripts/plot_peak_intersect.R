#!/usr/bin/env Rscript

### source: https://github.com/nf-core/chipseq/blob/master/bin/plot_peak_intersect.r
### own changes and adjustments for snakemake-workflow chipseq are marked with "# AVI: " in comment

#MIT License
#
#Copyright (c) Philip Ewels
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(UpSetR)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="Path to tab-delimited file containing two columns i.e sample1&sample2&sample3 indicating intersect between samples <TAB> set size.", metavar="path"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="Path to output file with '.pdf' extension.", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)){
    print_help(opt_parser)
    stop("Input file must be supplied.", call.=FALSE)
}
if (is.null(opt$output_file)){
    print_help(opt_parser)
    stop("Output pdf file must be supplied.", call.=FALSE)
}

OutDir <- dirname(opt$output_file)
if (file.exists(OutDir) == FALSE) {
    dir.create(OutDir,recursive=TRUE)
}

################################################
################################################
## PLOT DATA                                  ##
################################################
################################################

comb.dat <- read.table(opt$input_file,sep="\t",header=FALSE)
comb.vec <- comb.dat[,2]
comb.vec <- setNames(comb.vec,comb.dat[,1])
sets <- sort(unique(unlist(strsplit(names(comb.vec),split='&'))), decreasing = TRUE)

nintersects = length(names(comb.vec))
if (nintersects > 70) {
    nintersects <- 70
    comb.vec <- sort(comb.vec, decreasing = TRUE)[1:70]
    sets <- sort(unique(unlist(strsplit(names(comb.vec),split='&'))), decreasing = TRUE)
}

pdf(opt$output_file,onefile=F,height=10,width=20)

upset(
    fromExpression(comb.vec),
    nsets = length(sets),
    nintersects = nintersects,
    sets = sets,
    keep.order = TRUE,
    sets.bar.color = "#56B4E9",
    point.size = 3,
    line.size = 1,
    mb.ratio = c(0.55, 0.45),
    order.by = "freq",
    number.angles = 30,
    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.2)
)

dev.off()

################################################
################################################
################################################
################################################
