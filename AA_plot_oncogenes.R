#this script plots oncogenes against samples and specifies in what kind of amplicon the oncogene was amplified
#load in libraries
library(tidyverse)
library(gridExtra)
library(gplots)
library(gt)

#load in datatable from AA
H143rAscCLP4 <- read.delim("/home/jeannie/group/ampsuite/H143rAscCLP4/H143rAscCLP4_output/H143rAscCLP4_classification/H143rAscCLP4_result_table.tsv", sep="\t")
H143rAscCLP15 <- read.delim("/home/jeannie/group/ampsuite/H143-rAsc-P15/H143rAscP15_output/H143rAscP15_classification/H143rAscP15_result_table.tsv", sep="\t")
H143iOmeCLP5 <- read.delim("/home/jeannie/group/ampsuite/H143_iOme_CL_P5/H143iOmeCLP5_output/H143iOmeCLP5_classification/H143iOmeCLP5_result_table.tsv", sep="\t")
H143iOmeCLP16 <- read.delim("/home/jeannie/group/ampsuite/H143iOmeCLP16/H143iOmeCLP16_output/H143iOmeCLP16_classification/H143iOmeCLP16_result_table.tsv", sep="\t")

samples <- list(H143rAscCLP4, H143rAscCLP15, H143iOmeCLP5, H143iOmeCLP16)
sample_names <- c('H143rAscCLP4', 'H143rAscCLP15', 'H143iOmeCLP5', 'H143iOmeCLP16')
oncogenes <- c()
amplicon_levels <- c("ecDNA", "BFB", "Linear")

#create list of unique oncogenes
for (ii in 1:length(samples)) {
  oncogenes <- c(oncogenes, samples[[ii]]$Oncogenes)
}
cleaned_onco <- gsub("\\[|\\]|'", "", oncogenes)
split_onco <- unlist(strsplit(cleaned_onco, ",\\s*")) #splits the string on commas and whitespace, flattens
oncogenes <- unique(split_onco)

#make a dataframe to plot
#we want to make a df that has oncogenes on rows, samples on columns, and then the kind of amplicon as datavalue
plot_data <- data.frame(matrix(NA, nrow = length(samples), ncol = length(oncogenes)))
colnames(plot_data) <- oncogenes
rownames(plot_data) <- sample_names
for (ii in 1:length(samples)) {
  for (j in 1:nrow(samples[[ii]])) {
    onco <- samples[[ii]][j,]$Oncogenes
    cleaned_onco <- gsub("\\[|\\]|'", "", onco)
    split_onco <- unlist(strsplit(cleaned_onco, ",\\s*")) #splits the string on commas and whitespace, flattens
    onco <- split_onco #now this is a vector of the specific oncogenes expressed by this amplicon
    if (!is.null(onco)) {
      type <- samples[[ii]][j,]$Classification
      plot_data[sample_names[[ii]], onco] <- c(type)
      }
  }
}

amplicon_factors <- factor(plot_data[,], levels = amplicon_levels)
amplicon_factors <- addNA(amplicon_factors)

bar <- gt(plot_data, rownames_to_stub = T)
foo <- data_color(bar, palette = c("blue", "red", "green") , na_color = "white", levels = amplicon_factors)

show_col(col_factor("RdYlBu", levels = amplicon_factors)(amplicon_factors[1:3]))


ss <- tableGrob(plot_data)


grid.arrange(ss)


# script to get data from output folder or output file
# work on getting the DNA seq
#use pysam, .bam.bai files
#construct the sequence using .bedf ile and ref genome
#try to find it in the .fastq
