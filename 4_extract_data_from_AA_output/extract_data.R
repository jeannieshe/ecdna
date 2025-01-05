#!/usr/bin/env Rscript

# Part 0: Load in libraries ####
library(tidyverse)
library(gplots)
library(ggplot2)
library(gt)
library(readr)
library(optparse)
library(rstatix)
library(ggpubr)
library(dplyr)
library(ggbeeswarm)
library(RColorBrewer)
library(htmltools)

# Part 1: Parse inputs from command line ####

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="csv filename. find example in ./extract_data_example_sheet.csv", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-s", "--save"), type="boolean", default=NULL,
              help="save plots as .png", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one input file must be supplied.n", call.=FALSE)
}

# Part 2: Processing a file listing multiple samples ####
#(find example in ./extract_data_example_sheet.csv)
opt$file <- "/home/jeannie/group/jeannie/extract_data_08082024_run.csv"

## Part 2a: Create a table listing all ecDNA from samples ####
all_plots <- NULL
sample_names <- c()
if (!is.null(opt$file)) {
  amplicon_data <- data.frame()
  all_clean_data <- data.frame()
  lines <- read_lines(opt$file)
  
  for (line in lines) {
    #fields <- strsplit(line, ",")[[1]]
    foo <- str_split(line, "/")
    sample_name <- foo[[1]][length(foo[[1]])]
    sample_name <- substring(sample_name, 1, nchar(sample_name) - 7)
    sample_names <- c(sample_names, sample_name)
    table_dir <- paste0(sample_name, "_classification/", sample_name, "_result_table.tsv")
    full_dir <- paste0(line, "/", table_dir)
    result_table <- read.delim(full_dir, sep="\t")
    amplicon_data <- rbind(amplicon_data, result_table)
    
    clean_table <- result_table[c("Sample.name", "AA.amplicon.number", "Classification", "Oncogenes", "All.genes", "Feature.median.copy.number", 
                                  "Location", "Captured.interval.length", "ecDNA.context")]
    all_clean_data <- rbind(all_clean_data, clean_table)
  }
  
  ecDNA_table <- all_clean_data %>% filter(Classification == "ecDNA")
  ecDNA_table = subset(ecDNA_table, select = -c(Classification) )
}
#amplicon_data <- amplicon_data[(!is.na(amplicon_data$Classification)),] #the NA indicates no amplicons were identified in the sample
ecDNA_table$Amplicon.ID <- sprintf("%s_%s", ecDNA_table$Sample.name, ecDNA_table$AA.amplicon.number) 

#create a legible version and save it
ecDNA_table_clean <- select(ecDNA_table, -c(AA.amplicon.number))
colnames(ecDNA_table_clean) <- c("Sample name", "Oncogenes", "Median copy number", "Location", "Size", "ecDNA Context",
                                 "Amplicon ID", "Time to progression", "Post progression survival", "Platinum free interval",
                                 "Op1 Dissemination Index", "Op Dissemination Difference", "Time sampled", "Contains oncogenes", "Primary Therapy Outcome")
ecDNA_table_clean <- ecDNA_table_clean[, c("Sample name", "Amplicon ID", "Oncogenes", "Median copy number", "Location", "Size", "ecDNA Context",
                                           "Time to progression", "Post progression survival", "Platinum free interval",
                                           "Op1 Dissemination Index", "Op Dissemination Difference", "Time sampled", "Primary Therapy Outcome", "Contains oncogenes")]
ecDNA_table_gt <- gt(ecDNA_table_clean)
ecDNA_table_gt <- ecDNA_table_gt %>%
  data_color(
    columns = c("Median copy number", "Size"),
    colors = col_numeric(
      palette = brewer.pal(9, "GnBu"), # Use a Brewer palette
      domain = NULL, # Automatically detect the range of values
      na.color = "white" # Set NA values to white
    )
  ) %>%
  tab_options(
    table.font.size = px(12) # Set font size to 10 pixels
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "left", color = "grey", weight = px(1) # Vertical grid lines
    ),
    locations = cells_body(
      columns = everything()
    )
  )
print(ecDNA_table_gt)
html_file <- "ecDNA_table.html"
save_html(ecDNA_table_gt, file = html_file)

## Part 2b: Create table of oncogenes expressed in all amplicons. ####
#extract list of unique oncogenes from all_clean_data
oncogenes <- all_clean_data %>% pull(Oncogenes)
oncogenes <- oncogenes[oncogenes != "[]"]
cleaned_onco <- gsub("\\[|\\]|'", "", oncogenes)
split_onco <- unlist(strsplit(cleaned_onco, ",\\s*")) #splits the string on commas and whitespace, flattens
oncogenes <- unique(split_onco)

#create a df that has oncogenes on rows, samples on columns, and the kind of amplicon as datavalue
onco_table <- data.frame(matrix(NA, nrow = length(unique(all_clean_data$Sample.name)), ncol = length(oncogenes)))
colnames(onco_table) <- oncogenes
rownames(onco_table) <- unique(all_clean_data$Sample.name)
unique_samples <- unique(all_clean_data$Sample.name)

for (ii in seq_along(unique_samples)) {
  sample_name <- unique_samples[ii]
  curr_data <- all_clean_data[all_clean_data$Sample.name == sample_name, ]
  
  for (j in seq_len(nrow(curr_data))) {
    onco <- curr_data[j, ]$Oncogenes
    cleaned_onco <- gsub("\\[|\\]|'", "", onco)
    split_onco <- unlist(strsplit(cleaned_onco, ",\\s*")) # splits the string on commas and whitespace, flattens
    
    if (length(split_onco) > 0) {
      type <- curr_data[j, ]$Classification
      for (oncogene in split_onco) {
        onco_table[sample_name, oncogene] <- type
      }
    }
  }
}
onco_table <- t(onco_table)
onco_table[is.na(onco_table)] <- ""
onco_table <- as.data.frame(onco_table)

amplicon_factors <- c("BFB", "ecDNA", "Linear", "Complex-non-cyclic")
colors <- brewer.pal(4, "Pastel1")
display.brewer.pal(4, "Pastel1")
gt_table <- gt(onco_table, rownames_to_stub = TRUE)
gt_table <- gt_table %>%
  data_color(
    columns = everything(),
    colors = scales::col_factor(
      palette = colors,
      domain = amplicon_factors,
      na.color = "white"
    )
  ) %>%
  tab_options(
    table.font.size = px(12) # Set font size to 10 pixels
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "left", color = "grey", weight = px(1) # Vertical grid lines
    ),
    locations = cells_body(
      columns = everything()
    )
  )

print(gt_table)
html_file <- "onco_table.html"
save_html(gt_table, file = html_file)

##Part 2c: Create table of all genes ####
allgenes <- all_clean_data %>% pull(All.genes)
allgenes <- allgenes[allgenes != "[]"]
cleaned_all <- gsub("\\[|\\]|'", "", allgenes)
split_all <- unlist(strsplit(cleaned_all, ",\\s*")) #splits the string on commas and whitespace, flattens
allgenes <- unique(split_all)

all_genes_table <- data.frame(matrix(NA, nrow = length(unique(all_clean_data$Sample.name)), ncol = length(allgenes)))
colnames(all_genes_table) <- allgenes
rownames(all_genes_table) <- unique(all_clean_data$Sample.name)
unique_samples <- unique(all_clean_data$Sample.name)

for (ii in seq_along(unique_samples)) {
  sample_name <- unique_samples[ii]
  curr_data <- all_clean_data[all_clean_data$Sample.name == sample_name, ]
  
  for (j in seq_len(nrow(curr_data))) {
    genes <- curr_data[j, ]$All.genes
    cleaned_genes <- gsub("\\[|\\]|'", "", genes)
    split_genes <- unlist(strsplit(cleaned_genes, ",\\s*")) # splits the string on commas and whitespace, flattens
    
    if (length(split_genes) > 0) {
      type <- curr_data[j, ]$Classification
      for (gene in split_genes) {
        all_genes_table[sample_name, gene] <- type
      }
    }
  }
}

all_genes_table <- t(all_genes_table)
all_genes_table[is.na(all_genes_table)] <- ""
all_genes_table <- as.data.frame(all_genes_table)

amplicon_factors <- c("BFB", "ecDNA", "Linear", "Complex-non-cyclic")
colors <- brewer.pal(4, "Pastel1")
display.brewer.pal(4, "Pastel1")
gt_table <- gt(all_genes_table, rownames_to_stub = TRUE)
gt_table <- gt_table %>%
  data_color(
    columns = everything(),
    colors = scales::col_factor(
      palette = colors,
      domain = amplicon_factors,
      na.color = "white"
    )
  ) %>%
  tab_options(
    table.font.size = px(10) # Set font size to 10 pixels
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "left", color = "grey", weight = px(1) # Vertical grid lines
    ),
    locations = cells_body(
      columns = everything()
    )
  )
print(gt_table)
html_file <- "all_genes_table.html"
save_html(gt_table, file = html_file)


# Part 3: Statistical analysis of the data ####
#populate sample_data
variables <- c("Number.ecDNA", "Number.ecDNA.onco", "Time.to.progression", "Post.progression.survival", "Platinum.free.interval", "Op1.dissemination.idx", "Op.dissemination.difference")
sample_data <- data.frame(matrix(NA, nrow = length(unique(amplicon_data$Sample.name)), ncol = length(variables)))
rownames(sample_data) <- unique(amplicon_data$Sample.name)
colnames(sample_data) <- variables

sample_data$Time.to.progression <- c(NA, NA, NA, 236, 236, 236, 236, 670, 288, 288, 288, 300, 382, 
                                       382, 271, 271, 271, 271, 347, 288, 288, 288, 308, 534, 534, 
                                       845, 496, 496, NA, 377, 377, NA, 300)
sample_data$Post.progression.survival <- c(NA, NA, NA, 706, 706, 706, 706, 371, 151, 151, 151, 687,
                                             318, 318, 769, 769, 769, 769, NA, 287, 287, 287, 40, NA,
                                             NA, NA, NA, NA, NA, 299, 299, NA, 687)
sample_data$Platinum.free.interval <- c(NA, NA, NA, 81, 81, 81, 81, 511, 71, 71, 71, 126, 91, 91, 
                                          92, 92, 92, 92, 195, 83, 83, 83, 176, 393, 393, 692, 341, 
                                          341, NA, 230, 230, NA, 126)
sample_data$Op1.dissemination.idx <- c(NA, NA, NA, 9, 9, 9, 9, 15, 10, 10, 10, 7, 8, 8, 16, 
                                         16, 16, 16, 15, 15, 15, 15, 13, 11, 11, 12, 8, 8, NA, 
                                         10, 10, NA, 7)
sample_data$Op.dissemination.difference <- c(NA, NA, NA, -7, -7, -7, -7, -5, NA, NA, NA, 0, NA, 
                                               NA, -14, -14, -14, -14, NA, -7, -7, -7, NA, -5, 
                                               -5, NA, NA, NA, NA, -5, -5, NA, 0)  
sample_data$Outcome <- c(NA, NA, NA, "Complete", "Complete", "Complete", "Complete", "Partial", 
                         "Partial", "Partial", "Partial", "Complete", "Stable", "Stable", "Complete", "Complete", 
                         "Complete", "Complete", "Partial", "Partial", "Partial", "Partial", "Complete", 
                         "Complete", "Complete", "Partial", "Complete", "Complete", NA, 
                         "Partial", "Partial", NA, "Complete")

## Part 3a: Plot differences between H and D samples ####
# plot ecDNA count
ecDNA_counts <- table(ecDNA_table$Sample.name)
sample_data$Number.ecDNA <- 0
matching_indices <- match(rownames(sample_data), names(ecDNA_counts))
sample_data$Number.ecDNA[!is.na(matching_indices)] <- ecDNA_counts[matching_indices[!is.na(matching_indices)]]
sample_data$Sample.group <- ifelse(substr(rownames(sample_data), 1, 1) == "D", "D", "H")
sample_data$Time.sampled <- ifelse(substr(rownames(sample_data), 6, 6) == "p", "Primary", 
                                     ifelse(substr(rownames(sample_data), 6, 6) == "r", "Relapse", "Interval")) 

num_ecDNA <- ggplot(data = sample_data, aes(x = Sample.group, y = Number.ecDNA, fill = Sample.group)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Number of ecDNA identified in sample groups") +
  xlab("Sample Group") +
  ylab("Number of ecDNA identified by AmpliconArchitect") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("H", "D")), method = "wilcox.test", label = "p.signif")
ggsave("num_ecDNA.png", num_ecDNA, width = 6, height = 7, units = "in")

# plot ecDNA with oncogenes count
ecDNA_with_onco_counts <- ecDNA_table %>%
  filter(nchar(ecDNA_table$Oncogenes) > 2) %>% #if the list is more than "[]"
  count(Sample.name)
sample_data$Number.ecDNA.onco <- 0
matching_indices <- match(rownames(sample_data), ecDNA_with_onco_counts$Sample.name)
sample_data$Number.ecDNA.onco[!is.na(matching_indices)] <- ecDNA_with_onco_counts$n[matching_indices[!is.na(matching_indices)]]

num_ecDNA_with_onco <- ggplot(data = sample_data, aes(x = Sample.group, y = Number.ecDNA.onco, fill = Sample.group)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Number of ecDNA with oncogenes identified in sample groups") +
  xlab("Sample Group") +
  ylab("Number of ecDNA with oncogenes identified by AmpliconArchitect") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("H", "D")), method = "wilcox.test", label = "p.signif")
#facet_wrap(~variable)
ggsave("num_ecDNA_with_onco.png", num_ecDNA_with_onco, width = 6, height = 7, units = "in")

# plot ecDNA mean size (bp)
sample_data$Mean.ecDNA.size <- 0
ecDNA_mean_size <- aggregate(Captured.interval.length ~ Sample.name, data = ecDNA_table, FUN = mean)
sample_data$Mean.ecDNA.size[!is.na(matching_indices)] <- ecDNA_mean_size$Captured.interval.length[matching_indices[!is.na(matching_indices)]]

mean_ecDNA <- ggplot(data = sample_data, aes(x = Sample.group, y = Mean.ecDNA.size, fill = Sample.group)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Mean ecDNA size identified in sample groups") +
  xlab("Sample Group") +
  ylab("Mean ecDNA size (bp)") +
  theme_gray() +
  stat_compare_means(comparisons = list(c("H", "D")), method = "wilcox.test", label = "p.signif")
ggsave("mean_ecDNA.png", mean_ecDNA, width = 6, height = 7, units = "in")

# all_plots <- c(all_plots, num_ecDNA, num_ecDNA_with_onco, mean_ecDNA)


## Part 3b: Plot changes in copy number. ####
#first, across primary/relapse samples
#this is only the case for patient H092

foo <- c(9, 15)
bar <- c(9.086900, 9.919485)
cn_across_disease <- data.frame(foo, bar)
rownames(cn_across_disease) <- c("H092_pOme1-P9", "H092_rAsc-P15")
colnames(cn_across_disease) <- c("Passage", "Copy.Number")

across_disease <- ggplot(data = cn_across_disease, aes(x=Passage, y=Copy.Number)) + 
  geom_point() + 
  geom_line(aes(color = 'Patient H092'), linewidth=1) +
  geom_text(aes(label = c("Primary (pOme1-P9)", "Relapse (rAsc-P15)")), nudge_x = -0.2, nudge_y = -0.1, check_overlap=TRUE) + 
  ggtitle("Change in copy number of ecDNA over disease progression") +
  xlab("Passage Number") +
  ylab("Copy Number") +
  theme_gray() +
  xlim(0, 18) +
  ylim(7, 10) +
  labs(color = "Patient")

ggsave("CN_across_disease_progression.png", across_disease, width = 6, height = 7, units = "in")

#second, across passages
cn <- c(9.086900, 9.043967, 15.487663, 14.784765, 
        32.014240, 0, 9.006012, 12.606355, 
        0, 7.838868, 4.911985, 0, 
        5.475675, 0, 0, 5.203764)
passages <- c(9, 20, 9, 20, 
              9, 20, 6, 20,
              6, 20, 4, 15,
              7, 15, 7, 15)
ecDNA_id <- c(1, 1, 2, 2, 
              3, 3, 4, 4,
              5, 5, 6, 6,
              7, 7, 8, 8) #for plotting purposes
cn_across_passages <- data.frame(passages, cn, ecDNA_id)
rownames(cn_across_passages) <- c("H092_pOme1-P9_ecDNA1", "H092_pOme1-P20_ecDNA1", "H092_pOme1-P9_ecDNA2", "H092_pOme1-P20_ecDNA2", 
                                  "H092_pOme1-P6_ecDNA3", "H092_pOme1-P20_ecDNA3", "H180_p2Ome-P6_ecDNA1", "H180_p2Ome-P20_ecDNA1", 
                                  "H180_p2Ome-P6_ecDNA2", "H180_p2Ome-P20_ecDNA2", "H143_rAsc-P4_ecDNA1", "H143_rAsc-P15_ecDNA1",
                                  "H215_pAsc-P7_ecDNA1", "H215_pAsc-P15_ecDNA1", "H215_pAsc-P7_ecDNA2", "H215_pAsc-P15_ecDNA2")
colnames(cn_across_passages) <- c("Passage", "Copy.Number", "ecDNA.ID")
across_passages <- ggplot(data = cn_across_passages, aes(x=Passage, y=Copy.Number)) + 
  geom_point() + 
  #geom_text(label=rownames(cn_across_passages)) + 
  geom_line(data = filter(cn_across_passages, ecDNA_id == 1), aes(color = 'H092_pOme1_ecDNA1'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 2), aes(color = 'H092_pOme1_ecDNA2'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 3), aes(color = 'H092_pOme1_ecDNA3'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 4), aes(color = 'H180_p2Ome_ecDNA1'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 5), aes(color = 'H180_p2Ome_ecDNA2'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 6), aes(color = 'H143_rAsc_ecDNA1'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 7), aes(color = 'H215_pAsc_ecDNA1'), linewidth=1) +
  geom_line(data = filter(cn_across_passages, ecDNA_id == 8), aes(color = 'H215_pAsc_ecDNA2'), linewidth=1) +
  #scale_color_manual(values=c('orange')) +
  ggtitle("Change in copy number of ecDNA over passages") +
  xlab("Passage Number") +
  ylab("Copy Number") +
  theme_gray() +
  xlim(0, 21) +
  ylim(0, 33) +
  labs(color = "ecDNA")

ggsave("CN_across_passages.png", across_passages, width = 6, height = 7, units = "in")

## Part 3c: Plot changes per amplicon across the clinical statistics. ####
#populate the ecDNA_table
variables <- c("Time.to.progression", "Post.progression.survival", "Platinum.free.interval", 
               "Op1.dissemination.idx", "Op.dissemination.difference", "Time.sampled", "Outcome")
for (var in variables) {
  ecDNA_table[[var]] <- NA # Initialize the column with NA
  for (ii in 1:nrow(ecDNA_table)) {
    rowname <- ecDNA_table$Sample.name[ii]
    ecDNA_table[[var]][ii] <- sample_data[rownames(sample_data) == rowname, var]
  }
}
ecDNA_table$Contains.onco <- ifelse(nchar(ecDNA_table$Oncogenes) > 2, TRUE, FALSE)

#produce three plots for each variable
variables <- c("Time.to.progression", "Post.progression.survival", "Platinum.free.interval", 
               "Op1.dissemination.idx", "Op.dissemination.difference")
axes_labels <- c("Time to progression (days)", "Post progression survival (days)", "Platinum free interval (days)",
                 "Operation 1 Dissemination Index", "Dissemination Index Difference")

for (i in seq_along(variables)) {
  var <- variables[i]
  axis_label <- axes_labels[i]
  
  cnplot <- ggplot(data = ecDNA_table, aes_string(x = var, y = "Feature.median.copy.number")) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 1) +
    facet_wrap(~ Time.sampled) +
    ggtitle(sprintf("ecDNA Copy Number vs. %s", axis_label)) +
    xlab(axis_label) +
    ylab("Copy Number") +
    theme_gray()
  
  cn_oncoplot <- ggplot(data = filter(ecDNA_table, Contains.onco == TRUE), aes_string(x = var, y = "Feature.median.copy.number")) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 1) +
    facet_wrap(~ Time.sampled) +
    ggtitle(sprintf("ecDNA with Oncogene Copy Number vs. %s", axis_label)) +
    xlab(axis_label) +
    ylab("Copy Number") +
    theme_gray()
  
  stat_boxplot <- ggplot(data = ecDNA_table, aes_string(x = "Contains.onco", y = var, fill = "Contains.onco")) +
    geom_boxplot(alpha = 0.5) +
    ggtitle(sprintf("%s vs. ecDNA Content", axis_label)) +
    xlab("ecDNA Content") +
    ylab(axis_label) +
    theme_gray() + 
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), method = "wilcox.test", label = "p.signif") +
    scale_x_discrete(labels = c("TRUE" = "Oncogene(s) Present", "FALSE" = "Oncogene(s) Absent")) +
    scale_fill_manual(values = c("TRUE" = "blue2", "FALSE" = "deeppink2"), 
                      labels = c("TRUE" = "Oncogene(s) Present", "FALSE" = "Oncogene(s) Absent")) +
    labs(fill = "Oncogenes Present in ecDNA")
  
  ggsave(sprintf("%svsCN.png", var), cnplot, width = 10, height = 6, units = "in")
  ggsave(sprintf("%svsCNonco.png", var), cn_oncoplot, width = 10, height = 6, units = "in")
  ggsave(sprintf("%svsCNboxplot.png", var), stat_boxplot, width = 7, height = 6, units = "in")
  
}

sample_data$Contains.onco <- ifelse(sample_data$Number.ecDNA.onco > 0, TRUE, FALSE)
sample_data$Present.ecDNA <- ifelse(sample_data$Number.ecDNA > 0, TRUE, FALSE)
sample_data$Relapse.ecDNA <- ifelse(sample_data$Number.ecDNA > 0 & 
                                    sapply(strsplit(rownames(sample_data), "_"), function(x) substr(x[2], 1, 1)) == "r",
                                    TRUE, FALSE) #since its vectorized, need to use one & and sapply

#produce specific boxplots comparing samples with ecDNA under various filters
stat_boxplot <- ggplot(data = sample_data, aes(x = Contains.onco, y = Time.to.progression, color = Contains.onco)) +
  geom_beeswarm(alpha = 1) +
  ggtitle("Time to progression vs. ecDNA Content of Each Sample") +
  xlab("ecDNA Content") +
  ylab("Time to progression (days)") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("TRUE" = "Oncogene(s) Present", "FALSE" = "Oncogene(s) Absent")) +
  scale_color_manual(values = c("TRUE" = "blue2", "FALSE" = "deeppink2"), 
                    labels = c("TRUE" = "Oncogene(s) Present", "FALSE" = "Oncogene(s) Absent")) +
  labs(color = "ecDNA Content")

stat_boxplot <- ggplot(data = filter(sample_data, Time.sampled == "Relapse") , aes(x = Relapse.ecDNA, y = Time.to.progression, fill = Relapse.ecDNA)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Time to progression vs. ecDNA Presence in Relapse Samples") +
  xlab("ecDNA Presence in Relapse Samples") +
  ylab("Time to progression (days)") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("TRUE" = "ecDNA Present", "FALSE" = "ecDNA Absent")) +
  scale_fill_manual(values = c("TRUE" = "blue2", "FALSE" = "deeppink2"), 
                     labels = c("TRUE" = "ecDNA Present", "FALSE" = "ecDNA Absent")) +
  labs(fill = "ecDNA Presence in Relapse Samples")

stat_boxplot <- ggplot(data = sample_data, aes(x = Present.ecDNA, y = Time.to.progression, fill = Present.ecDNA)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Time to progression vs. ecDNA Presence in Samples") +
  xlab("ecDNA Presence") +
  ylab("Time to progression (days)") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("TRUE" = "ecDNA Present", "FALSE" = "ecDNA Absent")) +
  scale_fill_manual(values = c("TRUE" = "blue2", "FALSE" = "deeppink2"), 
                    labels = c("TRUE" = "ecDNA Present", "FALSE" = "ecDNA Absent")) +
  labs(fill = "ecDNA Presence in Samples")

#plot mean copy number of amplicons vs TTP
amplicon_data$Patient.ID <- ""
for (ii in 1:nrow(amplicon_data)) {
  amplicon_data$Patient.ID[ii] <- unlist(strsplit(amplicon_data$Sample.name[ii], "-"))[1]
}
variables <- c("Mean.CN")
clinical_data <- data.frame(matrix(NA, nrow = length(unique(amplicon_data$Patient.ID)), ncol = length(variables)))
rownames(clinical_data) <- unique(amplicon_data$Patient.ID)
colnames(clinical_data) <- variables
clinical_data$Mean.CN <- 0
mean_CN_per_patient <- aggregate(Feature.median.copy.number ~ Patient.ID, data = amplicon_data, FUN = mean)
matching_indices <- match(rownames(clinical_data), mean_CN_per_patient$Patient.ID)
clinical_data$Mean.CN[!is.na(matching_indices)] <- mean_CN_per_patient$Feature.median.copy.number[matching_indices[!is.na(matching_indices)]]

clinical_data$Time.to.progression <- NA
var <- c("Time.to.progression")
for (ii in 1:nrow(clinical_data)) {
  rowname <- rownames(clinical_data)[ii]
  matching_rows <- which(sapply(rownames(clinical_data), function(x) substr(unlist(strsplit(x, "_"))[2], 1, 1) == "r"))
  clinical_data$Time.to.progression[ii] <- sample_data[matching_rows[1], "Time.to.progression"]
}

meancnplot <- ggplot(data = clinical_data, aes(x = Time.to.progression, y = Mean.CN)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 1) +
  ggtitle("Mean Copy Number of All Amplicons per Patient vs. Time to Progression") +
  xlab("Time to Progression (days)") +
  ylab("Mean Copy Number") +
  theme_gray()

# plot mean copy number vs outcome
clinical_data$Outcome <- c(NA, NA, NA, "Complete", "Complete", "Complete", "Partial", 
                           "Partial", "Partial", "Complete", "Stable", "Stable", "Complete", "Complete", 
                           "Partial", "Partial", "Partial", "Complete", 
                           "Complete", "Complete", "Partial", "Complete", NA, 
                           "Partial", NA)
clinical_data$Time.sampled <- ifelse(substr(rownames(clinical_data), 6, 6) == "p", "Primary", 
                                     ifelse(substr(rownames(clinical_data), 6, 6) == "r", "Relapse", "Interval")) 

stat_boxplot <- ggplot(data = filter(clinical_data, !is.na(Outcome)), aes(x = Outcome, y = Mean.CN, fill = Outcome)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Mean CN vs. Patient Primary Therapy Response") +
  xlab("Patient Primary Therapy Response") +
  ylab("Mean Copy Number of All Amplicons per Sample") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("Partial", "Complete"), c("Complete", "Stable"), c("Partial", "Stable")), method = "wilcox.test", label = "p.signif") +
  facet_wrap(~factor(Time.sampled, c("Primary", "Interval", "Relapse"))) +
  labs(fill = "Patient Primary Therapy Response")

clinical_data$Mean.CN.onco <- 0
mean_CN_onco_per_patient <- aggregate(Feature.median.copy.number ~ Patient.ID, data = filter(amplicon_data, nchar(amplicon_data$Oncogenes) > 2), FUN = mean)
matching_indices <- match(rownames(clinical_data), mean_CN_onco_per_patient$Patient.ID)
clinical_data$Mean.CN.onco[!is.na(matching_indices)] <- mean_CN_onco_per_patient$Feature.median.copy.number[matching_indices[!is.na(matching_indices)]]

stat_boxplot <- ggplot(data = filter(clinical_data, !is.na(Outcome)), aes(x = Outcome, y = Mean.CN.onco, fill = Outcome)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Mean CN in Oncogene Containing ecDNA vs. Patient Primary Therapy Response") +
  xlab("Patient Primary Therapy Response") +
  ylab("Mean Copy Number of All Oncogene Containing ecDNA per Sample") +
  theme_gray() + 
  stat_compare_means(comparisons = list(c("Partial", "Complete"), c("Complete", "Stable"), c("Partial", "Stable")), method = "wilcox.test", label = "p.signif") +
  facet_wrap(~factor(Time.sampled, c("Primary", "Interval", "Relapse"))) +
  labs(fill = "Patient Primary Therapy Response")

#Part 4: Plot the limits of AmpliconArchitect ####
# ecDNA_size <- ggplot(ecDNA_table, aes(x=Captured.interval.length)) +
#   geom_histogram(color="black", fill="gray", bins = 5) +
#   ggtitle("Size of Identified ecDNA") +
#   xlab("ecDNA Size (bp)") +
#   ylab("Frequency") +
#   scale_x_continuous(name = "ecDNA Size (bp)", labels = label_number())
png("hist_ecDNAsize.png", height = 600, width = 700, units = "px")
hist_ecDNAsize <- hist(ecDNA_table$Captured.interval.length,
     xlab = "ecDNA Size (bp)",          
     main = "Size of Identified ecDNA",  
     col = "gray",                       
     border = "black") 
dev.off()

png("hist_ecDNACN.png", height = 600, width = 700, units = "px")
hist_ecDNACN <- hist(ecDNA_table$Feature.median.copy.number,
     xlab = "ecDNA Copy Number",          
     main = "Copy Number of Identified ecDNA",  
     col = "gray",                       
     border = "black")
dev.off()

print(mean(ecDNA_table$Captured.interval.length))
print(sum(ecDNA_table$Contains.onco))

#Part 5: Display subsets of tables for viewing purposes ####
selected_df <- amplicon_data %>%
  filter(grepl("180", Sample.name)) %>%
  select(Sample.name, Oncogenes, All.genes, Feature.median.copy.number, Location, Captured.interval.length) #%>%
  column_to_rownames(var = Sample.name)
  
selected_df <- amplicon_data %>%
  filter(grepl("180", Sample.name)) %>%
  select(Number.ecDNA) #%>%
  column_to_rownames(var = rowname)

#Part 6: Save all plots. ####

# save all plots as a single pdf **need to figure this out
  
pdf("allplots.pdf", onefile = TRUE)
for (plot in all_plots) {
    print(plot) #be sure that each plot is a ggplot object
}
dev.off()

ecDNA_table_rds <- saveRDS(ecDNA_table, file = "/isilonsund/NextGenSeqData/project-data/groups/wennerberglab/jeannie/ecDNA_table_data.rds")
onco_table_rds <- saveRDS(onco_table, file = "/isilonsund/NextGenSeqData/project-data/groups/wennerberglab/jeannie/onco_table.rds")
sample_data_rds <- saveRDS(sample_data, file = "/isilonsund/NextGenSeqData/project-data/groups/wennerberglab/jeannie/sample_data.rds")
amplicon_data_rds <- saveRDS(amplicon_data, file = "/isilonsund/NextGenSeqData/project-data/groups/wennerberglab/jeannie/amplicon_data.rds")
clinical_data_rds <- saveRDS(clinical_data, file = "/isilonsund/NextGenSeqData/project-data/groups/wennerberglab/jeannie/clinical_data.rds")

clinical_data_gt <- gt(clinical_data, rownames_to_stub = TRUE)
clinical_data_gt <- clinical_data_gt %>%
  tab_options(
    table.font.size = px(12) # Set font size to 10 pixels
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "left", color = "grey", weight = px(1) # Vertical grid lines
    ),
    locations = cells_body(
      columns = everything()
    )
  )
print(clinical_data_gt)
html_file <- "clinical_data.html"
save_html(clinical_data_gt, file = html_file)
