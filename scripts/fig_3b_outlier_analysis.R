library(ggplot2)
library(dplyr)
library(tidyr)

setwd("../")

## Read in clinicaldata. 
clinical = read.csv("data/from_scripts/merged_clinical.csv", 
                                                    stringsAsFactors = F, row.names = 1)
## read counts post-contamination filter
counts = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", 
                                                    row.names = 1, stringsAsFactors = F)
## subsetted samples to >1000 reads for unifrac
counts = counts[rowSums(counts) > 1000,]
clinical = clinical[rownames(counts),]

## get just ETA and OP samples
clinical_filter = clinical[clinical$SampleType %in% c("Oropharyngeal swab", "Endotracheal aspirate") &
                                                                    clinical$Study_group == "COVID",]
## add subject-time pair so we can match OP and ETA later
clinical_filter$Subject_Sample = paste(clinical_filter$SubjectID, clinical_filter$visit_int, sep = "_")

## 
ETA = clinical_filter[clinical_filter$SampleType == "Endotracheal aspirate",]

## get first ETA sample
#ETA$first_visit = ave(ETA$visit_int, ETA$SubjectID, FUN = min)
#ETA_first = ETA[ ETA$visit_int == ETA$first_visit, ]

## filter to only these subject-sample points
ETA_OP = clinical_filter[clinical_filter$Subject_Sample %in% ETA$Subject_Sample, ]

## filter to only those we have pairs for
samples = ETA_OP %>% group_by(Subject_Sample) %>% summarize(n = n()) %>%
    ungroup() %>% as.data.frame()
rownames(samples) = samples$Subject_Sample
pairs = samples[samples$n == 2,]

ETA_OP = ETA_OP[ ETA_OP$Subject_Sample %in% pairs$Subject_Sample,]

## only look at relatively common genera
percent = counts / rowSums(counts)
percent = percent[,colMeans(percent) > 0.001]


## ------------ FIND OUTLIERS ------------------------------------
## use binomial model, using counts data
counts_df = data.frame(sample = rownames(ETA_OP), counts[rownames(ETA_OP), colnames(percent)])
counts_tidy = gather(counts_df, "genus", "counts", -sample)
counts_tidy$SampleType = ETA_OP[counts_tidy$sample, "SampleType"]
counts_tidy$Subject_Sample = ETA_OP[counts_tidy$sample, "Subject_Sample"]

counts_tidy_ETA = counts_tidy[counts_tidy$SampleType == "Endotracheal aspirate",]
counts_tidy_OP = counts_tidy[counts_tidy$SampleType == "Oropharyngeal swab",]

## merge back
ETA_OP_counts = full_join(counts_tidy_ETA, counts_tidy_OP, by = c("Subject_Sample", "genus")) %>%
    rename("counts_ETA" = counts.x) %>%
    rename("counts_OP" = counts.y)

## filter out taxa
ETA_OP_counts = ETA_OP_counts[ETA_OP_counts$counts_ETA > 0 & ETA_OP_counts$counts_OP > 0,]

## add sample totals
ETA_OP_counts$total_ETA = ave(ETA_OP_counts$counts_ETA, ETA_OP_counts$sample.x, FUN = sum)
ETA_OP_counts$total_OP = ave(ETA_OP_counts$counts_OP, ETA_OP_counts$sample.x, FUN = sum)

## turn to percentage
ETA_OP_counts$percent_ETA = ETA_OP_counts$counts_ETA / ETA_OP_counts$total_ETA
ETA_OP_counts$percent_OP = ETA_OP_counts$counts_OP / ETA_OP_counts$total_OP

## add WHO score
ETA_OP_counts$WHO = clinical[ETA_OP_counts$sample.x, "Max.WHO.score"]
ETA_OP_counts$WHO_Subject_Sample = paste("WHO:", ETA_OP_counts$WHO, ETA_OP_counts$Subject_Sample) 

## outlier label
ETA_OP_counts$logFC = log2(ETA_OP_counts$percent_ETA) - log2(ETA_OP_counts$percent_OP)
ETA_OP_counts$outlier = ETA_OP_counts$logFC > 3 & ETA_OP_counts$percent_ETA >= .5
ETA_OP_counts$outlier_genus = ifelse(ETA_OP_counts$outlier, ETA_OP_counts$genus, "")

n_outliers = group_by(ETA_OP_counts, Subject_Sample) %>% summarize(n_outlier = sum(outlier))

## ---------------- PLOT ----------------------------
pdf("figures/OP_ETA/outlier_plot_genera.pdf", height = 15, width = 10)
ggplot(ETA_OP_counts, aes(y = 0.001 + percent_ETA, x = 0.001 + percent_OP, fill = outlier_genus)) + 
    theme_classic() +
    geom_point(shape = 21, size = 3) + 
#    geom_text(aes(label = outlier), nudge_x = .1, nudge_y = .1, size = 2) +
    scale_y_log10() + 
    scale_x_log10() + 
    facet_wrap(~Subject_Sample, ncol = 3) +
    ylab("ETA") + 
    xlab("OP") +
#   coord_fixed() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_fill_manual(values = c( "grey", "#c39bd3", "#aed6f1", "forestgreen", "#82e0aa", "#980d00", "#FFD700", 
                                  "#FF4500", "#bb00ff"))
#ggsave(file = paste(folder, "figures/OP_ETA_distance/outlier_plot.svg", sep = ""), plot = plt, height = 20, width = 10)
dev.off()

