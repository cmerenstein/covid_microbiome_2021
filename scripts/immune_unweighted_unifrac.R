library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ade4)

setwd("../")

## Read in clinicaldata. 
clinical = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## UniFrac analysis makes most sense on OP swabs
OP = clinical[clinical$SampleType == "Oropharyngeal swab" & clinical$Study_group == "COVID",]

## read in UniFrac
unifrac = as.matrix(readRDS("data/from_scripts/unweighted_unifrac.rds"))

## read in the messi data and filter columns
immune = read.csv("data/from_scripts/immune/messi_immune_matchedOropharyngeal swab.csv", row.names = 1)
immune = immune[,c("SubjectID", "component1", "component2", "first_sample")]

## filter out when immune compoenent is NA
immune = immune[ !(is.na(immune$component1)),]

## get sample IDs for these samples
sample_mapping = merge(OP[, c("SubjectID", "SampleID", "Date")], immune, 
                    by.x = c("SubjectID", "Date"), by.y = c("SubjectID", "first_sample"))

## Get matching matricies for immune and unifrac data 
unifrac_matrix = unifrac[sample_mapping$SampleID, sample_mapping$SampleID]
rownames(unifrac_matrix) = sample_mapping$SubjectID
colnames(unifrac_matrix) = sample_mapping$SubjectID

## get as distance for Mantel test
unifrac_dist = as.dist(unifrac_matrix)
immune_dist = dist(immune[,c("component1", "component2")])

mantel.rtest(unifrac_dist, immune_dist, nrepet = 9999)


pdf(paste("figures/immune/OP_dist_comparisons_unweighted.pdf"))
data.frame(unifrac_dist = as.vector(unifrac_dist), immune_dist = as.vector(immune_dist)) %>%
    filter(unifrac_dist != 0) %>%
    ggplot(aes(x = unifrac_dist, y = immune_dist)) +
        theme_classic() +
        geom_point(size = 4) +
        theme(text = element_text(size = 20)) +
        stat_smooth(color = "black", method = "lm") +
        ylab("Immune Distance") +
        xlab("Unweighted UniFrac Distance")

dev.off()





