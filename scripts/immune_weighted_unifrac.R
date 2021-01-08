library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ade4)

## Test the UniFrac distance between ETA OP samples

folder = "/home/cmeren/rotation_2020/COVID/"

## Read in clinicaldata. 
clinical = read.csv(paste(folder, "final_analysis/data/from_scripts/merged_clinical.csv", sep = ""), 
                                                    stringsAsFactors = F, row.names = 1)
## UniFrac analysis makes most sense on OP swabs
OP = clinical[clinical$SampleType == "Oropharyngeal swab",]

## read in UniFrac
unifrac = as.matrix(readRDS(paste(folder, "final_analysis/data/from_scripts/weighted_unifrac.rds", sep = "")))

## get first samples from this sample type
clinical_filter = clinical[clinical$SampleType == sample_type & clinical$Study_group == "COVID",]

## read in the messi data and filter columns
immune = read.csv("messi_immune_matchedOropharyngeal swab.csv", row.names = 1)
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







