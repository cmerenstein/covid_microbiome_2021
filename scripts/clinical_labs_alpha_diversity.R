library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

setwd("../")

## Read in clinicaldata. 
clinical = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## read counts post-contamination filter
counts = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", 
                                     row.names = 1, stringsAsFactors = F)
counts = counts[rownames(counts) %in% rownames(clinical),]

## subsetted samples to >1000 reads for unifrac
counts = counts[rowSums(counts) > 1000,]
clinical = clinical[rownames(counts),]

## read in the clinical labs data
labs = read.csv("data/raw/COVID_Metadata_22102020SW.csv")
labs$SubjectID = gsub("-", "", labs$CORE_SubjectID)
labs$Date = as.Date(labs$RESULT_DATE, format = "%m/%d/%Y")
labs$ORD_VALUE = as.numeric(labs$ORD_VALUE)

## Correlation with diversity
pvalues = list()
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab")){ 
    ## get first samples from this sample type
    clinical_filter = clinical[clinical$SampleType == sample_type & clinical$Study_group == "COVID",]

    ## split out each lab type
    for (test in unique(labs$RESULT_NAME)){
        lab = labs[labs$RESULT_NAME == test, c("SubjectID", "Date", "RESULT_NAME", "ORD_VALUE")]

        ## merge clinical and labs
        clinical_dates = clinical[,c("SubjectID", "Date")]
        dates = merge(clinical_dates, lab[,c("SubjectID", "Date")],
                                    by = "SubjectID", suffixes = c(".sample", ".lab"))
        dates$dif = as.Date(dates$Date.sample) - as.Date(dates$Date.lab)
        dates$min_dif = ave(abs(as.integer(dates$dif)), dates$SubjectID, FUN = min)

        ## get just the closest dates
        close_dates = dates[abs(as.integer(dates$dif)) == dates$min_dif,]

        ## filter to only within 7 days (in practice all that remain are w/in 4 days)
        close_dates = close_dates[close_dates$min_dif <= 7,]

        ## sometimes there are 2 closest samples, take the earlier one
        close_dates$first_sample = ave(close_dates$Date.sample, close_dates$SubjectID, FUN = min)
        close_dates$first_lab = ave(close_dates$Date.lab, close_dates$SubjectID, FUN = min)

        close_dates = unique(close_dates)

        ## add subject ID and get simpson diversity
        sample_mapping = unique( merge(clinical_filter[, c("SubjectID", "SampleID", "Date")], close_dates, 
                                by.x = c("SubjectID", "Date"), by.y = c("SubjectID", "Date.sample")))
        sample_mapping$simpson = diversity(counts[sample_mapping$SampleID, ], "simpson")
        
        ## add labs
        labs_matched = merge(sample_mapping, lab, by.x = c("SubjectID", "first_lab"), 
                                                   by.y = c("SubjectID", "Date"))

        ## Spearman Correlation
        print(paste(sample_type, test))
        print(cor.test(labs_matched$simpson, labs_matched$ORD_VALUE, method = "spearman"))
    }
}

    
