library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)

setwd("../")

## Read in clinicaldata. 
clinical = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## read unifrac data
UU = as.matrix(readRDS("data/from_scripts/unweighted_unifrac.rds"))
UU = UU[rownames(UU) %in% rownames(),]

## subsetted samples to >1000 reads for unifrac
clinical = clinical[rownames(UU),]

## read in the clinical labs data
labs = read.csv("data/raw/COVID_Metadata_22102020SW.csv")
labs$SubjectID = gsub("-", "", labs$CORE_SubjectID)
labs$Date = as.Date(labs$RESULT_DATE, format = "%m/%d/%Y")
labs$ORD_VALUE = as.numeric(labs$ORD_VALUE)

# -------------- UNWEIGHTED UNIFRAC PCOA AXES --------------------------------------------

## Correlation with UniFrac PCoA axes. Only OP swabs really make sense for unweighted
pvalues_unweighted = list()
sample_type = "Oropharyngeal swab" 
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

    ## get the pcoa axes
    set.seed(66) ## always set a lucky seed
    UU_filtered = UU[sample_mapping$SampleID, sample_mapping$SampleID]
    ordinates = pcoa(UU_filtered)
    vectors = ordinates$vectors
    pcoa_df = data.frame(SampleID = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"])

    ## add pcoa axes to sample matching data frame
    sample_mapping_pcoa = merge(sample_mapping, pcoa_df, by = "SampleID")

    ## add labs
    labs_matched = merge(sample_mapping_pcoa, lab, by.x = c("SubjectID", "first_lab"), 
                                               by.y = c("SubjectID", "Date"))

    ## Spearman Correlation
    spearman_1 = cor.test(labs_matched$Axis.1, labs_matched$ORD_VALUE, method = "spearman")
    spearman_2 = cor.test(labs_matched$Axis.2, labs_matched$ORD_VALUE, method = "spearman")
    pvalues_unweighted[[paste(test, sample_type, "axis1")]] <- data.frame(sample_type, test, 
                                            unifrac = "unweighted", pcoa = "Axis.1", 
                                            rho = spearman_1$estimate, pval = spearman_1$p.value)
                                                    
    pvalues_unweighted[[paste(test, sample_type)]] <- data.frame(sample_type, test, 
                                            unifrac = "unweighted", pcoa = "Axis.2",
                                            rho = spearman_2$estimate, pval = spearman_2$p.value)
}
cor_unweighted = do.call("rbind", pvalues_unweighted)

## ------------------- WEIGHTED UNIFRAC, USE BOTH NP AND OP ---------------
## read unifrac data
WU = as.matrix(readRDS("data/from_scripts/weighted_unifrac.rds"))
WU = WU[rownames(WU) %in% rownames(),]
    
## Correlation with UniFrac PCoA axes. 
pvalues_weighted = list()
for (sample_type in c( "Oropharyngeal swab" , "Nasopharyngeal swab")){
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

        ## get the pcoa axes
        set.seed(66) ## always set a lucky seed
        WU_filtered = WU[sample_mapping$SampleID, sample_mapping$SampleID]
        ordinates = pcoa(WU_filtered)
        vectors = ordinates$vectors
        pcoa_df = data.frame(SampleID = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                            Axis.2 = vectors[,"Axis.2"])

        ## add pcoa axes to sample matching data frame
        sample_mapping_pcoa = merge(sample_mapping, pcoa_df, by = "SampleID")

        ## add labs
        labs_matched = merge(sample_mapping_pcoa, lab, by.x = c("SubjectID", "first_lab"), 
                                                   by.y = c("SubjectID", "Date"))

        ## Spearman Correlation
        spearman_1 = cor.test(labs_matched$Axis.1, labs_matched$ORD_VALUE, method = "spearman")
        spearman_2 = cor.test(labs_matched$Axis.2, labs_matched$ORD_VALUE, method = "spearman")
        pvalues_weighted[[paste(test, sample_type, "axis1")]] <- data.frame(sample_type, test, 
                                                unifrac = "weighted", pcoa = "Axis.1", 
                                                rho = spearman_1$estimate, pval = spearman_1$p.value)
                                                        
        pvalues_weighted[[paste(test, sample_type)]] <- data.frame(sample_type, test, 
                                                unifrac = "weighted", pcoa = "Axis.2",
                                                rho = spearman_2$estimate, pval = spearman_2$p.value)
    }
}
cor_weighted = do.call("rbind", pvalues_weighted)

## all correlations and FDR correct
lab_correlations = rbind(cor_weighted, cor_unweighted)
lab_correlations$FDR = p.adjust(lab_correlations$pval)

