library(ape)
library(dplyr)
library(phyloseq)
library(vegan)

setwd("../")

clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)
clinical = clinical[clinical$Study_group != "Control",]
clinical$grouping = ifelse( is.na(clinical$Max.WHO.score), "non-COVID ICU",
                      ifelse( clinical$Max.WHO.score < 7, "COVID non ICU", "COVID ICU"))

## load unifrac distances
unifrac = as.matrix(readRDS("data/from_scripts/unweighted_unifrac.rds"))
unifrac = unifrac[rownames(unifrac) %in% rownames(clinical), colnames(unifrac) %in% rownames(clinical)]
clinical = clinical[rownames(unifrac),] ## unifrac has fewer samples cause of >1000 reads threshold

weighted_unifrac = as.matrix(readRDS("data/from_scripts/weighted_unifrac.rds"))
weighted_unifrac = weighted_unifrac[rownames(clinical),rownames(clinical)]

## First check by sample type
adonis(as.dist(unifrac) ~ SampleType, data = clinical, strata = clinical$SubjectID)
adonis(as.dist(weighted_unifrac) ~ SampleType, data = clinical, strata = clinical$SubjectID)

mean_pvalues = list()
## do one sample type at a time, UNWEIGHTED
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    clinical_filtered = clinical[clinical$SampleType == sample_type,]
    unifrac_filtered = unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]
    
    ## Randomly subsample one sample per patient
    pvals = sapply(seq(1:1000), function(i){
        subsampled = clinical_filtered %>% group_by(SubjectID) %>% sample_n(1) %>%
                                                     ungroup() %>% as.data.frame()
        unifrac_subsampled = unifrac_filtered[subsampled$SampleID, subsampled$SampleID]

        ## adonis() function computes PERMANOVA
        permanova = adonis(as.dist(unifrac_subsampled) ~ grouping, data = subsampled, permutations  = 5000)
        pval = permanova$aov.tab[1, "Pr(>F)"]
        return(pval)
    })
    pvals_stats = data.frame(unifrac_type = "unweighted", sample_type, mean_pvalue = mean(pvals), 
                              percent_lt_0.05 = (sum(pvals < .05) / length(pvals)))
    mean_pvalues[[paste("unweighted", sample_type, sep = "_")]] <- pvals_stats
}

## WEIGHTED
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    clinical_filtered = clinical[clinical$SampleType == sample_type,]
    weighted_unifrac_filtered = weighted_unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]

    ## Randomly subsample one sample per patient
    pvals = sapply(seq(1:1000), function(i){
        subsampled = clinical_filtered %>% group_by(SubjectID) %>% sample_n(1) %>%
                                                     ungroup() %>% as.data.frame()
        unifrac_subsampled = weighted_unifrac_filtered[subsampled$SampleID, subsampled$SampleID]

        ## adonis() function computes PERMANOVA
        permanova = adonis(as.dist(unifrac_subsampled) ~ grouping, data = subsampled, permutations  = 5000)
        pval = permanova$aov.tab[1, "Pr(>F)"]
        return(pval)
    })
    pvals_stats = data.frame(unifrac_type = "weighted", sample_type, mean_pvalue = mean(pvals), 
                              percent_lt_0.05 = (sum(pvals < .05) / length(pvals)))
    mean_pvalues[[paste("weighted", sample_type, sep = "_")]] <- pvals_stats

}

permanova_results = do.call("rbind", mean_pvalues)
write.csv(permanova_results, "data/from_scripts/permanova_results_covid_grouping.csv", row.names = F)


