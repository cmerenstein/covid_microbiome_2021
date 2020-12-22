library(ape)
library(dplyr)
library(phyloseq)
library(vegan)

setwd("../")

clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)
clinical = clinical[clinical$Study_group != "Control",]
clinical$grouping = ifelse( is.na(clinical$Max.WHO.score), "non-COVID ICU",
                      ifelse( clinical$Max.WHO.score < 7, "COVID non ICU",
                        ifelse( clinical$Max.WHO.score == 10, "Dead", "COVID ICU")))

## load unifrac distances
unifrac = as.matrix(readRDS("data/from_scripts/unweighted_unifrac.rds"))
unifrac = unifrac[rownames(unifrac) %in% rownames(clinical), colnames(unifrac) %in% rownames(clinical)]
clinical = clinical[rownames(unifrac),] ## unifrac has fewer samples cause of >1000 reads threshold

weighted_unifrac = as.matrix(readRDS("data/from_scripts/weighted_unifrac.rds"))
weighted_unifrac = weighted_unifrac[rownames(clinical),rownames(clinical)]

## First check by sample type
adonis(as.dist(unifrac) ~ SampleType, data = clinical, strata = clinical$SubjectID)
adonis(as.dist(weighted_unifrac) ~ SampleType, data = clinical, strata = clinical$SubjectID)

## make function for pairwise comparisons
pairwise_adonis = function(unifrac, grouping){
    groupings = unique(grouping)
    tested = character() ## so we don't repeat pairs

    ## we want all non-redundant pairs of tests
    pvals_list_1 = lapply(groupings, function(grouping_1){

        tested <<- c(tested, grouping_1) ## slow but clear, prevents repeated pairs

        other_groups = groupings[!(groupings %in% tested)] 
        pvals_list_2 = lapply(other_groups, function(grouping_2){
 
            ## filter to only the 2 we're comparing
            unifrac_groups = unifrac[grouping %in% c(grouping_1, grouping_2),
                                     grouping %in% c(grouping_1, grouping_2)]
            groupings_filtered = grouping[grouping %in% c(grouping_1, grouping_2)]

            permanova = adonis(as.dist(unifrac_groups) ~ groupings_filtered, permutations = 5000)
            pval = permanova$aov.tab[1, "Pr(>F)"]
            return(data.frame(pair = paste(grouping_1, grouping_2, sep = " | "), pval = pval))
        })
        ## merge pvals into data frame before returning
        return(do.call("rbind", pvals_list_2))
    })
    return(do.call("rbind", pvals_list_1))
}


mean_pvalues = list() ## list of all comparisons for FDR correction later
## do one sample type at a time, UNWEIGHTED
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    clinical_filtered = clinical[clinical$SampleType == sample_type,]
    unifrac_filtered = unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]
    
    ## Randomly subsample one sample per patient
    pvals_list = lapply(seq(1:1000), function(i){
        subsampled = clinical_filtered %>% group_by(SubjectID) %>% sample_n(1) %>%
                                                     ungroup() %>% as.data.frame()
        unifrac_subsampled = unifrac_filtered[subsampled$SampleID, subsampled$SampleID]

        ## pairwise permanova for each pair of gropus
        permanova = pairwise_adonis(unifrac_subsampled, subsampled$grouping)
        return(permanova)
    })
    pvals = do.call("rbind", pvals_list)
    means = group_by(pvals, pair) %>% summarize(mean_p = mean(pval)) %>% as.data.frame()
    means$sample_type = sample_type
    means$unifrac_type = "unweighted"

    ## return dataframe to master list of all comparisons (to join with weighted for FDR correction)
    mean_pvalues[[paste("unweighted", sample_type, sep = "_")]] <- means
}

## WEIGHTED
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    clinical_filtered = clinical[clinical$SampleType == sample_type,]
    weighted_unifrac_filtered = weighted_unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]

    ## Randomly subsample one sample per patient
    pvals_list = lapply(seq(1:1000), function(i){
        subsampled = clinical_filtered %>% group_by(SubjectID) %>% sample_n(1) %>%
                                                     ungroup() %>% as.data.frame()
        unifrac_subsampled = weighted_unifrac_filtered[subsampled$SampleID, subsampled$SampleID]

        ## pairwise permanova for each pair of gropus
        permanova = pairwise_adonis(unifrac_subsampled, subsampled$grouping)
        return(permanova)
    })
    pvals = do.call("rbind", pvals_list)
    means = group_by(pvals, pair) %>% summarize(mean_p = mean(pval)) %>% as.data.frame()
    means$sample_type = sample_type
    means$unifrac_type = "weighted"

    ## return dataframe to master list of all comparisons (to join with unweighted for FDR correction)
    mean_pvalues[[paste("weighted", sample_type, sep = "_")]] <- means

}

permanova_results = do.call("rbind", mean_pvalues)
permanova_results$FDR = p.adjust(permanova_results$mean_p, "BH")
write.csv(permanova_results, "data/from_scripts/permanova_results_pairwise.csv", row.names = F)


