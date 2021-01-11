library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

setwd("../")

## get metadata
clinical = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## get counts at the GENUS level
counts = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", row.names = 1, stringsAsFactors = F)

# order counts the same as clinicaldata
counts = counts[rownames(clinical),]

## FILTER TO > 1000 READS
counts = counts[rowSums(counts) > 1000,]
clinical = clinical[rownames(counts),] 

## get alpha diversity for each sample
simpson = diversity(counts, "simpson")

## differences between sample types
mean(simpson)
aggregate(simpson, by = list(clinical$SampleType), FUN = mean)
kruskal.test(simpson ~ clinical$SampleType)

## By sample type
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab", "Endotracheal aspirate")){

    ## first sample of this sample type
    clinical_filtered = clinical[clinical$SampleType == sample_type & !(is.na(clinical$Max.WHO.score)),]
    clinical_filtered$first_visit = ave(clinical_filtered$visit_int, clinical_filtered$SubjectID, FUN = min)
    first_visit = clinical_filtered[clinical_filtered$visit_int == clinical_filtered$first_visit,]

    ## add diversity
    first_visit$simpson = simpson[rownames(first_visit)]
    
    ## spearman correlation
    print(sample_type)
    print(cor.test(x = first_visit$Max.WHO.score, y = first_visit$simpson, method = "spearman"))

    ## plot, with some x jitter for clarity
    sample_type = gsub(" ", "_", sample_type)
    pdf(paste("figures/diversity/first_sample_simpson_", sample_type, ".pdf", sep = ""))
    p = ggplot(first_visit, aes(x = Max.WHO.score, y = simpson)) + 
        geom_point(size = 2) + 
        theme_classic() + 
        stat_smooth(color = "black", method = "lm") + 
        ylab("Simpson Diversity") + 
        xlab("Maximum WHO Score") +
        theme(text = element_text(size = 20))
    print(p)
    dev.off()
}


## Also test p value when using random subsampling, only for OP swab since that was the only significant one

## first sample of this sample type
pvals = sapply(seq(1, 1000), function(i){

    clinical_filtered = clinical[clinical$SampleType == "Oropharyngeal swab" & !(is.na(clinical$Max.WHO.score)),]
    subsampled = clinical_filtered %>% group_by(SubjectID) %>% sample_n(1) %>%
                                       ungroup() %>% as.data.frame()

    subsampled$simpson = simpson[subsampled$SampleID]
    spearman = cor.test(x = subsampled$simpson, y = subsampled$Max.WHO.score, method = "spearman")
    return(spearman$p.value)
})
mean(pvals)
    

## lastly, differences in COVID vs non COVID, for each sample type
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab", "Endotracheal aspirate")){
    
    ## filter to sample type 
    clinical_filtered = clinical[clinical$SampleType == sample_type,]
    simpson_filtered = simpson[clinical_filtered$SampleID]
    
    ## differences between covid / non-covid
    print(sample_type)
    print(aggregate(simpson_filtered, by = list(clinical_filtered$Study_group), FUN = mean))
    print(kruskal.test(simpson_filtered ~ clinical_filtered$Study_group))

}






