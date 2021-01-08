library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(stringr)
library(pheatmap)

setwd("../")

## get metadata
clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1)

## get counts at the GENUS levelf
counts = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", row.names = 1, stringsAsFactors = F)

# order counts the same as metadata
counts = counts[rownames(clinical),]

## FILTER TO > 1000 READS
counts = counts[rowSums(counts) > 1000,]
clinical = clinical[rownames(counts),]

## ------------- HEATMAP -----------------------------
percent = t(apply(counts, 1, function(r){ return(r / sum(r)) }))

## as pdf
for (sample_type in unique(clinical$SampleType)){
    
    ## filter to sample type, also to mean of > 0.01
    percent_filtered = percent[clinical$SampleType == sample_type, 
                                    colMeans(percent) >= 0.01]
 
    ## annotate heatmap by covid status
    covid = clinical[rownames(percent_filtered), c("Study_group", "Max.WHO.score", "Intubated.", 
                                                    "Age", "BMI.at.admission")]
    covid = as.data.frame(covid)
    rownames(covid) = rownames(percent_filtered)
    
    ## scale = "none" means color correponds to percent
    plot <- pheatmap(mat = percent_filtered, scale = "none",
                    annotation_row = covid, main = sample_type, fontsize = 12, 
                    labels_row = character(nrow(percent_filtered)))
    ## save plot
    sample_type = gsub(" ", "_", sample_type)
    pdf(paste( "figures/heatmaps/", sample_type, "_genus_heatmaps.pdf", sep = ""),
        height = 7.5, width = 7.5)
    print(plot)
    dev.off()
}

## -------------------- ETA heatmap for figure 3A ---------------
ETA = clinical[clinical$SampleType == "Endotracheal aspirate" & clinical$Study_group == "COVID",]
ETA_percent = percent[rownames(ETA),]
ETA_percent = ETA_percent[, colMeans(ETA_percent) >= 0.01]

ETA_annotation = ETA[, c("Outcome", "SubjectID"), drop = F]

pdf("figures/heatmaps/Endotracheal_aspirate_genus_heatmaps.pdf", height = 5, width = 8)
pheatmap(mat = t(ETA_percent), scale = "none", 
            annotation_col = ETA_annotation, cluster_cols = F, cluster_rows = F,
            treeheight_row = 0, treeheight_col = 0, fontsize_rows = 12,
            labels_col = ETA$SubjectID, fontsize_col = 8)
dev.off()

## Look at most common taxa to see which ones might be associated with severity
## Only look at the first sample from each patient for starters

spearman_correlation_list = list()
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab")){

    ## get first sample of this sample type COVID ONLY
    clinical_filtered = clinical[clinical$SampleType == sample_type & !(is.na(clinical$Max.WHO.score)),]
    clinical_filtered$first_visit = ave(clinical_filtered$visit_int, clinical_filtered$SubjectID, FUN = min)
    first_visit = clinical_filtered[clinical_filtered$visit_int == clinical_filtered$first_visit,]

    ## filter to sample type, remove "Unassigned"
    percent_filtered = percent[rownames(first_visit), colnames(percent) != "Unassigned"]
    ## filter by > 0.5% mean
    percent_filtered = percent_filtered[,colMeans(percent_filtered) >= 0.005] 
 
    who_score = clinical[rownames(percent_filtered), "Max.WHO.score"]
    for (genus in colnames(percent_filtered)){
        spearman = cor.test(y = percent_filtered[,genus], x = who_score, method = "spearman")
        df = data.frame(p = spearman$p.value, rho = spearman$estimate, sample_type, genus)
        spearman_correlation_list[[paste(genus, sample_type, sep = "_")]] <- df
    }
}
spearman_correlations = do.call("rbind", spearman_correlation_list)
spearman_correlations$fdr = p.adjust(spearman_correlations$p, "BH")

sig_cor = spearman_correlations[spearman_correlations$fdr < 0.05,]
saveRDS(sig_cor, "data/from_scripts/sig_genera_correlations_first_sample.rds")

## plot who score and significant correlations
pdf("figures/genus/significant_genera_first_sample.pdf", height = 6, width = 6)

par(mfrow = c(2, 2))

for (sig in rownames(sig_cor)){
    sample_type = sig_cor[sig, "sample_type"]
    genus = sig_cor[sig, "genus"]
    
    ## get first sample of this sample type
    clinical_filtered = clinical[clinical$SampleType == sample_type & !(is.na(clinical$Max.WHO.score)),]
    clinical_filtered$first_visit = ave(clinical_filtered$visit_int, clinical_filtered$SubjectID, FUN = min)
    first_visit = clinical_filtered[clinical_filtered$visit_int == clinical_filtered$first_visit,]
   
    ## get genus counts and who scores 
    genus_percent = percent[rownames(first_visit), genus]
    who_score = first_visit$Max.WHO.score
    fdr = as.character(round(sig_cor[sig, "fdr"], digits = 3))
    rho = as.character(round(sig_cor[sig, "rho"], digits = 3))
    
    ## plot with some x jitter
    plot(x = who_score + rnorm(length(who_score), 0, 0.1), y = genus_percent, ylab = "Percent genus first sample",
         xlab = "Max WHO score", main = paste(genus, "fdr:", fdr, "rho:", rho, "\n", sample_type))
    abline(lm(genus_percent ~ who_score))
    
}
dev.off()


