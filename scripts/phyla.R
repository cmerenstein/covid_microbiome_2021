library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(stringr)
library(pheatmap)

## R 4.0.2

setwd("../") ## set wd for access to data

## load clinical data
clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)

## get Phyla counts
counts = read.csv( "data/raw/counts_taxa_level_1.csv", row.names = 1, stringsAsFactors = F)

## filter to match clinical
counts = counts[rownames(clinical),]


## ------------- HEATMAP -----------------------------
percent = t(apply(counts, 1, function(r){ return(r / sum(r)) }))

## plot phylum heatmap for each sample type 
for (sample_type in unique(clinical$SampleType)){
    
    ## filter to sample type, also to mean of > 0.001
    percent_filtered = percent[clinical$SampleType == sample_type, 
                                    colMeans(percent) >= 0.001]
 
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
    sample_type_name = gsub(" ", "_", sample_type)
    pdf(paste("figures/heatmaps/", sample_type, "_phyla_heatmaps.pdf", sep = ""))
    print(plot)
    dev.off()
}

## -------------- BOXPLOT ----------------------------- 
top_taxa_list = c("p__Firmicutes",  "p__Bacteroidota", "p__Proteobacteria", "p__Actinobacteriota")

## Same figure as separate  images
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    ## filter to only common taxa and the right sample type
    percent_filter = percent[clinical$SampleType == sample_type,
                            colnames(percent) %in% top_taxa_list]

    ## get into tidy format for plotting
    percent_filter = as.data.frame(percent_filter)
    percent_filter$sample = rownames(percent_filter)
    percent_tidy = gather(percent_filter, "phyla", "percent", -sample)

    ## add who score
    percent_tidy$who_score = ifelse(is.na(clinical[percent_tidy$sample, "Max.WHO.score"]), 0,
                                    clinical[percent_tidy$sample, "Max.WHO.score"]) 
    
    plot <- ggplot(percent_tidy, aes( y = percent, x = who_score)) +
                theme_bw() + 
                geom_jitter(height = 0, width = .2) + 
                ggtitle(sample_type) + 
                facet_wrap(~phyla, ncol = 2) +
                theme(text = element_text(size = 20))
    ## save plot
    sample_type = gsub(" ", "_", sample_type)
    pdf(paste("figures/phyla_boxplots/" , sample_type, "_phyla_boxplots.pdf"))
    print(plot)
    dev.off()
}

## -------------- Broke down by critical/severe/dead -----------------------
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    ## filter to only common taxa and the right sample type
    percent_filter = percent[clinical$SampleType == sample_type,
                            colnames(percent) %in% top_taxa_list]

    ## get into tidy format for plotting
    percent_filter = as.data.frame(percent_filter)
    percent_filter$sample = rownames(percent_filter)

    ## filter clinical
    clinical_filter = clinical[percent_filter$sample,]

    ## in tidy format
    percent_tidy = gather(percent_filter, "phyla", "percent", -sample)
    percent_tidy$phyla = gsub("p__", "", percent_tidy$phyla) 

    ## grouping
    percent_tidy$grouping = ifelse(is.na(clinical_filter$Max.WHO.score), "non COVID",
                              ifelse(clinical_filter$Max.WHO.score < 7, "COVID mod/severe", 
                                ifelse(clinical_filter$Max.WHO.score == 10, "COVID Dead", "COVID critical")))
    percent_tidy$grouping = factor(percent_tidy$grouping, levels = c("non COVID", "COVID mod/severe", 
                                                                     "COVID critical", "COVID Dead"))
     
    plot <- ggplot(percent_tidy, aes( y = percent, x = grouping, color = grouping)) +
                theme_bw() +
                geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = .2) +
                ggtitle(sample_type) + facet_wrap(~phyla, ncol = 2) + 
                theme(text= element_text(size = 20)) + theme(axis.text.x = element_blank()) +
                scale_color_manual(values = c("forestgreen", "salmon", "orchid", "red4")) +
                ylab("Relative Abundance (Percent)") + xlab("COVID Status and Severity") +
                theme(legend.text = element_text(size = 14)) + theme(legend.title = element_blank())  +
                ylim(0, 1.1)
     #           geom_signif(comparisons = list(c("non COVID", "COVID mod/severe")), 
     #                           y = .8, tip_length = 0, color = "black", textsize = 3) +
     #           geom_signif(comparisons = list(c("non COVID", "COVID critical")), 
     #                           y = .9, tip_length = 0, color = "black", textsize = 3) +
     #           geom_signif(comparisons = list(c("non COVID", "COVID Dead")),
     #                           y = 1.2, tip_length = 0, color = "black", textsize = 3) +
     #           geom_signif(comparisons = list(c("COVID mod/severe", "COVID critical")),
     #                           y = 1, tip_length = 0, color = "black", textsize = 3) +
     #           geom_signif(comparisons = list(c("COVID mod/severe", "COVID Dead")),
     #                           y = 1.1, tip_length = 0, color = "black", textsize = 3) +
     #           geom_signif(comparisons = list(c("COVID critical", "COVID Dead")),
     #                           y = .8, tip_length = 0, color = "black", textsize = 3) 
    sample_type = gsub(" ", "_", sample_type)
    pdf(paste( "figures/phyla_boxplots/", sample_type, "_covid_grouping_boxplots_tall.pdf", sep = ""))
    print(plot)
    dev.off()
}

## --------------- Check significance with random sampling --------------
## record all p values for later FDR
pvalues_list = list()
for (sample_type in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    ## filter to only common taxa and the right sample type
    percent_filter = percent[clinical$SampleType == sample_type,
                            colnames(percent) %in% top_taxa_list]

    ## get into tidy format for plotting
    percent_filter = as.data.frame(percent_filter)
    percent_filter$sample = rownames(percent_filter)

    ## filter clinical
    clinical_filter = clinical[percent_filter$sample,]

    ## in tidy format
    percent_tidy = gather(percent_filter, "phyla", "percent", -sample)
    percent_tidy$phyla = gsub("p__", "", percent_tidy$phyla) 
    percent_tidy$SubjectID = clinical_filter[percent_filter$sample, "SubjectID"]

    ## gropuing
    percent_tidy$grouping = ifelse(is.na(clinical_filter$Max.WHO.score), "non COVID",
                              ifelse(clinical_filter$Max.WHO.score < 7, "COVID mod/severe", 
                                ifelse(clinical_filter$Max.WHO.score == 10, "COVID Dead", "COVID critical")))
    percent_tidy$grouping = factor(percent_tidy$grouping, levels = c("non COVID", "COVID mod/severe", 
                                                                     "COVID critical", "COVID Dead"))
    ## random sampleing for significance
    for (i in seq(0:1000)){
        subsample = percent_tidy %>% group_by(phyla, SubjectID) %>% sample_n(1) %>%
                                     ungroup() %>% as.data.frame()
        ## test each sample
        for (phylum in unique(subsample$phyla)) {
            phylum_percent = subsample[subsample$phyla == phylum,]
            p = kruskal.test(percent ~ grouping, data = phylum_percent)$p.value
            pvalues_list[[paste(i, phylum, sample_type)]] <- data.frame(phylum, p, sample_type)
        }
    }
}
pvalues = do.call("rbind", pvalues_list)

## get mean within each test, then do FDR correction
pvalues_mean = group_by(pvalues, phylum, sample_type) %>% summarize(mean_p = mean(p)) %>%
                ungroup() %>% as.data.frame()
pvalues_mean$FDR = p.adjust(pvalues_mean$mean_p)

## ---------- First sample only, WHO score ---------------
spearman_correlation_list = list()
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab")){

    ## get first sample of this sample type COVID ONLY
    clinical_filtered = clinical[clinical$SampleType == sample_type & !(is.na(clinical$Max.WHO.score)),]
    clinical_filtered$first_visit = ave(clinical_filtered$visit_int, clinical_filtered$SubjectID, FUN = min)
    first_visit = clinical_filtered[clinical_filtered$visit_int == clinical_filtered$first_visit,]

    ## filter to sample type, remove "Unassigned"
    percent_filtered = percent[rownames(first_visit), colnames(percent) != "Unassigned"]
    
    ## filter to 2% mean
    percent_filtered = percent_filtered[,colMeans(percent_filtered) >= 0.02]

    who_score = clinical[rownames(percent_filtered), "Max.WHO.score"]
    for (phylum in colnames(percent_filtered)){
        spearman = cor.test(y = percent_filtered[,phylum], x = who_score, method = "spearman")
        df = data.frame(p = spearman$p.value, rho = spearman$estimate, sample_type, phylum)
        spearman_correlation_list[[paste(phylum, sample_type, sep = "_")]] <- df
    }
}
spearman_correlations = do.call("rbind", spearman_correlation_list)
spearman_correlations$fdr = p.adjust(spearman_correlations$p, "BH")



