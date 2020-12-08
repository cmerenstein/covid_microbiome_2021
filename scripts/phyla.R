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

## as .png 
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
                    annotation_row = covid, main = sample_type, fontsize = 20, 
                    labels_row = character(nrow(percent_filtered)))
    ## save plot
    sample_type_name = gsub(" ", "_", sample_type)
    pdf(paste("figures/heatmaps/", sample_type, "_phyla_heatmaps.pdf", sep = ""))
    print(plot)
    dev.off()
}

## -------------- BOXPLOT ----------------------------- 
top_taxa_list = c("p__Firmicutes",  "p__Bacteroidota", "p__Proteobacteria", "p__Actinobacteriota")

## Same figure as separate PNG images
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

## -------------- Broke down by ICU vs non ICU -----------------------
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
    percent_tidy$grouping = ifelse(is.na(clinical_filter$Max.WHO.score), "non COVID ICU",
                              ifelse(clinical_filter$Max.WHO.score < 7, "COVID not ICU", 
                                ifelse(clinical_filter$Max.WHO.score == 10, "Dead", "COVID ICU")))
    percent_tidy$grouping = factor(percent_tidy$grouping, levels = c("non COVID ICU", "COVID not ICU", 
                                                                     "COVID ICU", "Dead"))
     
    plot <- ggplot(percent_tidy, aes( y = percent, x = grouping, color = grouping)) +
                theme_bw() +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(height = 0, width = .2) +
                ggtitle(sample_type) + 
                facet_wrap(~phyla, ncol = 2) + 
                theme(text= element_text(size = 20)) +
                theme(axis.text.x = element_blank()) +
                scale_color_manual(values = c("green", "salmon", "orchid", "red4"))

    sample_type = gsub(" ", "_", sample_type)
    pdf(paste( "figures/phyla_boxplots/", sample_type, "_covid_grouping_boxplots.pdf", sep = ""))
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

    ## grouping
    percent_tidy$grouping = ifelse(is.na(clinical_filter$Max.WHO.score), "non COVID ICU",
                              ifelse(clinical_filter$Max.WHO.score < 7, "COVID not ICU", 
                                ifelse(clinical_filter$Max.WHO.score == 10, "Dead", "COVID ICU")))
    percent_tidy$grouping = factor(percent_tidy$grouping, levels = c("non COVID ICU", "COVID not ICU", 
                                                                     "COVID ICU", "Dead"))
     
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





