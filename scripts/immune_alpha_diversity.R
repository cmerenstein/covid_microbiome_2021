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

## read in the messi data
immune = read.csv("data/raw/messi_immune.csv", stringsAsFactors = F)
immune$SubjectID = gsub("-", "", immune$CORE_ID)
immune$Date = as.Date(immune$date_received , format = "%m/%d/%Y")


## Correlation with diversity
pvalues = list()
for (sample_type in c("Oropharyngeal swab", "Nasopharyngeal swab")){ 
    ## get first samples from this sample type
    clinical_filter = clinical[clinical$SampleType == sample_type & clinical$Study_group == "COVID",]

    ## merge clinical and immune
    clinical_dates = clinical[,c("SubjectID", "Date")]
    dates = merge(clinical_dates, immune[,c("SubjectID", "Date")],
                                by = "SubjectID", suffixes = c(".sample", ".immune"))
    dates$dif = as.Date(dates$Date.sample) - as.Date(dates$Date.immune)
    dates$min_dif = ave(abs(as.integer(dates$dif)), dates$SubjectID, FUN = min)

    ## get just the closest dates
    close_dates = dates[abs(as.integer(dates$dif)) == dates$min_dif,]

    ## filter to only within 7 days (in practice all that remain are w/in 4 days)
    close_dates = close_dates[close_dates$min_dif <= 7,]

    ## sometimes there are 2 closest samples, take the earlier one
    close_dates$first_sample = ave(close_dates$Date.sample, close_dates$SubjectID, FUN = min)
    close_dates$first_immune = ave(close_dates$Date.immune, close_dates$SubjectID, FUN = min)
    
    ## add subject ID and get simpson diversity
    sample_mapping = unique( merge(clinical_filter[, c("SubjectID", "SampleID", "Date")], close_dates, 
                            by.x = c("SubjectID", "Date"), by.y = c("SubjectID", "Date.sample")))
    sample_mapping$simpson = diversity(counts[sample_mapping$SampleID, ], "simpson")
    
    ## add immune
    immune_matched = merge(sample_mapping, immune, by.x = c("SubjectID", "first_immune"), 
                                                   by.y = c("SubjectID", "Date"))
    ## save matched data frame
    write.csv(immune_matched, paste("data/from_scripts/immune/messi_immune_matched", sample_type, ".csv", sep = ""))

    ## plot immune features vs diversity
    pdf(paste("figures/immune/immune_features_simpson_diversity_", sample_type, ".pdf", sep = ""))
    for (feature in colnames(immune_matched)[12:25]){

        ## get correlation and plot
        df = data.frame(simpson = immune_matched$simpson,
                immune_feature = immune_matched[, feature])
        correlation = cor.test(df$simpson, df$immune_feature, method = "pearson")
        pvalues[[paste(sample_type, feature)]] <- data.frame(p = correlation$p.value,
                                                   test = paste(sample_type, feature))

        plot = ggplot(df, aes(x = simpson, y = immune_feature)) +
            theme_classic() + 
            geom_point() + 
            stat_smooth(method = "lm") +
            ylab(feature) + 
            ggtitle(paste(sample_type, feature, correlation$p.value))
        
        print(plot)    
    }
    dev.off()
} 
pvalues_df = do.call("rbind", pvalues)
pvalues_df$FDR = p.adjust(pvalues_df$p, "BH")

## ----------- LOAD JUST OP ---------------------------------------------------
OP_immune = read.csv("data/from_scripts/messi_immune_matchedOropharyngeal swab.csv", row.names = 1)

pdf(paste("figures/immune/OP_simpson_perc_CD4.pdf"))
ggplot(OP_immune, aes(x = simpson, y = percOf_Tcell_CD4)) +
    theme_classic() +
    geom_point(size = 4) +
    stat_smooth(color = "black", method = "lm") +
    ylab("Percent of T Cells CD4") + 
    theme(text = element_text(size = 20))

dev.off()






