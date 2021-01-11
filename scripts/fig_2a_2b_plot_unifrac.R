library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(ape)

setwd("../") ## set working directory to project root

## get COVID metadata
clinical = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)
clinical = clinical[clinical$Study_group != "Control",]

## WHO score of NA means non-COVID, 4-6 is hosptialized not intubated, 7-9 = intubated, 10 = dead
clinical$grouping = ifelse(is.na(clinical$Max.WHO.score), "non COVID", 
                      ifelse(clinical$Max.WHO.score < 7, "COVID mod/severe",
                            ifelse(clinical$Max.WHO.score == 10, "COVID Dead", "COVID critical")))

## load unifrac distances
unifrac = as.matrix(readRDS("data/from_scripts/unweighted_unifrac.rds"))
unifrac = unifrac[rownames(unifrac) %in% rownames(clinical), colnames(unifrac) %in% rownames(clinical)]
clinical = clinical[rownames(unifrac),] ## unifrac has fewer samples cause of >1000 reads threshold

weighted_unifrac = as.matrix(readRDS("data/from_scripts/weighted_unifrac.rds"))
weighted_unifrac = weighted_unifrac[rownames(clinical),rownames(clinical)]

## ------ UNWEIGHTED -----
for(SampleType in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    clinical_filtered = clinical[clinical$SampleType == SampleType,]
    unifrac_filtered = unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]

    ## PCoA
    set.seed(66) ## always set a lucky seed
    ordinates = pcoa(unifrac_filtered)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = factor(clinical_filtered[rownames(vectors),"grouping"], 
                                 levels = c("non COVID", "COVID critical", "COVID Dead", "COVID mod/severe")))
    ## save PCoA axes for correlation to taxa
    SampleType = gsub(" ", "_", SampleType)
    saveRDS(pcoa_df, paste( "data/from_scripts/pcoa/", SampleType, "_unweighted.rds", sep = ""))

    ## define centroids for plotting
    centroids = aggregate(cbind(Axis.1, Axis.2) ~ group, data = pcoa_df, mean)

    ## plot
    pdf(paste("figures/unifrac/", SampleType, "_unweighted_unifrac_pcoa.pdf", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = group)) + 
        theme_classic() +
        geom_point(size = 4, alpha = .5) +
        ggtitle(SampleType) + 
        scale_color_manual(values = c("forestgreen", "orchid", "red4", "salmon")) + 
        geom_point(data = centroids, aes(x = Axis.1, y = Axis.2, color = group), stroke = 3, shape = 4, size = 4) +
        scale_color_manual(values = c("forestgreen", "orchid", "red4", "salmon")) +
        theme(legend.position = "bottom")

    print(plot)
    dev.off()
}


## ------ WEIGHTED -----
for(SampleType in c("Endotracheal aspirate", "Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    clinical_filtered = clinical[clinical$SampleType == SampleType,]
    unifrac_filtered = weighted_unifrac[rownames(clinical_filtered), rownames(clinical_filtered)]

    ## PCoA
    set.seed(66)
    ordinates = pcoa(unifrac_filtered)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = factor(clinical_filtered[rownames(vectors),"grouping"], 
                                 levels = c("non COVID", "COVID mod/severe", "COVID critical", "COVID Dead")))
    ## save PCoA axes
    SampleType = gsub(" ", "_", SampleType)
    saveRDS(pcoa_df, paste("data/from_scripts/pcoa", SampleType, "_weighted.rds", sep = "")) 

    ## define centroids for plotting
    centroids = aggregate(cbind(Axis.1, Axis.2) ~ group, data = pcoa_df, mean)

    ## plot
    pdf(paste("figures/unifrac/", SampleType, "_weighted_unifrac_pcoa.pdf", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = group)) + 
        theme_classic() +
        geom_point(size = 4, alpha = .5) +
        ggtitle(SampleType) + 
        scale_color_manual(values = c("forestgreen", "orchid", "red4", "salmon")) + 
        geom_point(data = centroids, aes(x = Axis.1, y = Axis.2, color = group), stroke = 3, shape = 4, size = 4)
        scale_color_manual(values = c("forestgreen", "orchid", "red4", "salmon")) 

    print(plot)
    dev.off()
}

