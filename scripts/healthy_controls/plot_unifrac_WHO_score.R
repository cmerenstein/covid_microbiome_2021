library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(ape)

## THIS ANALYSIS DIDN"T GET INCLUDED IN THE PAPER AND THIS SCRIPT MIGHT NOT WORK RIGHT NOW
## BECAUSE IT WASN'T UPDATED AS THINGS CHANGED


setwd("../../") ## root of project

## get COVID metadata
covid_meta = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)
covid_meta = covid_meta[covid_meta$Study_group == "COVID",]

## get 454 control metadata
controls = read.csv("data/raw/454.healthy.csv", stringsAsFactors = F, row.names =1)
controls = controls[controls$respiratory_disease == "Healthy" & !(is.na(controls$respiratory_disease)),]
controls$Study_group = paste("Healthy", controls$smokingstatus, sep = "_")
controls$Max.WHO.score = NA

## combine metadata for sample types, covid/noncovid
sample_types = rbind( controls[,c("SampleType", "Study_group", "Max.WHO.score")],
                      covid_meta[,c("SampleType", "Study_group", "Max.WHO.score")])

## ------------ WEIGHTED ------------------------------------------
for(SampleType in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    samples = rownames(sample_types)[sample_types$SampleType == SampleType]
    
#   get unifrac and filter to sample type
    unifrac = readRDS("data/from_scripts/healthy_controls/weighted_unifrac.rds")
    unifrac = as.matrix(unifrac)

    ## filter
    unifrac_samples = samples[samples %in% rownames(unifrac)]
    unifrac = unifrac[unifrac_samples, unifrac_samples]

    ordinates = pcoa(unifrac)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = sample_types[rownames(vectors),"Study_group"], 
                        Max.WHO.score = sample_types[rownames(vectors), "Max.WHO.score"])

    ## plot
    SampleType = gsub(" ", "_", SampleType)
    pdf(paste("figures/454_healthy_controls/", SampleType, "COVID_weighted_unifrac_pcoa.pdf", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = Max.WHO.score)) + 
        theme_bw() +
        geom_point(shape = 21, size = 3) +
        scale_fill_gradient(low = "#f7c8c8", high = "#690000", na.value = "dodgerblue") +
        ggtitle("weighted unifrac")
    #    scale_shape_manual(values =  c(0, 1, 2, 8, 9, 11)) +
#        scale_color_manual(values = c("skyblue", "navy", "red", "forestgreen", "lightgreen"))  
    print(plot)
    dev.off()
}


## ------------ UNWEIGHTED ------------------------------------------
for(SampleType in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    samples = rownames(sample_types)[sample_types$SampleType == SampleType]
    
#   get unifrac and filter to sample type
    unifrac = readRDS("data/from_scripts/healthy_controls/unweighted_unifrac.rds")
    unifrac = as.matrix(unifrac)

    ## filter
    unifrac_samples = samples[samples %in% rownames(unifrac)]
    unifrac = unifrac[unifrac_samples, unifrac_samples]

    ordinates = pcoa(unifrac)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = sample_types[rownames(vectors),"Study_group"], 
                        Max.WHO.score = sample_types[rownames(vectors), "Max.WHO.score"])

    ## plot
    SampleType = gsub(" ", "_", SampleType)
    pdf(paste("figures/454_healthy_controls/", SampleType, "COVID_unweighted_unifrac_pcoa.pdf", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = Max.WHO.score)) + 
        theme_bw() +
        geom_point(shape = 21, size = 3) +
        scale_fill_gradient(low = "#f7c8c8", high = "#690000", na.value = "dodgerblue") +
        ggtitle("unweighted no threshold")
    #    scale_shape_manual(values =  c(0, 1, 2, 8, 9, 11)) +
#        scale_color_manual(values = c("skyblue", "navy", "red", "forestgreen", "lightgreen"))  
    print(plot)
    dev.off()
}


## ------------ UNWEIGHTED 1% THRESHOLD ------------------------------------------
for(SampleType in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    samples = rownames(sample_types)[sample_types$SampleType == SampleType]
    
#   get unifrac and filter to sample type
    unifrac = readRDS("data/from_scripts/healthy_controls/threshold_1percent_unifrac.rds")
    unifrac = as.matrix(unifrac)

    ## filter
    unifrac_samples = samples[samples %in% rownames(unifrac)]
    unifrac = unifrac[unifrac_samples, unifrac_samples]

    ordinates = pcoa(unifrac)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = sample_types[rownames(vectors),"Study_group"], 
                        Max.WHO.score = sample_types[rownames(vectors), "Max.WHO.score"])

    ## plot
    SampleType = gsub(" ", "_", SampleType)
    png(paste("figures/454_healthy_controls/", SampleType, "COVID_1percent_threshold_unifrac_pcoa.png", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = Max.WHO.score)) + 
        theme_bw() +
        geom_point(shape = 21, size = 3) +
        scale_fill_gradient(low = "#f7c8c8", high = "#690000", na.value = "dodgerblue") +
        ggtitle("unweighted 1% threshold")
    #    scale_shape_manual(values =  c(0, 1, 2, 8, 9, 11)) +
#        scale_color_manual(values = c("skyblue", "navy", "red", "forestgreen", "lightgreen"))  
    print(plot)
    dev.off()
}
