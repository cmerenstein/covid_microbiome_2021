library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(ape)

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
sample_types = rbind( controls[,c("SampleType", "Study_group", "Max.WHO.score", "SubjectID")],
                      covid_meta[,c("SampleType", "Study_group", "Max.WHO.score", "SubjectID")])

## 
pvalues_list = list()
for(SampleType in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    samples = rownames(sample_types)[sample_types$SampleType == SampleType]
    
    for (unifrac_type in c("unweighted", "weighted", "threshold_1percent")){
    #   get unifrac and filter to sample type
        unifrac = readRDS(paste("data/from_scripts/healthy_controls/", unifrac_type, "_unifrac.rds", sep = ""))
        unifrac = as.matrix(unifrac)

        ## filter
        unifrac_samples = samples[samples %in% rownames(unifrac)]
        unifrac = unifrac[unifrac_samples, unifrac_samples]
        meta_filtered = sample_types[unifrac_samples,] 

        ## get WHO score and distance to healthy for covid patients
        covid = meta_filtered[meta_filtered$Study_group == "COVID",]
        who_score = covid$Max.WHO.score
        patient = covid$SubjectID

        ## mean distance to healthy
        healthy = rownames(meta_filtered)[grepl("Healthy", meta_filtered$Study_group)]
        distance_to_healthy = rowMeans(unifrac[rownames(covid), healthy])
        mean_distance = ave(distance_to_healthy, patient, FUN = mean)
        MD_df = unique(data.frame(mean_distance, who_score, patient))
                
        ## spearman correlation between WHO score and distance to healthy
        spearman = cor.test(x = MD_df$who_score, y = MD_df$mean_distance, method = "spearman")
        test = paste(SampleType, unifrac_type, sep = "_") 
        pvalues_list[[test]] = data.frame(test, rho = spearman$estimate, p = spearman$p.value)

        ## plot 
        plot = ggplot(MD_df, aes(x = who_score, y = mean_distance)) + 
                theme_classic() + 
                geom_point(pch = 21, size = 3) + 
                stat_smooth(method = "lm", color = "black") +
                ggtitle(paste(SampleType, unifrac_type)) 

        st = gsub(" ", "_", SampleType)
        pdf(paste("figures/454_healthy_controls/COVID_distances/", st, "_", 
                    unifrac_type, "_unifrac.pdf", sep = ""))
        print(plot)
        dev.off()
    }
}
pvalues = do.call("rbind", pvalues_list)
pvalues$FDR = p.adjust(pvalues$p, method = "BH")
print(pvalues)
