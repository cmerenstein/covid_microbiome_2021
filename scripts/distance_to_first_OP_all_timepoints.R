library(ggplot2)
library(dplyr)
library(tidyr)

setwd("../")

## Read in merged clinical data
meta = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## weighted unifrac distance
unifrac = as.matrix(readRDS("data/from_scripts/weighted_unifrac.rds"))

## align samples
unifrac = unifrac[intersect(rownames(meta), rownames(unifrac)), intersect(rownames(meta), rownames(unifrac))] 
meta = meta[rownames(unifrac),] ## subsetted samples to >1000 reads for unifrac

## Subset to each sample type and do each one at a time
sample_type = "Oropharyngeal swab"

## filter to this sample type and just COVID
meta_filter = meta[meta$SampleType == sample_type & meta$Study_group != "Control",]
unifrac_filter = unifrac[rownames(meta_filter),rownames(meta_filter)]

## get the first sample date
meta_filter$Date = as.Date(meta_filter$Collection_date, format = "%m/%d/%Y") 
meta_filter$first_collection = ave(meta_filter$Date, meta_filter$SubjectID, FUN = min)
meta_filter$days_since_first = as.numeric(meta_filter$Date - meta_filter$first_collection)


## Get only in samples with multiple timepoints
n_samples = table(meta_filter$SubjectID)
multiple_samples = n_samples[n_samples > 1]
meta_filter = meta_filter[meta_filter$SubjectID %in% names(multiple_samples),]
unifrac_filter = unifrac_filter[rownames(meta_filter), rownames(meta_filter)]

## Get the ID of the first sample so we can match it by subject ID
first_sample = meta_filter[meta_filter$days_since_first == 0,]
first_sample_id = rownames(first_sample)
names(first_sample_id) = first_sample$SubjectID

## get distances
distances = lapply(rownames(meta_filter), function(s){
    first_sample_name = first_sample_id[meta_filter[s, "SubjectID"]]
    days_since_first = meta_filter[s, "days_since_first"]
    
    ## unifrac to first
    dist = unifrac_filter[first_sample_name, s]
    return(data.frame(days_since_first = days_since_first, weighted_unifrac = dist,
                        covid = meta_filter[s, "Study_group"],
                        patient = meta_filter[s, "SubjectID"]))
    
}) 
distances = do.call("rbind", distances)

## some control samples have 2 samples, but both at the same timepoint. Drop these
distances$max = ave(distances$days_since_first, distances$patient, FUN = max)
distances = distances[distances$max > 0,]

## global median
global_mean = mean(unifrac_filter)

## plot the distances over time
pdf("figures/timepoint/OP_time_since_first_COVID_only.pdf")
ggplot(distances, aes(x = days_since_first, y = weighted_unifrac, color = covid)) +
    geom_point() + 
    stat_smooth(span = .5) + 
    theme_bw() + 
    theme(text = element_text(size = 20)) + 
    ggtitle("Oropharyngeal swabs") +
    ylab("Weighted UniFrac Distance to First Sample") +
    xlab("Days Since First Sample") +
    geom_line(y = global_mean, color = "black")
dev.off()




