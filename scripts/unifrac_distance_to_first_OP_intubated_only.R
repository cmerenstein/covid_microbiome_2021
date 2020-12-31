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

## get just this sample type and the median distance between all samples
meta_filter = meta[meta$SampleType == sample_type & meta$Study_group != "Control",]

## filter COVID patients to only intubated patients, keep non-COVID patients
meta_filter = meta_filter[meta_filter$Intubated. == "yes" | is.na(meta_filter$Max.WHO.score),]

unifrac_filter = unifrac[rownames(meta_filter),rownames(meta_filter)]

## get the first sample date
meta_filter$Date = as.Date(meta_filter$Collection_date, format = "%m/%d/%Y") 
meta_filter$first_collection = ave(meta_filter$Date, meta_filter$SubjectID, FUN = min)
meta_filter$days_since_first = as.numeric(meta_filter$Date - meta_filter$first_collection)

## SUBSET TO TIMEPOINTS WE HAVE NON COVID FOR
meta_filter = meta_filter[meta_filter$days_since_first <=
            max( meta_filter[meta_filter$Study_group == "non COVID", "days_since_first"]), ]

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

## get the slope for each patient
slopes_list = list()
for (patient in unique(distances$patient)){
    days = distances[distances$patient == patient, "days_since_first"]
    wu = distances[distances$patient == patient, "weighted_unifrac"]
    covid = unique(distances[distances$patient == patient, "covid"])
    
    ## get the slope
    slope = coef(lm(wu ~ days))[2]

    ## some patients have multiple samples, but only at 1 time point. kick them out 
    if (!(is.na(slope))){
        slopes_list[[patient]] = data.frame(covid, patient, slope) 
    }
}
slopes = do.call("rbind", slopes_list)

## Compare the mean slopes between the two groups
kruskal.test(slopes$slope ~ slopes$covid)

## plot the slopes
pdf("figures/timepoint/OP_time_since_first_ICU_only.pdf")
ggplot(distances, aes(x = days_since_first, y = weighted_unifrac, color = covid)) +
    geom_point() + 
    stat_smooth(method = "lm") + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    ggtitle("Oropharyngeal swabs") +
    scale_color_manual(values = c("orangered2", "forestgreen")) +
    ylab("Weighted UniFrac Distance to First Sample") +
    xlab("Days Since First Sample")
dev.off()





