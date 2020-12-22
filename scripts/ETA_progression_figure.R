library(ggplot2)
library(dplyr)
library(tidyr)

setwd("../")

## Read in merged metadata table.
meta = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)
meta$Date = as.Date(meta$Date)
meta$tube_Date = as.Date(meta$tube_Date)
meta$symptom_onset = as.Date(meta$Symptom.Onset.Date, format = "%m/%d/%Y")

## add days from intubation and days from symptom onset
meta$days_from_intubation = as.integer(meta$Date - meta$tube_Date)
meta$days_from_SO = as.integer(meta$Date - meta$symptom_onset)

## Filter to just ETA samples, just COVID
## Removing one sample with preexisting intubation that has no intubation date
ETA = meta[meta$SampleType == "Endotracheal aspirate" & meta$Study_group == "COVID" & !(is.na(meta$days_from_intubation)),]

## read count data
counts = read.csv( "data/raw/counts_taxa_level_5.csv", row.names = 1, stringsAsFactors = F)
colnames(counts) = gsub("g__", "", colnames(counts))

## filter to >1000 reads
counts = counts[rowSums(counts) > 1000,]
ETA = ETA[rownames(ETA) %in% rownames(counts),]
counts = counts[rownames(ETA),]

## dominant taxa for each sample
dominant_genus = apply(counts, 1, function(r){
    return(colnames(counts)[which(r == max(r))])
}) %>% unlist()

## the perent that the dominant taxa makes up in each sample
dg_percent = apply(counts, 1, function(r){
    return(max(r) / sum(r))
}) %>% unlist()

## how many have a dominant taxa above 50%?
sum(dg_percent > .5) / length(dg_percent)

## combine for plotting
dg = data.frame(sample = names(dominant_genus), dominant_genus,
                patient = ETA[names(dominant_genus), "SubjectID"],
                sample_type = ETA[names(dominant_genus), "SampleType"],
                COVID = ETA[names(dominant_genus), "Study_group"],
                percent = dg_percent,
                WHO_score = ETA[names(dominant_genus), "Max.WHO.score"],
                days_from_intubation = ETA[names(dominant_genus), "days_from_intubation"])

## Get taxa tree to match phyla to genera, for more informative plotting
taxa = read.csv("data/SILVA/all_genus_tax_138.txt", stringsAsFactors = T) ## originally downloaded with Kraken2
taxa = taxa[!(taxa$genus %in% c("", "uncultured", "Incertae Sedis", "endosymbionts")) & 
                taxa$kingdom == "Bacteria",]
taxa$genus = make.names(taxa$genus)  ## necessary for when they're used as column names
taxa = unique(taxa) ## this has to happen after make.names because dumbness abounds
rownames(taxa) = taxa$genus

## add phyla to the genus names
dg$dominant_taxa_name = paste(taxa[dg$dominant_genus,"phyla"], dg$dominant_genus, sep = ": ")

## find most common taxa, then make "other" catagory
top_taxa = names(sort(table(dg$dominant_taxa_name), decreasing = T)[1:8]) ## idk why the figure has 11
dg$Most_Abundant_Taxon = ifelse(dg$dominant_taxa_name %in% top_taxa, dg$dominant_taxa_name, "other")

## clean of dg
dg$WHO_score = factor(dg$WHO_score)
dg$patient = as.character(dg$patient)


## plot
pdf( "figures/genus/ETA_samples_dominant_genus.pdf", height = 7, width = 7)
ggplot(dg, aes(x = days_from_intubation, y = percent, fill = Most_Abundant_Taxon)) + 
    theme_bw() + 
    geom_point(pch = 21, size =4) +
    geom_line(aes(group = patient), linetype = "dashed") +
    facet_wrap(~WHO_score, ncol = 1) +
    scale_fill_manual(values = c("lightblue", "dodgerblue", "lightgreen", "forestgreen", 
                                "pink", "red", "peachpuff1", "firebrick", "orange", "purple4",
                                "orchid2", "yellow", "brown"))  +
    theme(legend.text = element_text(size = 10))
dev.off()












