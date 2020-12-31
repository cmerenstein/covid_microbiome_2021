library(ggplot2)
library(dplyr)
library(tidyr)

setwd("../")

## Read in merged metadata table.
meta = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## add viruses to metadata
viruses = read.table("data/raw/small_dna_viruses.tsv", sep = "\t", row.names = 1, header = T)
meta$rv_presence = viruses[rownames(meta), "rv_presence"]
meta$ttv_presence = viruses[rownames(meta), "ttv_presence"]

## get the first sample date
meta$Date = as.Date(meta$Collection_date, format = "%m/%d/%Y")
meta$first_collection = ave(meta$Date, meta$SubjectID, FUN = min)
meta$days_since_first_sample = as.numeric(meta$Date - meta$first_collection)

## read count data
counts = read.csv( "data/raw/counts_taxa_level_5.csv", row.names = 1, stringsAsFactors = F)
colnames(counts) = gsub("g__", "", colnames(counts))

## filter to >1000 reads
counts = counts[rowSums(counts) > 1000,]
meta = meta[rownames(meta) %in% rownames(counts),]
counts = counts[rownames(meta),]

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
                patient = meta[names(dominant_genus), "SubjectID"],
                visit = meta[names(dominant_genus), "Visit"],
                sample_type = meta[names(dominant_genus), "SampleType"],
                COVID = meta[names(dominant_genus), "Study_group"],
                percent = dg_percent,
                WHO_score = meta[names(dominant_genus), "Max.WHO.score"],
                days_since_first_sample = meta[names(dominant_genus), "days_since_first_sample"],
                Redondo = meta[names(dominant_genus), "rv_presence"],
                Anello = meta[names(dominant_genus), "ttv_presence"])

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
top_taxa = names(sort(table(dg$dominant_taxa_name), decreasing = T)[1:11]) ## idk why the figure has 11
dg$Most_Abundant_Taxon = ifelse(dg$dominant_taxa_name %in% top_taxa, dg$dominant_taxa_name, "other")

## filter and clean of dg
dg_filter = dg[dg$sample_type %in% c("Nasopharyngeal swab", "Oropharyngeal swab", "Endotracheal aspirate") &
               dg$COVID != "Control",]
dg_filter$visit = factor(dg_filter$visit, levels = sort(unique(dg_filter$visit)))
dg_filter$WHO_score = factor(dg_filter$WHO_score)
dg_filter$sample_type = factor(dg_filter$sample_type, levels = c("Nasopharyngeal swab", 
                                            "Oropharyngeal swab", "Endotracheal aspirate"))
dg_filter$patient = as.character(dg_filter$patient)

## make one variable for virus presence
dg_filter$viruses = ifelse(dg_filter$Redondo == "Positive" & dg_filter$Anello == "Positive", "Both",
                    ifelse(dg_filter$Redondo == "Positive" & dg_filter$Anello == "Negative", "Redondo",
                    ifelse(dg_filter$Redondo == "Negative" & dg_filter$Anello == "Positive", "Anello", "Neither")))

## we only want the points with viruses to get the outline
dg_filter$stroke = ifelse(dg_filter$viruses == "Neither", 0, 2)

## plot
pdf( "figures/genus/most_abundant_taxon_w_viruses.pdf", height = 15, width = 10)
dg_filter %>% arrange(visit) %>%
    ggplot(aes(x = days_since_first_sample, y = patient, size = percent, fill = Most_Abundant_Taxon)) + 
    theme_bw() + 
    geom_point(pch = 21, alpha = .75, aes(color = viruses, stroke = stroke)) +
    facet_grid(rows = vars(WHO_score), cols = vars(sample_type), scale = "free_y", space = "free_y") +
    scale_fill_manual(values = c("lightblue", "dodgerblue", "lightgreen", "forestgreen", 
                                "pink", "red", "peachpuff1", "firebrick", "orange", "purple4",
                                "orchid2", "yellow", "brown")) + 
    scale_color_manual(values = c("salmon", "purple", "black", "cyan")) +
    theme(strip.background = element_rect(fill = "white")) +
    theme(panel.grid.minor = element_blank())
dev.off()


## Make a second copy of the plot that is easier to see, but drops 2 points
#pdf(paste(folder, "figures/commensal_viruses/most_abundant_taxon_w_viruses_40_days.pdf", sep = ""), 
#            height = 15, width = 15)
#ggplot(dg_filter, aes(x = days_since_first_sample, y = patient, size = percent, fill = Most_Abundant_Taxon)) + 
#    theme_bw() + 
#    geom_point(pch = 21, aes(color = viruses, stroke = stroke)) +
#    facet_grid(rows = vars(WHO_score), cols = vars(sample_type), scale = "free_y", space = "free_y") +
#    scale_fill_manual(values = c("lightblue", "dodgerblue", "lightgreen", "forestgreen", 
#                                "pink", "red", "peachpuff1", "firebrick", "orange", "purple4",
#                                "orchid2", "yellow", "brown")) + 
#    scale_color_manual(values = c("salmon", "purple", "black", "cyan")) + 
#    xlim(c(0, 40))
#dev.off()












