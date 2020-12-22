library(ggplot2)
library(dplyr)
library(stringr)
library(vegan)
library(tidyr)

setwd("../")

## viral load data FILTER TO ONLY OP NP
viral = read.csv("data/raw/SARS_CoV2_viral_load.csv")
viral = viral[viral$sample_type %in% c("OP", "NP", "NP-OP"),]

## clinical data
clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)
abvs = list("Nasopharyngeal swab" = "NP", "Oropharyngeal swab" = "OP", "Endotracheal aspirate" = "ETA")
clinical$sample_type_abbreviation = unlist(abvs[clinical$SampleType])

## WE ONLY WANT OP NP HERE
OP_NP = clinical[grepl("swab", clinical$SampleType),]

## get ID that can match between viral and OP_NP
OP_NP$subject_date = paste(substr(OP_NP$SubjectID, 5, 8), OP_NP$Collection_date, sep = "_")

#viral$sample_date = str_pad(viral$sample_date, 10, pad = "0")
viral$subject_date = paste( str_pad(viral$patient_id, width = 4, side = "left", pad = "0"),
                            viral$sample_date, sep = "_")

## virus with microbiome
virus_with_microbiome = viral[viral$subject_date %in% OP_NP$subject_date,]

## when there are multiple samples for a day, combine, take max
tidy = gather(virus_with_microbiome, "test", "copies_per_ul", `N1_copies_per_ul`, `N2_copies_per_ul`)
tidy$copies_per_ul = as.numeric(ifelse(tidy$copies_per_ul == "undetermined", 0, 
                                ifelse(tidy$copies_per_ul == "", 0, tidy$copies_per_ul)))
max_virus = group_by(tidy, subject_date) %>% summarize(max_copies_per_ul = max(copies_per_ul)) %>% 
    ungroup() %>% as.data.frame()

## merge the ones we have
OP_NP_filtered = OP_NP[OP_NP$subject_date %in% max_virus$subject_date,]
OP_NP_filtered = merge(OP_NP_filtered, max_virus[, c("subject_date", "max_copies_per_ul")],
                             by = "subject_date", all.x = T, all.y = F)

## OP swabs have been most informative so we're just going to use that for microbiome diversity
OP = OP_NP_filtered[OP_NP_filtered$SampleType == "Oropharyngeal swab",]
rownames(OP) = OP$SampleID

## get counts at the GENUS level
counts = read.csv("data/raw/counts_taxa_level_5.csv", row.names = 1, stringsAsFactors = F)

# order counts the same as OPdata
counts = counts[rownames(OP),]

## FILTER TO > 1000 READS
counts = counts[rowSums(counts) > 1000,]
OP = OP[rownames(counts),]

## get alpha diversity for each sample
simpson = diversity(counts, "simpson")

## correlation with viral load
cor.test(x = OP$max_copies_per_ul, y = simpson, method = "spearman")

## plot and plot dominant taxa
dominant_taxa = apply(counts, 1, function(r){return( colnames(counts)[r == max(r)])})

## have to put jitter in before taking the log so it only really effects the 0 samples 
data.frame(dominant_taxa, simpson_diversity = simpson, 
                log_viral_copies_w_jitter = log10(OP$max_copies_per_ul + 1 + rnorm(nrow(OP), 0, .3))) %>%
    ggplot(aes(x = log_viral_copies_w_jitter, y = simpson, color = dominant_taxa)) + 
    theme_bw() + 
    geom_point(size = 4) + 
    scale_color_manual(values = c("skyblue", "cyan", "navy", "lightgreen", "forestgreen", "salmon", "red",
                                   "purple1", "plum4", "orange3", "gold", "chocolate", "grey")) +
    theme(text = element_text(size = 20)) + 
    ggtitle("Oropharyngeal swabs")


## ---------- NP SWABS NOW -----------------------
NP = OP_NP_filtered[OP_NP_filtered$SampleType == "Nasopharyngeal swab",]
rownames(NP) = NP$SampleID

## get counts at the GENUS level
counts = read.csv("data/raw/counts_taxa_level_5.csv", row.names = 1, stringsAsFactors = F)

# order counts the same as NPdata
counts = counts[rownames(NP),]

## FILTER TO > 1000 READS
counts = counts[rowSums(counts) > 1000,]
NP = NP[rownames(counts),]

## get alpha diversity for each sample
simpson = diversity(counts, "simpson")

## correlation with viral load
cor.test(x = NP$max_copies_per_ul, y = simpson, method = "spearman")

## plot and plot dominant taxa
dominant_taxa = apply(counts, 1, function(r){return( colnames(counts)[r == max(r)])})


## have to put jitter in before taking the log so it only really effects the 0 samples
data.frame(dominant_taxa, simpson_diversity = simpson,
                log_viral_copies_w_jitter = log10(NP$max_copies_per_ul + 1 + rnorm(nrow(NP), 0, .3))) %>%
    ggplot(aes(x = log_viral_copies_w_jitter, y = simpson, color = dominant_taxa)) +
    theme_bw() +
    geom_point(size = 4) +
    scale_color_manual(values = c("cyan", "skyblue", "navy", "lightgreen", "forestgreen", "salmon", "red",
                                   "purple1", "plum4", "orange3", "gold", "chocolate","yellow", "linen", "grey")) +
    theme(text = element_text(size = 20)) +
    ggtitle("Nasopharyngeal swabs")






