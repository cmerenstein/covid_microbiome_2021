library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)

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
ETA = meta[meta$SampleType == "Endotracheal aspirate" & meta$Study_group == "COVID" ,]

## one patient with pre-existing tube, we'll use days from admission for this individual
ETA$days_from_admission = as.integer(ETA$Date - as.Date(ETA$Admission.Date, format = "%m/%d/%Y"))
ETA$days_from_intubation = ifelse(is.na(ETA$days_from_intubation), ETA$days_from_admission, ETA$days_from_admission)
ETA$tube_Date[is.na(ETA$tube_Date)] = as.Date(ETA$Admission.Date, format = "%m/%d/%Y")[is.na(ETA$tube_Date)]

## read count data
counts = read.csv( "data/from_scripts/counts_without_contaminants_genus.csv", row.names = 1, stringsAsFactors = F)

## filter to >1000 reads
counts = counts[rowSums(counts) > 1000,]
ETA = ETA[rownames(ETA) %in% rownames(counts),]
counts = counts[rownames(ETA),]

## transfer to percents
percent = counts / rowSums(counts)

## specifically check staph
staph = data.frame(sample = rownames(percent), staph = (percent[,"g__Staphylococcus"] > .3), 
            patient = substr(rownames(percent), 1, 8))
length(unique(staph[staph$staph, "patient"]))
length(unique(staph$patient))


coryn = data.frame(sample = rownames(percent), coryn = (percent[,"g__Corynebacterium"] > .3), 
                    patient = substr(rownames(percent), 1, 8))
length(unique(coryn[coryn$coryn, "patient"]))
length(unique(coryn$patient))



## common taxa plus potential pathogens list from Dr. Collman
top_taxa = names(tail(sort(colMeans(percent)), n  = 10))
top_taxa = c(top_taxa, "g__Pseudomonas", "o__Enterobacterales_unassigned")


## STACKED BAR FIGURE
tidy = gather(data.frame(sample = rownames(percent), percent), "taxa", "percent", -sample)
tidy$taxa_other = ifelse(tidy$taxa %in% top_taxa, tidy$taxa, "other")
taxa_collapsed = group_by(tidy, sample, taxa_other) %>% summarize(percent = sum(percent)) %>%
                    ungroup() %>% as.data.frame()

## add days from intubation and patient
taxa_collapsed$days_from_intubation = ETA[taxa_collapsed$sample, "days_from_intubation"]
taxa_collapsed$patient = ETA[taxa_collapsed$sample, "SubjectID"]

## ALSO INCLUDE CLINICAL CULTURE RESULTS
cultures = read.csv("data/raw/clinical_culture.csv", stringsAsFactors = F, fileEncoding = "latin1")
cultures = cultures[cultures$Subject_ID %in% ETA$SubjectID,]
isolates = read.csv("data/raw/clinical_isolates_identified.csv", 
                            fileEncoding = "latin1", stringsAsFactors = F, header = F)$V1

isolates_present = lapply(isolates, function(isolate){
   present = sapply(cultures$Culture_notes, function(note){
        return(grepl(isolate, note))})
    return(unname(present))})
isolates_df = as.data.frame(do.call("cbind", isolates_present))
colnames(isolates_df) = isolates

cultures = as.data.frame(cbind(cultures, isolates_df))

## match to the intubation date
intubation_dates = unique(ETA[,c("SubjectID", "tube_Date")])
rownames(intubation_dates) = intubation_dates$SubjectID
cultures$intubation_date = as.Date(intubation_dates[cultures$Subject_ID,"tube_Date"])
cultures$Date = as.Date(cultures$Culture_Date, format = "%m/%d/%Y")
cultures$days_from_intubation = as.integer(cultures$Date - cultures$intubation_date)

## gather in tidy formatting for plotting, drop unusued columns
cultures_clean = cultures[, !(colnames(cultures) %in% c("Culture_Date", "Culture_notes",
                                                         "intubation_date", "Date"))]
cultures_tidy = gather(cultures_clean, "Clinical_Culture", "Present", -Subject_ID, -days_from_intubation)
cultures_tidy = cultures_tidy[cultures_tidy$Present,]
cultures_tidy$patient = cultures_tidy$Subject_ID ## the lazy man's rename()

## need to do some funky stuff for the plotting
## give a dummy y value to them, between 1 and 1.3, this makes it more visible when there are
## multiple positive isolates at one point
isolate_levels = c("Aspergillus niger", "Mouth Flora", "Staphylococcus epidermidis", "Enterococcus faecalis",
                    "Klebsiella Enterobacter aerogenes", "Streptococcus constellatus", "Escherichia coli",
                    "Staphylococcus aureus", "Streptococcus anginosus", "Pseudomonas putida",
                    "Yeast", "Stenotrophomonas maltophilia", "Enterobacter cloacae", 
                    "Pseudomonas aeruginosa", "Klebsiella pneumoniae", "No growth")
cultures_tidy$percent = 1 + (as.numeric(factor(cultures_tidy$Clinical_Culture, levels = isolate_levels)) / 26)


## ALSO include outcome
ETA$last_date = as.Date(ETA$Discharged.Date..Alive., format = "%m/%d/%Y")
ETA$last_date[is.na(ETA$last_date)] = as.Date(ETA$Deceased.Date, format = "%m/%d/%Y")[ETA$Outcome == "dead"]
patient_last_date = unique(ETA[,c("SubjectID", "last_date", "Outcome", "tube_Date")])
patient_last_date$days_from_intubation = as.integer(patient_last_date$last_date - as.Date(patient_last_date$tube_Date))

## rename for plotting
patient_last_date$Outcome = ifelse(patient_last_date$Outcome == "dead", "Dead", "Discharged")
patient_last_date$patient = patient_last_date$SubjectID

## set colors based on Dr. Collman's annotation of possible pathogens
taxa_collapsed$taxa_other = factor(taxa_collapsed$taxa_other, levels = c("g__Pseudomonas",
                            "g__Mycoplasma", "g__Stenotrophomonas", "g__Enterococcus", "g__Staphylococcus",
                            "f__Enterobacteriaceae_unassigned", "o__Enterobacterales_unassigned","g__Streptococcus" ,
                            "g__Alloprevotella", "g__Corynebacterium", "g__Fusobacterium", "g__Prevotella",
                            "other"))            

colors11 = c("#980d00", "#ff6e60", "#FF4500", "#FFA500", "#FFD700", "#c39bd3", "#bb00ff", "#ff6dc1",
             "#2980b9", "#aed6f1", "#196f3d", "#82e0aa", "grey")
            

pdf("figures/ETA/ETA_bar_plot.pdf", height = 10, width = 10)
ggplot(taxa_collapsed, aes(x = days_from_intubation)) + 
    theme_classic() + 
    geom_vline(data = patient_last_date, size = 1, linetype = "longdash", 
                        aes(xintercept = days_from_intubation, color = Outcome)) +
    geom_bar(stat = "identity", color = "black", width = 3, aes(y = percent, fill = taxa_other)) + 
    geom_point(data = cultures_tidy, fill = "grey", 
                        aes(y = percent, x = days_from_intubation, shape = Clinical_Culture)) +
    ylim(c(0, 1.6)) +
    facet_wrap(~ patient, ncol = 3) + 
    theme(legend.position = "right", legend.text = element_text(size = 8)) +
    scale_fill_manual(values = colors11)  +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 21, 7, 8, 9, 10, 11, 15, 16, 17)) +
    geom_hline(yintercept = 1.0) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
    theme(strip.text.x = element_text(size = 8))
dev.off()








