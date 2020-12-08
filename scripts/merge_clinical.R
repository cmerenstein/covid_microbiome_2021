library(dplyr)
library(tidyr)
library(stringr)

## R version 4.0.2
## Set working directory one above
setwd("../")

## Read in metadata for each sample. Only has basic data, must be merged with clinical data
meta = read.csv("data/raw/COVID_only_metadata.csv", stringsAsFactors = F)

## read in phenotype and clincial data
pheno = read.csv("data/raw/COVID_Metadata_11302020SW_v2.csv", quote = "\"",
                    stringsAsFactors = F, fileEncoding = "latin1")
pheno = pheno[,!(apply(pheno, 2, function(c){ return(all(is.na(c))) }))] ## some columns are all NA
pheno = pheno[pheno$Subject.ID != "",]
pheno$Subject.ID = gsub("-", "", pheno$Subject.ID)

## merge the two phenotype data tables
clinical = merge(meta, pheno, by.x = "SubjectID", by.y = "Subject.ID")
rownames(clinical) = clinical$SampleID

## get Phyla counts
counts = read.csv("data/raw/counts_taxa_level_1.csv", row.names = 1, stringsAsFactors = F)

## filter to match clinical, removes a few samples that were enrolled but subsequently deemed ineligible
counts = counts[rownames(clinical),]

## filter out samples with < 1000 reads
counts = counts[rowSums(counts) > 1000,]
counts = counts[,colSums(counts) > 0]

## filter clinical
clinical = clinical[rownames(counts),]

## It'll be useful to include visit as an integer, to calculate the first and last timepoint
clinical$visit_int = as.numeric(substr(clinical$Visit, 2, 2))

## there's one sample we need to remove from the clinical  cause it's misannotated
clinical = clinical[!(clinical$SampleType == "Endotracheal aspirate" & clinical$Intubated. == "no"),]

## 2 further samples that were excluded because of missannotated
## these both were not intubated, but had samples labled as ETA, which need to be removed
clinical = clinical[ !(clinical$SubjectID == "CORE189" & clinical$SampleType == "Endotracheal aspirate"),]
clinical = clinical[ !(clinical$SubjectID == "CORE264" & clinical$SampleType == "Endotracheal aspirate"),]

## add column to indicate if this sample is at or after date of intubation
clinical$Date = as.Date(clinical$Collection_date, format = "%m/%d/%Y")
clinical$tube_Date = as.Date(clinical$Intubation_date, format = "%m/%d/%Y")
clinical$after_intubation = clinical$Date >= clinical$tube_Date
clinical$after_intubation = ifelse(clinical$Intubation_date == "pre-existing tracheostomy",
                                    TRUE, clinical$after_intubation)

write.csv(clinical, "data/from_scripts/merged_clinical.csv")



