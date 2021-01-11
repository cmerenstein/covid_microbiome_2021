library(ape)
library(phyloseq)
library(plyr)

setwd("../../") ## root of this project

## build phylogeny from taxonomy
taxa = read.csv("data/SILVA/all_genus_tax_138.txt", stringsAsFactors = T)
taxa = taxa[!(taxa$genus %in% c("", "uncultured", "Incertae Sedis", "endosymbionts")) & taxa$kingdom == "Bacteria",]
taxa$genus = make.names(taxa$genus)  ## necessary for when they're used as column names
taxa = unique(taxa) ## this has to happen after make.names

## make as phylo class
taxa$genus = as.factor(taxa$genus) ## make.names coerced genus to a character
phylo = as.phylo.formula(~kingdom/phyla/class/order/family/genus, taxa)

## reads genus
covid = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", check.names = T, row.names = 1)
colnames(covid) = gsub("g__", "", colnames(covid))
covid = covid[grepl("COVID\\.", rownames(covid)), ] ## make sure there aren't any controls

## filter out things that aren't in the taxa tree
covid = covid[rowSums(covid) > 1000,
                    colnames(covid) %in% as.character(taxa$genus) & colSums(covid) > 10]

## load healthy controls
controls = read.csv("data/raw/454_healthy_controls.csv", row.names = 1)
colnames(controls) = gsub("g__", "", colnames(controls))
controls = controls[grepl("OP", rownames(controls)) | grepl("NP", rownames(controls)),]

## merge 2 data sets, need to add sample
combined = plyr::rbind.fill(data.frame(sample = rownames(covid), covid), 
                            data.frame(sample = rownames(controls), controls))

combined[is.na(combined)] = 0
rownames(combined) = combined$sample
combined = combined[,2:ncol(combined)]

## turn into percentages to make 454 data more comparable
combined_percent = combined / rowSums(combined)

## make a "tax table" for phyloseq object
taxa_strings = as.matrix(apply(taxa, 2, as.character))
rownames(taxa_strings) = taxa_strings[,"genus"]
tax_table_object = tax_table(taxa_strings)

## as phyloseq 
otu_object = otu_table(combined_percent, taxa_are_rows = F)
phyloseq_object = phyloseq(otu_object, phylo, tax_table_object)

## need to give placeholder edge lengths, of 0.01, to calculate unifrac
phy_tree(phyloseq_object)$edge.length = rep(0.01, nrow(phy_tree(phyloseq_object)$edge))

set.seed(87)
unifrac = UniFrac(phyloseq_object, weighted = F)
saveRDS(unifrac, "data/from_scripts/healthy_controls/unweighted_unifrac.rds")

weighted_unifrac = UniFrac(phyloseq_object, weighted = T)
saveRDS(weighted_unifrac, "data/from_scripts/healthy_controls/weighted_unifrac.rds")

## also try unweighted unifrac but with a threshold of 1%
phyloseq_object_threshold = phyloseq_object
otu_table(phyloseq_object_threshold)[otu_table(phyloseq_object_threshold) < 0.01] = 0
threshold_unifrac = UniFrac(phyloseq_object_threshold, weighted = F)
saveRDS(threshold_unifrac, "data/from_scripts/healthy_controls/threshold_1percent_unifrac.rds")


