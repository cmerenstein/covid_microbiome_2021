library(ape)
library(phyloseq)

setwd("../")

## need to read in metadata because of the few samples that got dropped
meta = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)

## build phylogeny from taxonomy
taxa = read.csv("data/SILVA/all_genus_tax_138.txt", stringsAsFactors = T)
taxa = taxa[!(taxa$genus %in% c("", "uncultured", "Incertae Sedis", "endosymbionts")) & taxa$kingdom == "Bacteria",]
taxa$genus = make.names(taxa$genus)  ## necessary for when they're used as column names
taxa = unique(taxa) ## this has to happen after make.names because dumbness abounds

## make as phylo class
taxa$genus = as.factor(taxa$genus) ## make.names coerced genus to a character
phylo = as.phylo.formula(~kingdom/phyla/class/order/family/genus, taxa)

## reads genus
genus = read.csv("data/from_scripts/counts_without_contaminants_genus.csv", check.names = T, row.names = 1)
colnames(genus) = gsub("g__", "", colnames(genus))
colnames(genus) = gsub("f__", "", colnames(genus))
colnames(genus) = gsub("_unassigned", "", colnames(genus))

## filter and make matrix
genus = genus[rownames(genus) %in% rownames(meta),]
genus_mat = as.matrix(genus)

## filter out things that aren't in the taxa tree
genus_mat = genus_mat[rowSums(genus_mat) > 1000,
                    colnames(genus_mat) %in% as.character(taxa$genus) & colSums(genus_mat) > 10]

## make a "tax table" for phyloseq object
taxa_strings = as.matrix(apply(taxa, 2, as.character))
rownames(taxa_strings) = taxa_strings[,"genus"]
tax_table_object = tax_table(taxa_strings)

## as phyloseq 
otu_object = otu_table(genus_mat, taxa_are_rows = F)
phyloseq_object = phyloseq(otu_object, phylo, tax_table_object)

## Calculating UniFrac based on taxonomic tree, rather than an alignment of our ASVs, which
## may be more accurate than estimating based off of short read sequencing. Each branch length
## gets the same weight at each taxonomic level
phy_tree(phyloseq_object)$edge.length = rep(0.01, nrow(phy_tree(phyloseq_object)$edge))

set.seed(87)
unifrac = UniFrac(phyloseq_object, weighted = F)
saveRDS(unifrac, "data/from_scripts/unweighted_unifrac.rds")

weighted_unifrac = UniFrac(phyloseq_object, weighted = T)
saveRDS(weighted_unifrac, "data/from_scripts/weighted_unifrac.rds")
