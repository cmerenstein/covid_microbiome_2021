library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

setwd("../")

# get sample processing data
meta = read.csv("data/raw/metadata_with_dna_conc.csv", row.names = 1, stringsAsFactors = F)

## counts at the genus level
counts = read.csv("data/raw/counts_taxa_level_5.csv", row.names = 1, stringsAsFactors = F)
counts = counts[rowSums(counts) > 1000,]

## match rows
meta = meta[rownames(counts),]

## convert to percents, put in data frame
percent = t(apply(counts, 1, function(r){return( r / sum(r))}))

## identify taxa associated with DNA conc
genus_to_test = colnames(percent)[colMeans(percent) > 0.001]
correlation_list = lapply(genus_to_test, function(genus) {
        
    ## remove zeros
    dna_conc = meta[rownames(percent), "library_conc_ng_ul"]
    genus_percent = percent[,genus]
    dna_conc = dna_conc[genus_percent != 0]
    genus_percent = genus_percent[genus_percent != 0]

    correlation = cor.test(y = genus_percent, x = dna_conc, method = "pearson")
    return(data.frame(genus = genus, rho = correlation$estimate, p = correlation$p.value))
})
correlation = do.call("rbind", correlation_list)
correlation$fdr = p.adjust(correlation$p, "BH")

tidy = gather(data.frame(sample = rownames(percent), percent[, genus_to_test]), "genus", "percent", -sample)
tidy$genus = gsub("X", "", tidy$genus)
tidy$dna_conc = meta[tidy$sample, "library_conc_ng_ul"]
tidy$sample_type = meta[tidy$sample, "sample_type"]

## remove 0s
tidy = tidy[tidy$percent != 0,]

pdf("figures/contamination.pdf", height = 15, width = 8)
tidy[tidy$genus %in% correlation[correlation$fdr < 0.05 & correlation$rho < 0, "genus"], ] %>%
    ggplot(aes(x = dna_conc, y = percent, color = sample_type)) +
    theme_classic() +
    geom_point() +
    facet_wrap(~genus, ncol = 3) +
    scale_y_log10() + 
    scale_x_log10()
dev.off()


## using these figures, we select contaminants
contaminants = c("g__Ralstonia", "g__Sphingomonas", "g__Chloroplast", "g__Massilia",
                 "g__Methylobacterium.Methylorubrum", "g__Aenigmarchaeales", "g__Acinetobacter",
                 "f__uncultured_uncultured", "f__Sphingomonadaceae_unassigned", "f__Comamonadaceae_unassigned")

write.table(contaminants, "data/from_scripts/contaminant_genera.csv", row.names = F, col.names = F)

read.table("data/from_scripts/contaminant_genera.csv")$V1
counts_without_contaminants = counts[, !(colnames(counts) %in% contaminants)]

write.csv(counts_without_contaminants, "data/from_scripts/counts_without_contaminants_genus.csv")
