library(ggplot2)
library(dplyr)
library(tidyr)

setwd("../")

## read in feature table
ft = t(read.table("data/raw/reformated_feature_table.txt", row.names = 1, header = T, sep = "\t"))
ft = ft[rowSums(ft) > 1000,]

percent = ft / rowSums(ft)

meta = read.csv("data/raw/metadata_with_dna_conc.csv", row.names = 1, stringsAsFactors = F)
meta = meta[rownames(percent),]

ASV_to_test = colnames(percent)[colMeans(percent) > 0.0001]
correlation_list = lapply(ASV_to_test, function(ASV) {

    correlation = cor.test(y = log(percent[, ASV] + .0001), x = log(meta[rownames(percent), "library_conc_ng_ul"] + .0001),
                            method = "pearson")
    return(data.frame(ASV = ASV, mean_abundance = mean(percent[,ASV]), r = correlation$estimate, p = correlation$p.value))
})
correlation = do.call("rbind", correlation_list)
correlation$adjusted_p = p.adjust(correlation$p, "bonferroni")

## write contaminants
write.csv(correlation, "data/from_scripts/contaminants_ASV.csv", row.names = F)


tidy = gather(data.frame(sample = rownames(percent), percent[, ASV_to_test]), "ASV", "percent", -sample)
tidy$ASV = gsub("X", "", tidy$ASV)
tidy$dna_conc = meta[tidy$sample, "library_conc_ng_ul"]

pdf("figures/contamination_ASV.pdf", height = 20, width = 8)
tidy[tidy$ASV %in% correlation[correlation$adjusted_p < 0.05 & correlation$r < 0, "ASV"], ] %>%
    ggplot(aes(x = log(0.001 + dna_conc), y = (0.001 + percent))) + 
    theme_classic() + 
    geom_point() + 
    facet_wrap(~ASV, ncol = 3)  
dev.off()


