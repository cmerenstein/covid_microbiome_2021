library(tidyverse)
library(polyafit) 

setwd("../")

## Load ASV feature table
asv = t(read.table("data/raw/reformated_feature_table.txt",
                sep = "\t", row.names = 1, header = T))

## filter out likely contaminant ASVs
contaminants = read.csv("data/from_scripts/contaminants_ASV.csv")
likely_contaminants = contaminants[contaminants$adjusted_p < 0.1 & contaminants$r < 0, "ASV"]
asv = asv[ ,!(colnames(asv) %in% likely_contaminants)]

## Load clinical table
clinical = read.csv("data/from_scripts/merged_clinical.csv", row.names = 1, stringsAsFactors = F)

## get just ETA and OP samples
clinical_filter = clinical[clinical$SampleType %in% c("Oropharyngeal swab", "Endotracheal aspirate") &
                                                                    clinical$Study_group == "COVID",]

## add subject-time pair so we can match OP and ETA later
clinical_filter$Subject_Sample = paste(clinical_filter$SubjectID, clinical_filter$visit_int, sep = "_")

## filter to only these subject-sample points
ETA = clinical_filter[clinical_filter$SampleType == "Endotracheal aspirate",]
ETA_OP = clinical_filter[clinical_filter$Subject_Sample %in% ETA$Subject_Sample, ]

## filter to only those we have pairs for
samples = ETA_OP %>% group_by(Subject_Sample) %>% summarize(n = n()) %>%
    ungroup() %>% as.data.frame()
rownames(samples) = samples$Subject_Sample
pairs = samples[samples$n == 2,]

ETA_OP = ETA_OP[ ETA_OP$Subject_Sample %in% pairs$Subject_Sample,]

## filter counts to just these samples
ETA_OP_counts = asv[rownames(ETA_OP),]

## get p values for each taxa for each sample
pvalues_list = list()
for (Subject_Sample in unique(ETA_OP$Subject_Sample)){

    ## get corresponding OP and ETA samples
    ETA_sample = ETA_OP[ETA_OP$Subject_Sample == Subject_Sample &
                        ETA_OP$SampleType == "Endotracheal aspirate", "SampleID"]
    OP_sample = ETA_OP[ETA_OP$Subject_Sample == Subject_Sample &
                        ETA_OP$SampleType == "Oropharyngeal swab", "SampleID"] 
    
    ## get counts, filter to at least 5 counts
    ETA_counts = ETA_OP_counts[ETA_sample,]
    OP_counts = ETA_OP_counts[OP_sample,]

    filter = (ETA_counts + OP_counts) >= 5
    ETA_counts = ETA_counts[filter]
    OP_counts = OP_counts[filter]

    ## generate p values using the optim_polya function
    pfit = optim_polya(rbind(ETA_counts, OP_counts))
    pvals = 1 - ppolya_marginal(ETA_counts, pfit$par, log.p = FALSE)
    
    ## save p values and asv names
    pvals_df = data.frame(Subject_Sample, asv = names(ETA_counts), pvals, theta = (1 / sum(pfit$par)))
    pvalues_list[[Subject_Sample]] <- pvals_df
}

## bind, FDR adjust, write
pvalues_df = do.call("rbind", pvalues_list)
pvalues_df$FDR = p.adjust(pvalues_df$pvals, method = "BH")
write.csv(pvalues_df, "data/from_scripts/pvalues_ASVs.csv", row.names = F)

## get sample names for pairs
sample_names = merge(ETA_OP[ETA_OP$SampleType == "Endotracheal aspirate", c("Subject_Sample", "SampleID")],
                     ETA_OP[ETA_OP$SampleType == "Oropharyngeal swab", c("Subject_Sample", "SampleID")],
                     by = "Subject_Sample")
colnames(sample_names) = c("Subject_Sample", "ETA_Sample", "OP_Sample")

## get percentages
percent = asv[unique(rownames(ETA_OP)),] / rowSums(asv[unique(rownames(ETA_OP)),])
percent_df = data.frame(sample = rownames(percent), percent)
percent_tidy = gather(percent_df, "asv", "percent", -sample)

## -------------------- plot and label outliers --------------
asv_percentages = merge(pvalues_df, sample_names, by = "Subject_Sample")

## include log of theta
asv_percentages$theta = round(log10(asv_percentages$theta), digits = 2)
asv_percentages$Subject_Sample_theta = paste(asv_percentages$Subject_Sample, " (",
                                                 asv_percentages$theta, ")", sep = "")

## ETA
asv_percentages = merge(asv_percentages, percent_tidy, by.x = c("ETA_Sample", "asv"), by.y = c("sample", "asv"))
asv_percentages = rename(asv_percentages, ETA_percent = percent)

## OP
asv_percentages = merge(asv_percentages, percent_tidy, by.x = c("OP_Sample", "asv"), by.y = c("sample", "asv"))
asv_percentages = rename(asv_percentages, OP_percent = percent)

## add taxonomic information
taxonomy = read.csv("data/raw/taxonomy.tsv", sep = "\t")
taxonomy = separate(taxonomy, "Taxon", into = c("domain", "phylum", "class", "order",
                        "family", "genus"), sep = "; ")
rownames(taxonomy) = taxonomy$Feature.ID

asv_percentages$genus = taxonomy[asv_percentages$asv, "genus"]
asv_percentages$family = taxonomy[asv_percentages$asv, "family"]
asv_percentages$order = taxonomy[asv_percentages$asv, "order"]

## add outlier column, for plotting, just classify at the highest level available, the lazy ifelse way
asv_percentages$outlier = ifelse(asv_percentages$FDR < 0.1 & 
                            asv_percentages$ETA_percent > asv_percentages$OP_percent,
                              ifelse( !(is.na(asv_percentages$genus)), asv_percentages$genus, 
                                ifelse( !(is.na(asv_percentages$family)), asv_percentages$family,
                                                                            asv_percentages$order)), "")
    
## plot
pdf("figures/OP_ETA/ASV_outlier_plot.pdf", height = 15, width = 10)
ggplot(asv_percentages, aes(x = (0.001 + OP_percent), y = (0.001 + ETA_percent))) + 
    theme_classic() + 
    geom_point(shape = 21, size = 3, aes(fill = outlier)) + 
#    geom_text(aes(label = outlier), nudge_x = .6, nudge_y = .1, size = 3) +
    scale_y_log10() + 
    scale_x_log10() + 
    facet_wrap(~Subject_Sample, ncol = 3) + 
    ylab("ETA") + 
    xlab("OP") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_fill_manual(values = c("grey", "#c39bd3", "orange", "#aed6f1", "navy", "#82e0aa", 
                                     "#FFD700", "pink", "cyan"))
dev.off()


























