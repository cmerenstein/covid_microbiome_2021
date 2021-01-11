library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(ape)
library(vegan)

setwd("../../") ## root of project

## get COVID study metadata
covid_meta = read.csv("data/from_scripts/merged_clinical.csv", stringsAsFactors = F, row.names = 1)

## get 454 control metadata
controls = read.csv("data/raw/454.healthy.csv", stringsAsFactors = F, row.names =1)
controls = controls[controls$respiratory_disease == "Healthy" & !(is.na(controls$respiratory_disease)),]
controls$Study_group = "Healthy"
controls$Max.WHO.score = NA

## combine metadata for sample types, covid/noncovid
sample_types = rbind( controls[,c("SampleType", "Study_group", "Max.WHO.score")],
                      covid_meta[,c("SampleType", "Study_group", "Max.WHO.score")])


## ------------ UNWEIGHTED 1% THRESHOLD ------------------------------------------
for(SampleType in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    # filter to sample type
    samples = rownames(sample_types)[sample_types$SampleType == SampleType]
    
#   get unifrac and filter to sample type
    unifrac = readRDS("data/from_scripts/healthy_controls/threshold_1percent_unifrac.rds")
    unifrac = as.matrix(unifrac)

    ## filter
    unifrac_samples = samples[samples %in% rownames(unifrac)]
    unifrac = unifrac[unifrac_samples, unifrac_samples]

    ordinates = pcoa(unifrac)
    vectors = ordinates$vectors
    pcoa_df = data.frame(samples = rownames(vectors), Axis.1 = vectors[,"Axis.1"],
                        Axis.2 = vectors[,"Axis.2"], 
                        group = sample_types[rownames(vectors),"Study_group"], 
                        Max.WHO.score = sample_types[rownames(vectors), "Max.WHO.score"])
    
    ## define centroids for plotting
    centroids = aggregate(cbind(Axis.1, Axis.2) ~ group, data = pcoa_df, mean)

    ## plot
    SampleType = gsub(" ", "_", SampleType)
    pdf(paste("figures/454_healthy_controls/", SampleType, "COVID_status_unweighted.pdf", sep = ""))
    plot = ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = group)) + 
        theme_classic() +
        geom_point(alpha = .5, size = 4) +
        scale_color_manual(values = c("orangered2", "navy", "forestgreen")) +
        geom_point(data = centroids, shape = 4, size = 4, stroke = 3,
                                     aes(x = Axis.1, y = Axis.2, color = group)) +
        scale_color_manual(values = c("orangered2", "navy", "forestgreen")) 
    print(plot)
    dev.off()
}

## ---------------------- PERMANOVA via ADONIS ----------------------------------------
## make function for pairwise comparisons
pairwise_adonis = function(unifrac, grouping){
    groupings = unique(grouping)
    tested = character() ## so we don't repeat pairs

    ## we want all non-redundant pairs of tests
    pvals_list_1 = lapply(groupings, function(grouping_1){

        tested <<- c(tested, grouping_1) ## slow but clear, prevents repeated pairs

        other_groups = groupings[!(groupings %in% tested)]
        pvals_list_2 = lapply(other_groups, function(grouping_2){

            ## filter to only the 2 we're comparing
            unifrac_groups = unifrac[grouping %in% c(grouping_1, grouping_2),
                                     grouping %in% c(grouping_1, grouping_2)]
            groupings_filtered = grouping[grouping %in% c(grouping_1, grouping_2)]

            permanova = adonis(as.dist(unifrac_groups) ~ groupings_filtered, permutations = 10000)
            pval = permanova$aov.tab[1, "Pr(>F)"]
            return(data.frame(pair = paste(grouping_1, grouping_2, sep = " | "), pval = pval))
        })
        ## merge pvals into data frame before returning
        return(do.call("rbind", pvals_list_2))
    })
    return(do.call("rbind", pvals_list_1))
}


## do one sample type at a time, UNWEIGHTED
for (sample_type in c("Nasopharyngeal swab", "Oropharyngeal swab")){
    unifrac = readRDS("data/from_scripts/healthy_controls/threshold_1percent_unifrac.rds")
    unifrac = as.matrix(unifrac)

    ## filter
    sample_types = sample_types[rownames(sample_types) %in% rownames(unifrac),]

    ## get sample type
    grouping_filtered = sample_types[sample_types$SampleType == sample_type,]
    unifrac_filtered = unifrac[rownames(grouping_filtered), rownames(grouping_filtered)]

    ## pairwise permanova for each pair of gropus
    permanova = pairwise_adonis(unifrac_filtered, grouping_filtered$Study_group)

    print(sample_type)
    print(permanova)
}







