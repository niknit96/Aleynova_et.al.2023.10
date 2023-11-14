library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

dir = as.character(read.table("dir.txt"))
setwd(dir)

asv_meta = read_qza("./16s/dada2/FeatureTable[Frequency]_16s.qza") # Table of metagenome data reads for each ASV by samples
asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)

tax_meta = read_qza("./16s/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_16s.qza") # Taxonomy table for metagenome data
tax_meta = parse_taxonomy(tax_meta$data, trim_extra=FALSE)
tax_meta[is.na(tax_meta)] <- "uncultured"

# Swap values ​​from Swap_unidentified in tax_meta to previous taxon level
Swap_unidentified = c("uncultured", "unidentified", "metagenome", "bacteriap25", "Unknown")
for(Swap in Swap_unidentified) {
    for(j in 1:7) {
        for (i in 1:nrow(tax_meta)) {
            if(grepl(Swap, tax_meta[i,j])) {
                tax_meta[i,j] = tax_meta[i,j-1] }
                else { tax_meta[i,j] = tax_meta[i,j]
            }       
        }
    }
}


#

# Filtering metagenome data from non-significant taxa
tax_meta = filter(tax_meta, 
Phylum != "d__Bacteria" & 
Phylum != "Unassigned" &
Kingdom != "d__Eukaryota" & 
Phylum != "k__Fungi" & 
Kingdom != "d__Archaea" &
Kingdom != "k__Viridiplantae" &
Genus != "g__Mitochondria" &
Genus != "g__Chloroplast")
#

tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")

tax_meta[grep("f__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("f__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("o__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("o__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("c__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("c__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("p__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("p__", tax_meta[,"Genus"]), "Genus"])
tax_meta[,"Genus"] = str_replace_all(tax_meta[,"Genus"], "(\\d+)$", function(x) as.numeric(x) + 1)

sampledata_meta = read.table(file="./Metagenome metadata/sampledata_16s_meta.txt", row.names = 2, header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome

sampledata_meta$Plant_organ = paste(sampledata_meta$Plant, sampledata_meta$Organ_material, sep=" ")

V_amur = c("Gh","M","S-Va","P-1","P-2","P-3","P-4","P-5","P-6","Kh-1","Kh-2")
cultivars = c("Pr-St","Alfa","Ad","Muk")

V_amur_organ = c("Gh Leaf","Gh Stem","M Leaf","M Stem","S-Va Leaf","S-Va Stem","P-1 Leaf","P-1 Stem",
"P-2 Leaf","P-2 Stem","P-3 Leaf","P-3 Stem","P-4 Leaf","P-4 Stem","P-5 Leaf","P-5 Stem",
"P-6 Leaf","P-6 Stem","Kh-1 Leaf","Kh-1 Stem","Kh-2 Leaf","Kh-2 Stem")
cultivars_organ = c("Pr-St Leaf","Pr-St Stem","Alfa Leaf","Alfa Stem",
"Ad Leaf","Ad Stem","Muk Leaf","Muk Stem")

sampledata_meta_2 = sampledata_meta

# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_meta)


tax_meta = tax_table(as.matrix(tax_meta))
asv_meta = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq = phyloseq(asv_meta, tax_meta, sampledata_meta)
#

physeq_read_counts = physeq
physeq_read_counts = merge_samples(physeq_read_counts, "Plant", fun = mean)
physeq_read_counts_organ = physeq
physeq_read_counts_organ = merge_samples(physeq_read_counts_organ, "Plant_organ", fun = mean)

physeq_for_div = physeq

physeq_organ = physeq

physeq = tax_glom(physeq, taxrank="Genus", NArm=TRUE) # Merge ASVs into genus level taxa

physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100) # Convert read counts to relative abundance (%)

physeq = merge_samples(physeq, "Plant", fun = mean) # Merge samples by plant
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)
physeq_for_info = physeq

physeq_main = filter_taxa(physeq, function(x) max(x) > 0.1, TRUE) # Genus level taxa with relative abundance > 0.1%

physeq_other = filter_taxa(physeq, function(x) max(x) < 0.1, TRUE) # Genus level taxa with relative abundance < 0.1%
physeq_other = merge_taxa(physeq_other, taxa_names(physeq_other))
taxa_names(physeq_other) = "Other"

# phyloseq object where genus level taxa with a relative abundance of less than 0.1% are placed in the "other" category
physeq = merge_phyloseq(physeq_main, physeq_other)
othertaxa = rbind(rank_names(tax_table(physeq)), "Other")
colnames(othertaxa) = othertaxa[1,]
othertaxa = othertaxa[-1,]
othertaxa = t(as.matrix(othertaxa))
row.names(othertaxa) = "Other"
tax_table(physeq)["Other",] <- tax_table(othertaxa)
#

###### Organ
physeq_organ = tax_glom(physeq_organ, taxrank="Genus", NArm=TRUE) # Merge ASVs into genus level taxa

physeq_organ = transform_sample_counts(physeq_organ, function(x) x / sum(x) * 100) # Convert read counts to relative abundance (%)

physeq_organ = merge_samples(physeq_organ, "Plant_organ", fun = mean) # Merge samples by plant
physeq_organ = transform_sample_counts(physeq_organ, function(x) x / sum(x) * 100)

physeq_organ_main = filter_taxa(physeq_organ, function(x) max(x) > 0.1, TRUE) # Genus level taxa with relative abundance > 0.1%

physeq_organ_other = filter_taxa(physeq_organ, function(x) max(x) < 0.1, TRUE) # Genus level taxa with relative abundance < 0.1%
physeq_organ_other = merge_taxa(physeq_organ_other, taxa_names(physeq_organ_other))
taxa_names(physeq_organ_other) = "Other"

# phyloseq object where genus level taxa with a relative abundance of less than 0.1% are placed in the "other" category
physeq_organ = merge_phyloseq(physeq_organ_main, physeq_organ_other)
othertaxa = rbind(rank_names(tax_table(physeq_organ)), "Other")
colnames(othertaxa) = othertaxa[1,]
othertaxa = othertaxa[-1,]
othertaxa = t(as.matrix(othertaxa))
row.names(othertaxa) = "Other"
tax_table(physeq_organ)["Other",] <- tax_table(othertaxa)
#
##########

# phyloseq object for beta diversity analysis
physeq_for_div = tax_glom(physeq_for_div, taxrank="Genus", NArm=TRUE)
physeq_for_beta = prune_taxa(taxa_names(physeq_main), physeq_for_div)
total = median(sample_sums(physeq_for_beta))
standf = function(x, t=total) round(t * (x / sum(x))) # for transformation for the even sample depth
physeq_for_beta = transform_sample_counts(physeq_for_beta, standf)
#

physeq_for_alfa = physeq_for_div # phyloseq object for alpha diversity analysis

# phyloseq object for network analysis
physeq_16s = physeq_for_div

## End



physeq_info_before_filtration_V_amur = filter_taxa(prune_samples(V_amur, physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_Pr_St = filter_taxa(prune_samples("Pr-St", physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_Alfa = filter_taxa(prune_samples("Alfa", physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_Ad = filter_taxa(prune_samples("Ad", physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_Muk = filter_taxa(prune_samples("Muk", physeq_for_info), function(x) max(x) > 0, TRUE)

physeq_info_after_filtration_V_amur = filter_taxa(prune_samples(V_amur, physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_Pr_St = filter_taxa(prune_samples("Pr-St", physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_Alfa = filter_taxa(prune_samples("Alfa", physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_Ad = filter_taxa(prune_samples("Ad", physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_Muk = filter_taxa(prune_samples("Muk", physeq_main), function(x) max(x) > 0, TRUE)


Taxa_of_genus_level = rbind(
paste("Taxa of genus level before filtration in V. amurensis samples:", length(taxa_names(physeq_info_before_filtration_V_amur))),
paste("Taxa of genus level before filtration in Vitis cv. Prairie Star samples:", length(taxa_names(physeq_info_before_filtration_Pr_St))),
paste("Taxa of genus level before filtration in V. labrusca × V. riparia cv. Alfa samples:", length(taxa_names(physeq_info_before_filtration_Alfa))),
paste("Taxa of genus level before filtration in V. vinifera × V. amurensis cv. Adele samples:", length(taxa_names(physeq_info_before_filtration_Ad))),
paste("Taxa of genus level before filtration in V. riparia × V. vinifera cv. Mukuzani:", length(taxa_names(physeq_info_before_filtration_Muk))),
paste("Taxa of genus level with relative abundance > 0.1% in V. amurensis samples:", length(taxa_names(physeq_info_after_filtration_V_amur))),
paste("Taxa of genus level with relative abundance > 0.1% in Vitis cv. Prairie Star samples:", length(taxa_names(physeq_info_after_filtration_Pr_St))),
paste("Taxa of genus level with relative abundance > 0.1% in V. labrusca × V. riparia cv. Alfa samples:", length(taxa_names(physeq_info_after_filtration_Alfa))),
paste("Taxa of genus level with relative abundance > 0.1% in V. vinifera × V. amurensis cv. Adele samples:", length(taxa_names(physeq_info_after_filtration_Ad))),
paste("Taxa of genus level with relative abundance > 0.1% in V. riparia × V. vinifera cv. Mukuzani:", length(taxa_names(physeq_info_after_filtration_Muk)))
)

write.table(Taxa_of_genus_level, sep = "\t", row.names = FALSE, quote = FALSE, file="./Taxa of genus level (16s).txt")


save(physeq, physeq_organ, physeq_read_counts_organ, physeq_for_div, V_amur, cultivars, V_amur_organ, cultivars_organ, physeq_for_beta, physeq_for_alfa, physeq_main, sampledata_meta_2, physeq_read_counts, dir,  file = "./16s/physeq.RData")

save(physeq_16s, file = "./16s/physeq_16s_for_network.RData")