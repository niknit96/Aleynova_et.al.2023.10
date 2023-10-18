library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

dir = as.character(read.table("dir.txt"))
setwd(dir)

asv_meta = read_qza("./ITS/dada2/FeatureTable[Frequency]_ITS.qza") # Table of metagenome data reads for each ASV by samples
asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)

tax_meta = read_qza("./ITS/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_ITS.qza") # Taxonomy table for metagenome data
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
Kingdom != "k__Rhizaria" &
Kingdom != "k__Protista" &
Kingdom != "k__Metazoa" &
Kingdom != "k__Alveolata" &
Genus != "g__Chloroplast")
#

tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")

sampledata_meta = read.table(file="./Metagenome metadata/sampledata_its_meta.txt", row.names = 1, header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome

V_amur = c("Gh","M","S-Va","P-1","P-2","P-3","P-4","P-5","P-6","Kh-1","Kh-2")
cultivars = c("Pr-St","Alfa","Ad","Muk")

sampledata_meta_2 = sampledata_meta

# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_meta)
tax_meta = tax_table(as.matrix(tax_meta))
asv_meta = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq = phyloseq(asv_meta, tax_meta, sampledata_meta)
#

physeq_read_counts = physeq
physeq_read_counts = merge_samples(physeq_read_counts, "Plant", fun = mean)

physeq_for_div = physeq

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

# phyloseq object for beta diversity analysis
physeq_for_div = tax_glom(physeq_for_div, taxrank="Genus", NArm=TRUE)
physeq_for_beta = prune_taxa(taxa_names(physeq_main), physeq_for_div)
total = median(sample_sums(physeq_for_beta))
standf = function(x, t=total) round(t * (x / sum(x))) # for transformation for the even sample depth
physeq_for_beta = transform_sample_counts(physeq_for_beta, standf)
#

physeq_for_alfa = physeq_for_div # phyloseq object for alpha diversity analysis


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

write.table(Taxa_of_genus_level, sep = "\t", row.names = FALSE, quote = FALSE, file="./Taxa of genus level (ITS).txt")


save(physeq, physeq_for_div, V_amur, cultivars, physeq_for_beta, physeq_for_alfa, physeq_main, sampledata_meta_2, physeq_read_counts, dir,  file = "./ITS/physeq.RData")