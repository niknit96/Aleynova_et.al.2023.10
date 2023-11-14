library("tidyverse")
library ("phyloseq")
library("vegan")


dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./16s/physeq.RData")

sampledata_meta_2[sampledata_meta_2[] == "Vitis Elmer Swenson 2-7-13 cv. Prairie Star "] = "Pr-St"
sampledata_meta_2[sampledata_meta_2[] == "V. labrusca × V. riparia cv. Alfa"] = "Alfa"
sampledata_meta_2[sampledata_meta_2[] == "V. vinifera × V. amurensis cv. Adele "] = "Ad"
sampledata_meta_2[sampledata_meta_2[] == "V. riparia × V. vinifera cv. Mukuzani "] = "Muk"

# sampledata_meta_2$Species = paste(sampledata_meta_2$Species, sampledata_meta_2$Organ_material, sep=" ")

sampledata_meta_2$Sample = rownames(sampledata_meta_2)

# Figure 5a. Shannon′s alpha diversity boxplot (16s)
alpha_diversity = estimate_richness(physeq_for_alfa, split = TRUE, measures = NULL)
alpha_diversity$Sample = rownames(alpha_diversity)
alpha_diversity = left_join(sampledata_meta_2, alpha_diversity, by = "Sample")


ggplot(data = alpha_diversity, aes(x = Species, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	scale_x_discrete(limits = c(cultivars, "V. amurensis")) +
	# scale_x_discrete(limits = c(cultivars, 'V. amurensis Leaf', 'V. amurensis Stem')) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=30, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("Figure 5a. Shannon′s alpha diversity boxplot (16s).png", width = 8, height = 10)
#


# Pairwise Wilcoxon rank sum test of Vitis microbial communities (16s)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Species, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test of Vitis microbial communities (16s)")
Shannon
