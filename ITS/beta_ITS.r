library("tidyverse")
library ("phyloseq")
library("vegan")
# library("RVAideMemoire")

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq.RData")


sampledata_meta_2$Sample = row.names(sampledata_meta_2)
physeq_for_beta = psmelt(physeq_for_beta)

physeq_for_beta_Species = physeq_for_beta

physeq_for_beta_Species = physeq_for_beta_Species %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq_for_beta_Species = pivot_wider(physeq_for_beta_Species, id_cols = Genus, names_from = "Sample", values_from = "Abundance")
physeq_for_beta_Species = as.data.frame(physeq_for_beta_Species)
row.names(physeq_for_beta_Species) = physeq_for_beta_Species[,1]
physeq_for_beta_Species = physeq_for_beta_Species[,-1]
physeq_for_beta_Species = t(physeq_for_beta_Species)

vdist = vegdist(physeq_for_beta_Species, method="bray")

vdist_NMDS = metaMDS(physeq_for_beta_Species, distance = "bray", k = 2, maxit = 999,  trymax = 500, wascores = TRUE)
NMDS = as.data.frame(scores(vdist_NMDS, display="site"))
NMDS$Sample = rownames(NMDS)

NMDS = left_join(NMDS, sampledata_meta_2, by ="Sample")

mycolors = c("#FF1202", "#401CE6", "#11a40e", "#ffc30f", "#e952fa")

# Figure 5d. NMDS ordination plot of Vitis communities (ITS)
Vitis = c("Vitis Elmer Swenson 2-7-13 cv. Prairie Star ",
"V. labrusca × V. riparia cv. Alfa",
"V. vinifera × V. amurensis cv. Adele ",
"V. riparia × V. vinifera cv. Mukuzani ",
"V. amurensis")

NMDS$Species = factor(NMDS$Species, levels = Vitis)
NMDS$Plant = factor(NMDS$Plant, levels = c(cultivars, V_amur))

ggplot(data = NMDS, aes(x = NMDS1, y = NMDS2, color = Species, shape = Plant)) + 
	geom_point(size=5) +
	scale_colour_manual(values = mycolors) +
	scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8)) +
	theme_bw() +
	guides(shape = guide_legend(order = 1), 
		color = guide_legend(order = 2)) +
	theme(legend.text = element_text(size = 15, face = "italic"),
		legend.title = element_text(size = 20)
	)
ggsave("Figure 5d. NMDS ordination plot of Vitis communities (ITS).png", width = 14, height = 9)
#


beta_Plant = adonis2(physeq_for_beta_Species ~ Species + Location + Plant, data = NMDS, permutations = 999, method = "bray")
as.data.frame(beta_Plant)
write.table(as.data.frame(beta_Plant), sep = "\t", quote = FALSE, file="Supporting information Table S6. PERMANOVA results (ITS data).txt")
