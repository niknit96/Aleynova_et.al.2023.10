library("tidyverse")
library ("phyloseq")
library ("RColorBrewer")
library("ComplexHeatmap")
library("circlize")
library("ggdendro")


colors_16s = c("#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fe4365", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#ecd078", "#d95b43", "#c02942", "#542437", "#53777a", 
"#4ecdc4", "#c7f464", "#ff6b6b", "#c44d58", "#774f38", "#e08e79", "#f1d4af", 
"#c5e0dc", "#018BFF", "#01FDFF", "#01FF61", "#F9FF01", "#FFD101", "#C08227", 
"#C02727", "#65216B", "#C177C7", "#EBBAE3", "#BAEBC4", "#397DD7", "#6FE8D9")

colors_its = c("#018BFF", "#01FDFF", "#01FF61", "#F9FF01", "#FFD101", "#C08227", 
"#C02727", "#65216B", "#C177C7", "#EBBAE3", "#BAEBC4", "#397DD7", "#6FE8D9", 
"#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fe4365", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#ecd078", "#d95b43", "#c02942", "#542437", "#53777a", 
"#4ecdc4", "#c7f464", "#ff6b6b", "#c44d58", "#774f38", "#e08e79", "#f1d4af", 
"#c5e0dc")


mycolors = colors_16s

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./16s/physeq.RData")

# For read counts above bars in bar plots
physeq_read_counts = psmelt(physeq_read_counts) %>%
	select(Class, Sample, Abundance) %>%
	group_by(., Class, Sample) %>%
	summarise(., Abundance = sum(Abundance))
For_barplot = pivot_wider(physeq_read_counts, id_cols = "Class", names_from = "Sample", values_from = "Abundance")
For_barplot = For_barplot[,-1]
For_barplot = colSums(For_barplot)
For_barplot = t(as.matrix(For_barplot))
row.names(For_barplot) = "Sum"
For_barplot = as.data.frame(For_barplot)
For_barplot = pivot_longer(For_barplot, cols = c(colnames(For_barplot)), names_to = "Sample", values_to = "Abundance")
colnames(For_barplot) = c("Sample", "Sum")
For_barplot[nrow(For_barplot)+1, ] = list("V. amurensis", sum(For_barplot[match(V_amur, For_barplot$Sample),"Sum"]))
For_barplot[,2] = apply(For_barplot[,2], 1, function(x) formatC(x, big.mark=",", format = "d"))
#

physeq = psmelt(physeq)
physeq = physeq %>%
	select(Class, Sample, Abundance) %>%
	group_by(., Class, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq = pivot_wider(physeq, id_cols = "Class", names_from = "Sample", values_from = "Abundance")
physeq = as.data.frame(physeq)
row.names(physeq) = physeq[,1]
physeq = physeq[,-1]

# For classes ordering and coloring in bar plots
physeq_factors = physeq
physeq_factors$Mean = rowMeans(physeq_factors[])
physeq_factors$Class = row.names(physeq_factors)
physeq_factors = physeq_factors[,c("Class", "Mean")]
physeq_factors = pivot_longer(physeq_factors, cols = !Class, names_to = "Sample", values_to = "Abundance")
physeq_factors = physeq_factors[order(physeq_factors$Abundance),]
physeq_factors = physeq_factors$Class
physeq_factors = physeq_factors[!physeq_factors == "Other"]
physeq_factors = c("Other", physeq_factors)
physeq_factors = cbind(physeq_factors, c("#C5C6D0", rev(mycolors[1:length(physeq_factors)-1])))
colnames(physeq_factors) = c("Class", "color")
#

Barplot = physeq[, names(physeq) %in% c(cultivars, V_amur)]

Barplot['V. amurensis'] = rowMeans(Barplot[,names(Barplot) %in% V_amur])
Barplot = Barplot[, c(cultivars, 'V. amurensis')]

Barplot$Class = row.names(Barplot)
Barplot = pivot_longer(Barplot, cols = !Class, names_to = "Sample", values_to = "Abundance")
Barplot = Barplot %>%
	filter(Abundance > 0)

# UPGMA clustering for x axis in bar plots
tree = Barplot %>%
	select(Class, Sample, Abundance) %>%
	filter(Sample %in% cultivars)
	tree = pivot_wider(tree, id_cols = "Sample", names_from = "Class", values_from = "Abundance")
	tree = as.data.frame(tree)
row.names(tree) = tree[,1]
tree = tree[,-1]
tree[is.na(tree[])] = 0
tree = hclust(dist(tree), method = "average")
tree = dendro_data(tree)
scale = 0.1
#

Above_barplot = filter(For_barplot, Sample %in% c(tree$labels$label, "V. amurensis"))
Above_barplot[] = Above_barplot[match(c(tree$labels$label, "V. amurensis"), Above_barplot$Sample),]
for(i in 1:nrow(Above_barplot)) {
	if(i %% 2 == 0){Above_barplot[i,"For_sum"] = 102}
	else {Above_barplot[i,"For_sum"] = 102}
}
Barplot = left_join(Barplot, Above_barplot, by = "Sample") 

Barplot$Sample = factor(Barplot$Sample, levels = c(tree$labels$label, "V. amurensis"))
Barplot_factors = physeq_factors[is.element(physeq_factors[,"Class"], unique(Barplot$Class)),]
Barplot$Class = factor(Barplot$Class, levels = Barplot_factors[, "Class"])

# Figure 1(a) Cultivars and V. amurensis 16s metagenome bar plot
ggplot(Barplot) +
	geom_col(aes(x = Sample, y = Abundance, fill = Class), colour = "black", linetype = "solid", size = 0.1) +
#	xlab("Sample ID") +
	scale_fill_manual(values = Barplot_factors[, "color"]) +
	guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
	geom_segment(
    	data = tree$segments,
    	aes(x = x, y = -y * scale, xend = xend, yend = -yend * scale)) +
	theme_bw() +
	theme(
	legend.text = element_text(size=30),
	legend.title = element_text(size=30),
	axis.text.x = element_text(size=40, angle=-30, vjust = 0.5, face = "italic"),
	axis.text.y = element_text(size=30),
	strip.text = element_text(size=30),
	axis.title.x = element_text(size=30),
	axis.title.y = element_text(size=30)
	) +
	labs(x = "Plant", y = "Relative abundance", fill = "Class") +
	scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    geom_text(aes(x = Sample, y = For_sum, label = Sum), inherit.aes = FALSE, size = 10) +
	geom_text(aes(x = Sample, y = Abundance, 
		label = ifelse(Abundance > 2, scales::percent(Abundance, scale = 1, accuracy = 1), ""), 
		group = Class), 
		inherit.aes = FALSE, position = position_stack(vjust = 0.5), 
		size = 10, col = "black")

ggsave("Figure 1(a) Cultivars and V. amurensis 16s metagenome bar plot.png", width = 25, height = 20)
#