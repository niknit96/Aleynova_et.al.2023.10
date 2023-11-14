library("tidyverse")
library ("phyloseq")
library("ComplexHeatmap")
library("circlize")


dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./16s/physeq.RData")

physeq_main = psmelt(physeq_main)

pre_upset_data = pivot_wider(physeq_main, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data[, names(pre_upset_data) %in% c(cultivars, V_amur)]

upset_data['V. amurensis'] = rowMeans(upset_data[,names(upset_data) %in% V_amur])
upset_data = upset_data[, names(upset_data) %in% c(cultivars, 'V. amurensis')]

# upset_data['V. amurensis Leaf'] = rowMeans(upset_data[,names(upset_data) %in% grep("Leaf", V_amur, value=TRUE)])
# upset_data['V. amurensis Stem'] = rowMeans(upset_data[,names(upset_data) %in% grep("Stem", V_amur, value=TRUE)])
# upset_data = upset_data[, names(upset_data) %in% c(cultivars, 'V. amurensis Leaf', 'V. amurensis Stem')]

# Figure 2. Heatmap for cultivars and V. amurensis (16s)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.1)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.1)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.1)] <- "<0.1%"

png("Figure 2. Heatmap for cultivars and V. amurensis (16s).png",  width = 17, height = 10, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	# column_title = "Top 10 most represented taxa of genus level for each V. amurensis plant",
	# column_title_gp = gpar(fontsize = 20, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = c(cultivars, 'V. amurensis'),
	# column_order = c(cultivars, 'V. amurensis Leaf', 'V. amurensis Stem'),
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 20),
	column_names_gp = gpar(fontsize = 23, fontface = "italic"),
	column_names_rot = -30,
	column_names_centered = FALSE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 12))
})

draw(Hmap)

dev.off()
#