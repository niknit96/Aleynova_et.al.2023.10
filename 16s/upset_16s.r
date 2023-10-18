library ("RColorBrewer")
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
# V_amur = V_amur[V_amur %in% colnames(upset_data)]

upset_data['V. amurensis'] = rowMeans(upset_data[,names(upset_data) %in% V_amur])
upset_data = upset_data[, names(upset_data) %in% c(cultivars, 'V. amurensis')]

upset_data[] = ifelse(upset_data[] > 0, 1, 0) # Transform data to upset format
upset_data = upset_data[!apply(upset_data == 0, 1, all),]

# A table that shows which plants contain taxa from Upset plots (for V. amurensis)
upset_data_for_suppl = upset_data[, c(cultivars, 'V. amurensis')]
upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)
for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]
upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}
Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Supporting information 3. Intersections in 16s metagenome data of cultivars and V. amurensis plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 1b. Upset diagram for V. amurensis plants (16s)
png("Figure 1b. Upset diagram for cultivars and V. amurensis plants (16s).png",  width = 13, height = 8, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 20, fontface="italic"),
	set_order = c(cultivars, 'V. amurensis'), 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()


#######