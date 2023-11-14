library("tidyverse")
library ("phyloseq")
library("NetCoMi")
library("ggraph")
library("tidygraph")
library("SpiecEasi")

dir = as.character(read.table("dir.txt"))
setwd(dir)

load("./16s/physeq_16s_for_network.RData")
load("./ITS/physeq_its_for_network.RData")
# load("./16s/physeq.RData")
# load("./ITS/physeq.RData")

sampledata = read.table(file="./Metagenome metadata/sampledata_16s_meta.txt", row.names = 1, header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome



physeq_16s = merge_samples(physeq_16s, "Sample", fun = mean)
physeq_16s = tax_glom(physeq_16s, taxrank="Genus", NArm=TRUE)
physeq_16s_filt = filter_taxa(physeq_16s, function(x) sum(x > 3) > (0.6*length(x)), TRUE)
physeq_16s = prune_taxa(taxa_names(physeq_16s_filt), physeq_16s)
sample_data(physeq_16s) = sample_data(sampledata)

physeq_its = merge_samples(physeq_its, "Sample", fun = mean)
physeq_its = tax_glom(physeq_its, taxrank="Genus", NArm=TRUE)
physeq_its_filt = filter_taxa(physeq_its, function(x) sum(x > 3) > (0.6*length(x)), TRUE)
physeq_its = prune_taxa(taxa_names(physeq_its_filt), physeq_its)
sample_data(physeq_its) = sample_data(sampledata)

taxa_names(physeq_16s) = paste0("bacteria-", taxa_names(physeq_16s))
taxa_names(physeq_its) = paste0("fungi-", taxa_names(physeq_its))

tax = as.data.frame(rbind(tax_table(physeq_16s), tax_table(physeq_its)))


physeq = multi.spiec.easi(list(physeq_16s, physeq_its), 
                                 method='mb', nlambda=40, 
                                 lambda.min.ratio=1e-2, 
                                 pulsar.params = list(thresh = 0.05))

physeq = SpiecEasi::symBeta(SpiecEasi::getOptBeta(physeq), mode = "ave")

physeq = as.matrix(physeq)


taxnames = c(tax_table(physeq_16s)[,"Genus"], tax_table(physeq_its)[,"Genus"])
colnames(physeq) <- rownames(physeq) <- taxnames
diag(physeq) = 1

net_genus = netConstruct(data = physeq, 
                                 dataType = "condDependence", 
                                 sparsMethod = "none")

graph = igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)


corr = net_genus$edgelist1

colnames(corr) = c("from_name", "to_name", "asso", "diss", "adja")

graph = as_tbl_graph(graph)

graph = graph %>%
    activate(nodes) %>%
    mutate(Genus = name) %>%
    left_join(tax, by="Genus") %>%
    activate(edges) %>%
    mutate(from_name = .N()$name[from], to_name = .N()$name[to]) %>%
    left_join(corr, by=c("from_name", "to_name")) %>%
    filter(asso != "NA") %>%
    mutate(pair_asso_pos = if_else(asso>0, paste(.N()$Class[from], .N()$Class[to], sep="_"), NA), 
        pair_asso_neg = if_else(asso<0, paste(.N()$Class[from], .N()$Class[to], sep="_"), NA), 
        Association = if_else(asso>0, "positive", "negative")) %>%
    activate(nodes) %>%
    mutate(amplicon = if_else(name %in% c(tax_table(physeq_16s)[,"Genus"]), "Bacteria", "Fungi"),
        Hub_score = centrality_hub()) %>%
        activate(edges) %>%
        mutate(pair_asso_pos_amplicon = if_else(asso>0, paste(.N()$amplicon[from], .N()$amplicon[to], sep="_"), NA),
            pair_asso_neg_amplicon = if_else(asso<0, paste(.N()$amplicon[from], .N()$amplicon[to], sep="_"), NA))


colors_16s = c("#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fe4365", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#ecd078", "#d95b43")

colors_its = c("#018BFF", "#01FDFF", "#01FF61", "#F9FF01", "#FFD101", "#C08227", 
"#C02727", "#65216B", "#C177C7", "#EBBAE3")



graph_edges = graph %>%
    activate(edges) %>%
    data.frame()


graph_nodes = graph %>%
    activate(nodes) %>%
    data.frame()


nodes = matrix(data = 1:length(unique(graph_nodes$Class))*length(unique(graph_nodes$Class)),
       nrow = length(unique(graph_nodes$Class)),
       ncol = length(unique(graph_nodes$Class)))
nodes = as.data.frame(nodes)
rownames(nodes) = unique(graph_nodes$Class)
colnames(nodes) = unique(graph_nodes$Class)

nodes_class_pos = sapply(c(1:nrow(nodes)), function(x) sapply(c(1:ncol(nodes)), function(y) 
    length(grep(paste(rownames(nodes)[x], colnames(nodes)[y], sep="_"), graph_edges$pair_asso_pos)) + 
    length(grep(paste(colnames(nodes)[y], rownames(nodes)[x], sep="_"), graph_edges$pair_asso_pos))
))
nodes_class_pos = as.data.frame(nodes_class_pos)
rownames(nodes_class_pos) = rownames(nodes)
colnames(nodes_class_pos) = colnames(nodes)

nodes_class_pos$Class = rownames(nodes_class_pos)
nodes_class_pos = as.tibble(nodes_class_pos)
nodes_class_pos = pivot_longer(nodes_class_pos, cols=!Class, names_to="Class2", values_to="Positive")

nodes_class_neg = sapply(c(1:nrow(nodes)), function(x) sapply(c(1:ncol(nodes)), function(y) 
    length(grep(paste(rownames(nodes)[x], colnames(nodes)[y], sep="_"), graph_edges$pair_asso_neg)) + 
    length(grep(paste(colnames(nodes)[y], rownames(nodes)[x], sep="_"), graph_edges$pair_asso_neg))
))

nodes_class_neg = as.data.frame(nodes_class_neg)
rownames(nodes_class_neg) = rownames(nodes)
colnames(nodes_class_neg) = colnames(nodes)

nodes_class_neg$Class = rownames(nodes_class_neg)
nodes_class_neg = as.tibble(nodes_class_neg)

nodes_class_neg = pivot_longer(nodes_class_neg, cols=!Class, names_to="Class2", values_to="Negative")

####

associations_class = full_join(nodes_class_pos, nodes_class_neg, by=c("Class", "Class2"))
associations_class$Class = factor(associations_class$Class, colnames(nodes))
associations_class$Class2 = factor(associations_class$Class2, rev(colnames(nodes)))



nodes = matrix(data = 1:length(unique(graph_nodes$amplicon))*length(unique(graph_nodes$amplicon)),
       nrow = length(unique(graph_nodes$amplicon)),
       ncol = length(unique(graph_nodes$amplicon)))
nodes = as.data.frame(nodes)
rownames(nodes) = unique(graph_nodes$amplicon)
colnames(nodes) = unique(graph_nodes$amplicon)


nodes_amplicon_pos = sapply(c(1:nrow(nodes)), function(x) sapply(c(1:ncol(nodes)), function(y) 
    length(grep(paste(rownames(nodes)[x], colnames(nodes)[y], sep="_"), graph_edges$pair_asso_pos_amplicon)) + 
    length(grep(paste(colnames(nodes)[y], rownames(nodes)[x], sep="_"), graph_edges$pair_asso_pos_amplicon))
))
nodes_amplicon_pos = as.data.frame(nodes_amplicon_pos)
rownames(nodes_amplicon_pos) = rownames(nodes)
colnames(nodes_amplicon_pos) = colnames(nodes)

nodes_amplicon_pos$amplicon = rownames(nodes_amplicon_pos)
nodes_amplicon_pos = as.tibble(nodes_amplicon_pos)
nodes_amplicon_pos = pivot_longer(nodes_amplicon_pos, cols=!amplicon, names_to="amplicon2", values_to="Positive")

nodes_amplicon_neg = sapply(c(1:nrow(nodes)), function(x) sapply(c(1:ncol(nodes)), function(y) 
    length(grep(paste(rownames(nodes)[x], colnames(nodes)[y], sep="_"), graph_edges$pair_asso_neg_amplicon)) + 
    length(grep(paste(colnames(nodes)[y], rownames(nodes)[x], sep="_"), graph_edges$pair_asso_neg_amplicon))
))

nodes_amplicon_neg = as.data.frame(nodes_amplicon_neg)
rownames(nodes_amplicon_neg) = rownames(nodes)
colnames(nodes_amplicon_neg) = colnames(nodes)

nodes_amplicon_neg$amplicon = rownames(nodes_amplicon_neg)
nodes_amplicon_neg = as.tibble(nodes_amplicon_neg)

nodes_amplicon_neg = pivot_longer(nodes_amplicon_neg, cols=!amplicon, names_to="amplicon2", values_to="Negative")

associations_amplicon = full_join(nodes_amplicon_pos, nodes_amplicon_neg, by=c("amplicon", "amplicon2"))


associations_class = mutate(associations_class, Kingdom=if_else(Class2 %in% tax_table(physeq_16s)[,"Class"], "Bacteria", "Fungi"))


ggplot(associations_class) +
	geom_col(aes(x = Class, y = Positive, fill = Class2), size = 0.2) + 
    geom_col(aes(x = Class, y = -Negative, fill = Class2), size = 0.2) +
    geom_hline(aes(yintercept = 0), color = "black") +
	guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    scale_fill_manual(values = rev(c(colors_16s,colors_its))) +
    	labs(x = "Class", y = 'Associations ("+" - positive, "-" - negative)', fill = "Class") +
	scale_y_continuous(limits = c(-30,70), n.breaks=10) +
    geom_text(aes(x = Class, y = Positive, label = ifelse(Positive>2,Positive, "")), inherit.aes = FALSE, size = 4, position = position_stack(vjust = 0.5)) +
    geom_text(aes(x = Class, y = -Negative, label = ifelse(Negative<c(-2),Negative, "")), inherit.aes = FALSE, size = 4, position = position_stack(vjust = 0.5)) +
theme_bw() +
	theme(
	legend.text = element_text(size=15),
	axis.text.x = element_text(size=15, angle=-90, vjust = 1, face = "italic", color="black"),
	axis.title.x = element_text(size=10),
	)

ggsave("Figure 6c.png", width = 10, height = 8, dpi=400, limitsize=FALSE)

associations_amplicon$amplicon2 = factor(associations_amplicon$amplicon2, c("Fungi", "Bacteria"))

ggplot(associations_amplicon) +
	geom_col(aes(x = amplicon, y = Positive, fill = amplicon2), colour = "black", linetype = "solid", size = 0.1) + 
    geom_col(aes(x = amplicon, y = -Negative, fill = amplicon2), colour = "black", linetype = "solid", size = 0.1) +
	guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    geom_hline(aes(yintercept = 0), color = "black") +
    scale_fill_manual(values = rev(c("#4285f4", "#f25022"))) +
    	labs(x = "Kingdom", y = 'Associations ("+" - positive, "-" - negative)', fill = "Kingdom") +
	scale_y_continuous(limits = c(-60,180), n.breaks=12) +
    geom_text(aes(x = amplicon, y = Positive, label = Positive), inherit.aes = FALSE, size = 4, position = position_stack(vjust = 0.5)) +
    geom_text(aes(x = amplicon, y = -Negative, label = Negative), inherit.aes = FALSE, size = 4, position = position_stack(vjust = 0.5)) +
theme_bw() +
	theme(
legend.text = element_text(size=15),
	axis.text.x = element_text(size=15, angle=0, vjust = 1, face = "italic", color="black"),
	axis.title.x = element_text(size=10),
	)

ggsave("Figure 6b.png", width = 6, height = 8, dpi=300, limitsize=FALSE)


graph = graph %>%
    activate(nodes) %>%
    arrange(amplicon, Class, desc(Hub_score)) %>%
    mutate(Class = factor(Class,unique(associations_class$Class)), name = factor(name,name),
    id = row_number(),
    angle = 90 - 360 * id / length(name),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle+180, angle))


ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(x=x, xend=xend, y=y, yend=yend, edge_colour = Association)) +
  geom_node_point(aes(x=x, y=y, size = Hub_score, color=Class, shape=Kingdom), alpha=0.8) +
  geom_node_text(aes(x=x*1.05, y=y*1.05, label = name, angle = angle, hjust=hjust, size=0.04, fontface=if_else(Hub_score>as.numeric(quantile(Hub_score, probs=0.90)),"bold","plain"))) + 
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_color_manual(values = c(colors_16s,colors_its)) +
  theme_classic() + th_no_axes() + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"),
    legend.box.margin=unit(c(0,0,0,0),"cm"),
    aspect.ratio=1,
    legend.position="bottom", legend.box = "vertical") 

ggsave("Figure 6a.png", width = 13, height = 13, dpi=500, limitsize=FALSE)

