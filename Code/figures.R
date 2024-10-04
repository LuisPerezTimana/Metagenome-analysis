
# IMPORT LIBRARIES
library(biomformat)
library(tidyverse)
library(microeco)
library(magrittr)
library(ape)
library(seqinr)
library(gridExtra)

# IMPORT DATASET
mt <- readRDS("../Data/metagenome.Rds")

# Calculate metrics
mt$cal_abund()
mt$cal_alphadiv()
mt$cal_betadiv()

setwd("../Results")

##### FIGURES
## RELATIVE ABUNDANCE
# ABUNDANCE BY SPECIE
t1 <- trans_abund$new(dataset = mt, taxrank = "Species", ntaxa = 20)
g1 <- t1$plot_bar(bar_full = FALSE, use_alluvium = TRUE, clustering = T, 
                  xtext_angle = 90, , legend_text_italic = T,
                  color_values = RColorBrewer::brewer.pal(8, "Set2"), facet = "Group")

# HEATMAP OF ABUNDANCE BY GENUS
t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 40)
g3 <- t1$plot_heatmap(facet = "Group", withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10), xtext_angle = 90)
g3 <- g3 + theme(axis.text.y = element_text(face = 'italic'))

# SAVE FIGURES
ggsave(filename = "Figure 1a.tiff", plot = g1, dpi = 300, height = 10, width = 20)
ggsave(filename = "Figure 1c.tiff", plot = g3, dpi = 300,  height = 10, width = 15)


## ALPHA INDEX DIVERSITY
# BOXPLOTS OF THE INDEX FILTERED BY SIGNIFICANCE
t1 <- trans_alpha$new(dataset = mt, group = "Group")
t1$cal_diff(method = "wilcox")
t1$res_diff %<>% subset((Significance %in% "***") | (Significance %in% "**") | (Significance %in% "*"))

# FOR EACH INDEX
index <- unique(t1$res_diff$Measure)
list_plots = list()
for (i in 1:length(index)){
  list_plots[[i]] <- t1$plot_alpha(measure = index[i], shape = "Group", add = "jitter", add_sig_label = "P.adj")
}

# JOIN THE PLOTS
fig2 <- gridExtra::grid.arrange(list_plots[[1]], list_plots[[2]], list_plots[[3]],
                        list_plots[[4]], list_plots[[5]], list_plots[[6]],
                        list_plots[[7]], list_plots[[8]], list_plots[[9]], ncol = 3, nrow = 3)

# SAVE FIGURE
ggsave(filename = "Figure 2.tiff", plot = fig2, dpi = 300, height = 15, width = 20)


## BETA INDICES
## BRAY INDEX
t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "bray")
# BOXPLOT OF INTRAGROUP DISTANCES
t1$cal_group_distance(within_group = T)
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff %<>% subset((Significance %in% "***") | (Significance %in% "**") | (Significance %in% "*"))
g1 <- t1$plot_group_distance(add = "mean", add_sig_label = "P.adj")

# BOXPLOT OF INTERGROUP DISTANCES
t1$cal_group_distance(within_group = F)
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff %<>% subset((Significance %in% "***") | (Significance %in% "**") | (Significance %in% "*"))
g2 <- t1$plot_group_distance(add = "mean", add_sig_label = "P.adj")

## JACCARD INDEX
t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "jaccard")
# BOXPLOT OF INTRAGROUP DISTANCES
t1$cal_group_distance(within_group = T)
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff %<>% subset((Significance %in% "***") | (Significance %in% "**") | (Significance %in% "*"))
g3 <- t1$plot_group_distance(add = "mean", add_sig_label = "P.adj")

# BOXPLOT OF INTERGROUP DISTANCES
t1$cal_group_distance(within_group = F)
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff %<>% subset((Significance %in% "***") | (Significance %in% "**") | (Significance %in% "*"))
g4 <- t1$plot_group_distance(add = "mean", add_sig_label = "P.adj")

# JOIN THE PLOTS
fig3 <- gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)

# SAVE THE FIGURE
ggsave(filename = "Figure 3.tiff", plot = fig3, dpi = 300, height = 15, width = 20)


### PCOA AND CLUSTERING

# BRAY INDEX PCO
t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "bray")
t1$cal_ordination(method = "PCoA")
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))

# BRAY INDEX CLUSTERING
g2 <- t1$plot_clustering(group = "Group", replace_name = c("Group"))

# JACCARD INDEX PCO
t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "jaccard")
t1$cal_ordination(method = "PCoA")
g3 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))

# JACCARD INDEX CLUSTERING
g4 <- t1$plot_clustering(group = "Group", replace_name = c("Group"))

# JOIN THE PLOTS
fig4 <- gridExtra::grid.arrange(g1, g3, g2, g4, ncol = 2, nrow = 2)
# SAVE THE IMAGE
ggsave(filename = "Figure 4.tiff", plot = fig4, dpi = 300, height = 15, width = 15)

### LDA AND RANDOM FOREST ANALYSIS
## LDA FOR SPECIES
# JOIN THE SIMILAR GROUPS (G1-G2 AND G3-G4)
mt$sample_table$Group2 <- c(rep("G1-G2", 20), rep("G3-G4", 12))
t1 <- trans_diff$new(dataset = mt, method = "lefse", group = "Group2", alpha = 0.05, taxa_level = "Species")

# FILTER TAXA
lista <- c()
for (i in t1$res_diff$Taxa){
  specie <- unlist(strsplit(i, split = "s__"))
  lista <- c(lista, specie[2])

}

t1$res_diff$Taxa <- lista
g1 <- t1$plot_diff_bar(threshold = 0, width = 0.8) + theme(axis.text.y = element_text(face = "italic"))

## MEAN DECREASING GINI BY RANDOM FOREST
t1 <- trans_diff$new(dataset = mt, method = "rf", group = "Group2", taxa_level = "Species")

# FILTER TAXA
lista <- c()
for (i in t1$res_diff$Taxa){
  specie <- unlist(strsplit(i, split = "s__"))
  lista <- c(lista, specie[2])
  
}

t1$res_diff$Taxa <- lista
g2 <- t1$plot_diff_bar(use_number = 1:50, width = 0.8)  + theme(axis.text.y = element_text(face = "italic"))


## ABUNDANCE OF POTENTIAL BIOMARKERS
t1 <- trans_diff$new(dataset = mt, method = "KW_dunn", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
g3 <- t1$plot_diff_abund(use_number = 1:20, add_sig = F, coord_flip = F, keep_prefix = F)  + theme(axis.text.x = element_text(face = "italic"))

# JOIN THE PLOTS
fig5 <- gridExtra::grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 2),c(3, 3)))

# SAVE THE FIGURE
ggsave(filename = "Figure 5.tiff", plot = fig5, dpi = 300, height = 15, width = 15)

### FUNTIONAL ANALYSIS
t2 <- trans_func$new(mt)
t2$cal_spe_func(prok_database = "FAPROTAX")
t2$cal_spe_func_perc()

## FUNTIONAL HEATMAP
g1 <- t2$plot_spe_func_perc()

## BOXPLOT OF PRODUCTION AND CONSUME FUNCTIONS
t1 <- trans_func$new(mt)
# Select NJC19 database
t1$cal_spe_func(prok_database = "NJC19")
# get the trait percentage data
t1$cal_spe_func_perc()

mt2 <- mt$clone()
# inset the trait percentage result into taxa_abund of microtable object
mt2$taxa_abund$Trait <- as.data.frame(t(t1$res_spe_func_perc))
# use trans_abund to plot
t1 <- trans_abund$new(dataset = mt2, taxrank = "Trait", ntaxa = 10, use_percentage = FALSE)
g2 <- t1$plot_box(group = "Group", xtext_angle = 30) + ylab("Relative population abundance (%)") + theme(axis.text.x = element_text(size = 13))

# JOIN THE PLOTS
fig6 <- gridExtra::grid.arrange(g1, g2, layout_matrix = rbind(c(1, 1),
                                                              c(2, 2)))
# SAVE THE FIGURE
ggsave(filename = "Figure 6.tiff", plot = fig6, dpi = 300, height = 17, width = 17)


