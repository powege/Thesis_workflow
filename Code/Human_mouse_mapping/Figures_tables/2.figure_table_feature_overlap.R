rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)

### SET VARS
m_matrix_file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/m.mus_GRC38_v101_whole_genome_features_multicell_overlap.matrix"
h_matrix_file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/h.sap_GRC38_v101_whole_genome_features_multicell_overlap.matrix"
m_plot_file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_m.mus_GRC38_v101_whole_genome_features_multicell_overlap.jpg"
h_plot_file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_h.sap_GRC38_v101_whole_genome_features_multicell_overlap.jpg"


### IMPORT 
m_matrix <- fread(m_matrix_file)
h_matrix <- fread(h_matrix_file)

### FORMAT
m_matrix <- as.matrix(m_matrix)
row.names(m_matrix) <- colnames(m_matrix)
m_matrix <- setNames(melt(m_matrix), c('annotation_A', 'annotation_B', '%_A_overlap_B'))
m_matrix$`%_A_overlap_B` <- round(m_matrix$`%_A_overlap_B`, 1)

h_matrix <- as.matrix(h_matrix)
row.names(h_matrix) <- colnames(h_matrix)
h_matrix <- setNames(melt(h_matrix), c('annotation_A', 'annotation_B', '%_A_overlap_B'))
h_matrix$`%_A_overlap_B` <- round(h_matrix$`%_A_overlap_B`, 1)

### PLOT 
m_ggheatmap <- ggplot(m_matrix, aes(annotation_B, annotation_A, fill = `%_A_overlap_B`)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#D55E00", high = "#56B4E9", mid = "white", 
                       midpoint = 50, limit = c(0,100), space = "Lab", 
                       name="A overlap B (%)") +
  xlab("B") +
  ylab("A") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(annotation_B, annotation_A, label = `%_A_overlap_B`, fontface=2), color = "black", size = 5) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
    text = element_text(size = 14)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
m_ggheatmap 

h_ggheatmap <- ggplot(h_matrix, aes(annotation_B, annotation_A, fill = `%_A_overlap_B`)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#D55E00", high = "#56B4E9", mid = "white", 
                       midpoint = 50, limit = c(0,100), space = "Lab", 
                       name="A overlap B (%)") +
  xlab("B") +
  ylab("A") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(annotation_B, annotation_A, label = `%_A_overlap_B`, fontface=2), color = "black", size = 5) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
    text = element_text(size = 14)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
h_ggheatmap 

### EXPORT
ggsave(filename = m_plot_file, plot = m_ggheatmap, height = 7, width = 9)
ggsave(filename = h_plot_file, plot = h_ggheatmap, height = 7, width = 9)
  






