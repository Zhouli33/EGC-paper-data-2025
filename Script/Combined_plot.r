library(ggplot2)
library(dplyr)
library(patchwork)


df <- read.csv("groupby_group_freq.csv")
df$new_celltype <- factor(df$new_celltype, levels = rev(c("Epithelial","Plasma_cells","B_cells","T_NK_cells","Mast_cell","Myeloid","Fibroblast","Endothelial")))
df$group <- factor(df$group, levels = c("NAG","CAG_IM","TNM_1","TNM_2","TNM_3","TNM_4"))

colours=c("Epithelial" = "#7fc97f" ,"T_NK_cells"="#beaed4","Plasma_cells"="#fdc086","Fibroblast"="#386cb0","B_cells"="#f0027f","Mast_cell"="#a34e3b","Myeloid"="#666666","Endothelial"="#1b9e77")

#dot plot
dot_plot <- ggplot(df, aes(x = group, y = new_celltype, size = freq, color = new_celltype)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  labs(x = "", y = "Cell proportion") +
  theme(axis.text.x = element_text(color = "black",angle = 90, vjust = 0.5, hjust=1),
  		  axis.text.y = element_text(color = "black"),
  		  panel.grid = element_blank(), 
		    axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
		    panel.border = element_rect(color = "black", fill = NA, size = 1),
		    legend.key = element_blank(),
		    legend.position = "left",
		    plot.margin = margin()  # 去除绘图margin
	    )+
  scale_color_manual(values = colours)+
  guides(color = "none", size = guide_legend(title = "Cell proportion")) 


# Stacked bar plot
data <- df[, c("group", "new_celltype", "cell_number", "total")]
replace_map <- c("B_cells" = "Lymphoid",
                 "T_NK_cells" = "Lymphoid",
                 "Mast_cell" = "Myeloid",
                 "Myeloid" = "Myeloid",
                 "Fibroblast" = "Stromal",
                 "Endothelial" = "Stromal",
                 "Epithelial" = "Epithelial",
                 "Plasma_cells" = "Plasma_cells")
data$new_celltype <- as.character(data$new_celltype)
data$new_celltype <- replace(data$new_celltype, data$new_celltype %in% names(replace_map), replace_map[data$new_celltype])

merged_data <- data %>%
  group_by(group, new_celltype) %>%
  summarize(
    cell_number = sum(cell_number),
    total = first(total), 
  ) %>%
  ungroup() %>%
  mutate(freq = cell_number / total) %>%
  ungroup()
  merged_data <- as.data.frame(merged_data)
merged_data$new_celltype <- factor(merged_data$new_celltype,levels=c("Epithelial","Plasma_cells","Lymphoid","Myeloid","Stromal"))
palettes <- c("#7fc97f","#fdc086","#FE9D97","#736A2B","#5B5BE6")
write.csv(merged_data,"merged_data.csv",quote=F,row.names=F)

stacked_bar_plot <- ggplot(merged_data, aes(x = group, y = cell_number, fill = new_celltype)) +
  geom_bar(stat = "identity",width=0.5) +
  scale_fill_manual(values = palettes) +
  labs( y = "Number of cells")+
  theme(panel.grid = element_blank(), 
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_blank(),
				axis.title.x = element_blank(),
                axis.text.y = element_text(colour = "black"),
                panel.border = element_blank(),
                axis.line = element_line(colour = "black"),
				axis.ticks.length.x = unit(0, "mm"),# 无缝拼接 ticks去除的关键，ticks的绘图区域调为0
				plot.margin = margin(),
				legend.position = "left",
				legend.key.size= unit(0.4, "cm"))+
  scale_y_continuous(expand=c(0,0))

# Horizontal bar plot
ordered_levels <- c("Epithelial","Plasma_cells","B_cells","T_NK_cells","Mast_cell","Myeloid","Fibroblast","Endothelial")
horizontal_bar_plot <- ggplot(df, aes(x = cell_number, y = new_celltype, fill = new_celltype)) +
  geom_bar(stat = "identity",width=0.5) +
  scale_fill_manual(values = colours,breaks = ordered_levels) +
  labs(x = "Number of cells", y = "",fill = "Metacluster") +
  theme(
		panel.grid = element_blank(), 
    	panel.background = element_rect(fill = "white"),
		axis.text.y = element_blank(),
		axis.line = element_line(colour = "black"),
		axis.ticks.y = element_blank(),
		plot.margin = margin(),
		axis.ticks.length.y = unit(0, "mm")
  )+
  scale_x_continuous(expand=c(0,0))

#拼图
combined_plot <- (stacked_bar_plot/dot_plot)+plot_layout(heights = c(0.4, 1))|(plot_spacer()/horizontal_bar_plot)+plot_layout(heights = c(0.4, 1)) #plot_layout 调整图形长宽比例

ggsave("combined_plot.pdf",combined_plot,height=6,width=10,bg="white")
ggsave("combined_plot.png",combined_plot,height=6,width=10,bg="white")
