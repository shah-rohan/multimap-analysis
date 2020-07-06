library(ggplot2)
library(dplyr)
library(reshape2)

full_data <- read.table("hg38_6mil_200bp_ID_UMAP50_golden_iter0-scored-unscored_iter1-scored-unscored_map1-scored-unscored.tab")
map_reads <- subset(full_data, select = 5:12)
rm(full_data)
colnames(map_reads) = c("UMAP50", "gold", "i0s", "i0u", "i1s", "i1u", "m1s", "m1u")

map_bins <- seq(0, 1.01, by=0.01)
map_reads <- map_reads %>%
  mutate(UMAPbin = cut(UMAP50, map_bins, include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01)))
reads_by_umap <- map_reads %>%
  group_by(UMAPbin) %>%
  summarise_at(vars(UMAP50:m1u), mean)

plot_rbu <- melt(reads_by_umap, id.vars = "UMAPbin", measure.vars = 3:9)
plot_rbu$UMAP <- as.numeric(levels(plot_rbu$UMAPbin))[plot_rbu$UMAPbin]
ggplot(data = plot_rbu, aes(x = UMAP, y = value, color = variable, group = variable)) + theme_classic() + 
  geom_line() + scale_x_continuous(breaks = seq(0,1,0.1))

gold_reads_by_umap <- map_reads %>%
  filter(gold>0) %>%
  group_by(UMAPbin) %>%
  summarise_at(vars(UMAP50:m1u), mean)

plot_gold_rbu <- melt(gold_reads_by_umap, id.vars = "UMAPbin", measure.vars = 3:9)
plot_gold_rbu$UMAP <- as.numeric(levels(plot_gold_rbu$UMAPbin))[plot_gold_rbu$UMAPbin]
ggplot(data = plot_gold_rbu, aes(x = UMAP, y = value, color = variable, group = variable)) + theme_classic() + 
  geom_line() + scale_x_continuous(breaks = seq(0,1,0.1))
