#Loading in the libraries, obviously
library(ggplot2)
library(dplyr)
library(reshape2)

#Loading the interval read depth data, selecting just the read depths, and then removing the full data to save RAM
full_data <- read.table("hg38_6mil_200bp_ID_UMAP50_golden_iter0-scored-unscored_iter1-scored-unscored_map1-scored-unscored.tab")
map_reads <- subset(full_data, select = 5:12)
colnames(map_reads) = c("UMAP50", "gold", "i0s", "i0u", "i1s", "i1u", "m1s", "m1u")
rm(full_data)

#Creating a range by which to bin UMAP values, creating bin factor with cut, and computing mean of each column for that bin
map_bins <- seq(0, 1.01, by=0.01)
map_reads <- map_reads %>%
  mutate(UMAPbin = cut(UMAP50, map_bins, include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01)))
reads_by_umap <- map_reads %>%
  group_by(UMAPbin) %>%
  summarise_at(vars(UMAP50:m1u), mean)

#Creating a melted dataframe for easy plotting
plot_rbu <- melt(reads_by_umap, id.vars = "UMAPbin", measure.vars = 3:9)
plot_rbu$UMAP <- as.numeric(levels(plot_rbu$UMAPbin))[plot_rbu$UMAPbin]

ggplot(data = plot_rbu, aes(x = UMAP, y = value, color = variable, group = variable)) + theme_classic() + 
  geom_line() + scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 50), expand = c(0,0)) + labs(xaxis = "Average UMAP Score", yaxis = "Average Read Depth")

#Select only rows where gold>0
map_gold <- map_reads %>%
  filter(gold>0)

#Computing mean of each column for each UMAP bin only for rows where gold>0
gold_reads_by_umap <- map_gold %>%
  group_by(UMAPbin) %>%
  summarise_at(vars(UMAP50:m1u), mean)

#Creating melted dataframe for easy plotting
plot_gold_rbu <- melt(gold_reads_by_umap, id.vars = "UMAPbin", measure.vars = 3:9)
plot_gold_rbu$UMAP <- as.numeric(levels(plot_gold_rbu$UMAPbin))[plot_gold_rbu$UMAPbin]

ggplot(data = plot_gold_rbu, aes(x = UMAP, y = value, color = variable, group = variable)) + theme_classic() + 
  geom_line() + scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 50), expand = c(0,0)) + labs(xaxis = "Average UMAP Score", yaxis = "Average Read Depth")

#Sorting each column of map_gold to generate a percentile plot (from gold_sub, with only ~6000 points so my GPU survives)
gold_sort <- map_gold %>% select(gold:m1u) %>% mutate_each(list(~sort))
gold_sub <- gold_sort %>%
  select(gold:m1u) %>%
  mutate(ptile = row_number()/length(gold)*100) %>%
  filter(row_number() %% 1000 == length(gold) %% 1000)

