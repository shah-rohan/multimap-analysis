library(data.table)
library(tidyverse)

k51_data = fread("data/hg38_6mil_200bp_ID_UMAP50_golden_iter0-iter1-map1_scored-unscored.tab")
colnames(k51_data) = c("chr", "start", "stop", "ID", "UMAP","gold", "iter0", "iter0u", "k51MM", "i1u", "uni", "uniu")
data = k51_data %>% select(ID, UMAP, gold, k51MM, uni)
rm(k51_data)

k101_data = fread("data/hg38_6mil_200bp_ID_UMAP50_golden_k101.tab")
colnames(k101_data) = c("chr", "start", "stop", "ID", "UMAP", "gold", "k101MM")
data = data %>% bind_cols(., k101_data$k101MM)
colnames(data) = c("ID", "UMAP", "gold", "k51mm", "uni", "k101mm")
rm(k101_data)

map_bins <- seq(0, 1.01, by=0.01)
map_reads <- data %>%
	mutate(UMAPbin = cut(UMAP, map_bins, include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01)))
map_reads <- map_reads %>%
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin])

#Only selecting for those regions that have a golden read depth > 0
map_gold <- map_reads %>%
	filter(gold > 0)

#Count the number of elements in each of the UMAPbins of map_gold
UMAP_bins = map_gold %>% 
	group_by(UMAPbin) %>% 
	count()
ggplot(data = UMAP_bins, aes(x = UMAPbin, y = n)) + geom_line() +
	theme_classic() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(limits = c(1e3, 1e7), expand = c(0,0)) +
	labs(x = "UMAP Score", y = "Number of Regions")
rm(UMAP_bins)

#Computing the mean of each read depth column binned by UMAP, then plotting
reads_by_umap <- map_gold %>%
	group_by(UMAPbin) %>%
	summarise_at(vars(gold:k101mm), mean)
plot_reads_by_umap <- reads_by_umap %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_reads_by_umap, aes(x = UMAPbin, y = value, color = variable)) + geom_line() +
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 50), expand = c(0,0)) +
	labs(x = "Average UMAP Score", y = "Average Read Depth")

#Keep my workspace clean!
rm(reads_by_umap, plot_reads_by_umap)

#Sorting each column of map_gold to generate a percentile plot (from gold_sub, with only ~6000 points so my GPU survives)
gold_sort <- map_gold %>%
	select(gold:k101mm) %>%
	mutate_each(list(~sort(.)))
gold_sub <- gold_sort %>%
	mutate(ptile = row_number()/length(gold)*100) %>%
	filter(row_number() %% 1000 == 0)
plot_gold_sub <- gold_sub %>%
	gather(variable, value, -ptile)
ggplot(data = plot_gold_sub, aes(x = ptile, y = value, color = variable)) + geom_line() + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 80), expand = c(0,0)) + 
	labs(x = "Percentile", y = "Read Depth")

mae_by_UMAP = map_gold %>%
	group_by(UMAPbin) %>%
	transmute_at(vars(k51mm:k101mm), ~abs(.-gold)) %>%
	summarize_all(mean)
plot_mae_umap <- mae_by_UMAP %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_mae_umap, aes(x = UMAPbin, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	labs(xaxis = "Average UMAP Score", yaxis = "Excess Read Depth")
rm(mae_by_UMAP, plot_mae_umap)