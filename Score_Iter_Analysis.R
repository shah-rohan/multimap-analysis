#Loading in the libraries, obviously
library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

#Loading the interval read depth data, selecting just the read depths, and then removing the full data to save RAM
full_data <- read.table("hg38_6mil_200bp_ID_UMAP50_golden_iter0-scored-unscored_iter1-scored-unscored_map1-scored-unscored.tab")
map_reads <- subset(full_data, select = 5:12)
colnames(map_reads) = c("UMAP50", "gold", "i0s", "i0u", "i1s", "i1u", "m1s", "m1u")
rm(full_data)

#Creating a range by which to bin UMAP values and creating bins with cut
map_bins <- seq(0, 1.01, by=0.01)
map_reads <- map_reads %>%
	mutate(UMAPbin = cut(UMAP50, map_bins, include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01)))
map_reads <- map_reads %>%
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin])

#Only selecting for those regions that have a golden read depth > 0
map_gold <- map_reads %>%
	filter(gold > 0)

#Computing the mean of each read depth column binned by UMAP, then plotting
reads_by_umap <- map_gold %>%
	group_by(UMAPbin) %>%
	summarise_at(vars(i0s:m1u), mean)
plot_reads_by_umap <- reads_by_umap %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_reads_by_umap, aes(x = UMAPbin, y = value, color = variable)) + geom_line() +
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 50), expand = c(0,0)) +
	labs(x = "Average UMAP Score", y = "Average Read Depth")

#Sorting each column of map_gold to generate a percentile plot (from gold_sub, with only ~6000 points so my GPU survives)
gold_sort <- map_gold %>%
	select(gold:m1u) %>%
	mutate_each(list(~sort))
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

#Analyze excess read depth over map-1
map_excess <- map_gold %>%
	transmute(UMAP50, UMAPbin, i0sdiff = i0s - m1s, i1sdiff = i1s - m1s, i0udiff = i0u - m1u, i1udiff = i1u - m1u)

map_excess_by_umap <- map_excess %>%
	group_by(UMAPbin) %>%
	summarise_at(vars(2:5), mean)

plot_excess_rbu <- melt(map_excess_by_umap, id.vars = "UMAPbin")
plot_excess_rbu$UMAP <- as.numeric(levels(plot_excess_rbu$UMAPbin))[plot_excess_rbu$UMAPbin]
ggplot(data = plot_excess_rbu, aes(x = UMAP, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	labs(xaxis = "Average UMAP Score", yaxis = "Excess Read Depth")

map_excess_sorted <- map_excess %>% select(i0sdiff:i1udiff) %>% mutate_each(list(~sort))
excess_sub <- map_excess_sorted %>% mutate(ptile = row_number()/length(i0sdiff)*100) %>% filter(row_number() %% 1000 == 0)

plot_excess_ptile <- melt(excess_sub, id.vars = "ptile")
ggplot(data = plot_excess_ptile, aes(x = ptile, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 60), expand = c(0,0))

UMAP_bins = map_gold %>% group_by(UMAPbin) %>% count()
UMAP_bins$UMAPbin = as.numeric(levels(UMAP_bins$UMAPbin))[UMAP_bins$UMAPbin]
ggplot(data = UMAP_bins, aes(x = UMAPbin, y = n)) + geom_line() + theme_classic() + scale_y_log10()

plt_gold = gold_sort %>% filter(row_number() %% 28540 == 0) %>% mutate(ptile = seq(0.5, 99.5, 0.5))
ggplot(data = plt_gold, aes(x = m1s, y = i0s, color = ptile)) + geom_point() + 
	scale_color_gradientn(colors = wes_palette("Zissou1", 100, type="continuous"), name = "Percentile", limits = c(0, 100), position = "left") + 
	theme_classic() + scale_x_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	labs(x = "Map 1 Scored", y = "Iteration 1 Scored") + geom_abline(slope = 1)
ggsave(filename = "QQ-M1S-I1S.eps", width = 3, height = 3, units = "in")

ggplot(data = plt_gold, aes(x = i1s, y = gold, color = ptile)) + geom_point() + 
	scale_color_gradientn(colors = wes_palette("Zissou1", 100, type="continuous"), name = "Percentile", limits = c(0, 100), position = "left") + 
	theme_classic() + scale_x_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	labs(x = "Iteration 1 Scored", y = "Golden") + geom_abline(slope = 1)
ggsave(filename = "QQ-I1S-Golden.eps", width = 3, height = 3, units = "in")

ggplot(data = plt_gold, aes(x = m1s, y = gold, color = ptile)) + geom_point() + 
	scale_color_gradientn(colors = wes_palette("Zissou1", 100, type="continuous"), name = "Percentile", limits = c(0, 100), position = "left") + 
	theme_classic() + scale_x_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	labs(x = "Map 1 Scored", y = "Golden") + geom_abline(slope = 1)
ggsave(filename = "QQ-M1S-Golden.eps", width = 3, height = 3, units = "in")

ggplot(data = plt_gold, aes(x = i0s, y = gold, color = ptile)) + geom_point() + 
	scale_color_gradientn(colors = wes_palette("Zissou1", 100, type="continuous"), name = "Percentile", limits = c(0, 100), position = "left") + 
	theme_classic() + scale_x_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	labs(x = "Iteration 0 Scored", y = "Golden") + geom_abline(slope = 1)

ggplot(data = plt_gold, aes(x = i1u, y = i1s, color = ptile)) + geom_point() + 
	scale_color_gradientn(colors = wes_palette("Zissou1", 100, type="continuous"), name = "Percentile", limits = c(0, 100), position = "left") + 
	theme_classic() + scale_x_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	labs(x = "Iteration 0 Scored", y = "Iteration 1 Scored") + geom_abline(slope = 1)

mae_by_UMAP = map_gold %>% group_by(UMAPbin) %>% transmute_at(vars(i0s:m1u), ~abs(.-gold)) %>% summarize_all(mean)
plot_mae_umap <- melt(mae_by_UMAP, id.vars = "UMAPbin")
plot_mae_umap$UMAP <- seq(0, 1, 0.01)
ggplot(data = plot_mae_umap, aes(x = UMAP, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	labs(xaxis = "Average UMAP Score", yaxis = "Excess Read Depth")
