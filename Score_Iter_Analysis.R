#Loading in the libraries, obviously
library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

#Loading the interval read depth data, selecting just the read depths, and then removing the full data to save RAM
full_data <- read.table("data/hg38_6mil_200bp_ID_UMAP50_golden_iter0-iter1-map1_scored-unscored.tab")
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
	summarise_at(vars(i0s:m1u), mean)
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

#Keep my workspace clean!
rm(gold_sub, plot_gold_sub)

#Analyze excess read depth over map-1
map_excess <- map_gold %>%
	transmute(UMAPbin, i0sdiff = i0s - m1s, i1sdiff = i1s - m1s,
		i0udiff = i0u - m1u, i1udiff = i1u - m1u)
map_excess_by_umap <- map_excess %>%
	group_by(UMAPbin) %>%
	summarise_all(mean)
plot_excess_rbu <- map_excess_by_umap %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_excess_rbu, aes(x = UMAPbin, y = value, color = variable)) + geom_line() + 
	theme_classic() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	labs(x = "Average UMAP Score", y = "Excess Read Depth")

#Keep my workspace clean!
rm(map_excess_by_umap, plot_excess_rbu)

#Sorting each column to generate a percentile plot of excesses
map_excess_sorted <- map_excess %>% 
	select(i0sdiff:i1udiff) %>% 
	mutate_each(list(~sort))
excess_sub <- map_excess_sorted %>%
	mutate(ptile = row_number()/length(i0sdiff)*100) %>%
	filter(row_number() %% 1000 == 0)
plot_excess_ptile <- excess_sub %>%
	gather(variable, value, -ptile)
ggplot(data = plot_excess_ptile, aes(x = ptile, y = value, color = variable)) + geom_line() +
	theme_classic() + 
	scale_x_continuous(limits = c(0, 100), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 60), expand = c(0,0)) +
	labs(x = "Percentile", y = "Excess Read Depth")

#Keep my workspace clean!
rm(excess_sub, plot_excess_ptile)

#Define the color palette that I'll be using so I don't have to keep dealing with it.
scl_pal <- wes_palette("Zissou1", type="continuous")

#Generate several qq plots of different samples
plt_gold = gold_sort %>%
	filter(row_number() %% ceiling(length(gold)/100) == 0) %>%
	mutate(ptile = 1:99)

#QQ plot of Map 1 scored vs Iteration 1 Scored
ggplot(data = plt_gold, aes(x = m1s, y = i1s, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "Map 1 Scored", y = "Iteration 1 Scored") 
ggsave(filename = "QQ-M1S-I1S.eps", width = 3, height = 3, units = "in")

#QQ Plot of Iteration 1 scored vs golden
ggplot(data = plt_gold, aes(x = i1s, y = gold, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "Iteration 1 Scored", y = "Golden") 
ggsave(filename = "QQ-I1S-Golden.eps", width = 3, height = 3, units = "in")

#QQ Plot of Map 1 scored vs golden
ggplot(data = plt_gold, aes(x = m1s, y = gold, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "Map 1 Scored", y = "Golden") 
ggsave(filename = "QQ-M1S-Golden.eps", width = 3, height = 3, units = "in")

#Keep my workspace clean!
rm(plt_gold)

#Generate plots of mean absolute error by UMAP
mae_by_UMAP = map_gold %>%
	group_by(UMAPbin) %>%
	transmute_at(vars(i0s:m1u), ~abs(.-gold)) %>%
	summarize_all(mean)
plot_mae_umap <- mae_by_UMAP %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_mae_umap, aes(x = UMAPbin, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	labs(xaxis = "Average UMAP Score", yaxis = "Excess Read Depth")
rm(mae_by_UMAP, plot_mae_umap)

dens = density(map_gold$UMAP50, adjust = 2, from=0, to = 1, n = 200, bw = "nrd")
dens2 = as.data.frame(dens$x) %>% bind_cols(as.data.frame(dens$y))
colnames(dens2) = c("x", "y")
ggplot(dens2, aes(x = x, y = y)) + geom_line() + scale_y_log10()

#Compute Mean absolute error stratified by golden read depth
print(map_reads %>% select(gold:m1u) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
print(map_reads %>% select(gold:m1u) %>% filter(gold < 30) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
print(map_reads %>% select(gold:m1u) %>% filter(gold < 40 & gold >=30) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
print(map_reads %>% select(gold:m1u) %>% filter(gold < 50 & gold >=40) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
print(map_reads %>% select(gold:m1u) %>% filter(gold < 60 & gold >=50) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
print(map_reads %>% select(gold:m1u) %>% filter(gold >= 60) %>% mutate(across(-1, ~abs((.-gold)/mean(gold)))) %>% summarize_all(mean))
