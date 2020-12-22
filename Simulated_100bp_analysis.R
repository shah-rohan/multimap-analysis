library(tidyverse)

fulldata100bp = read.table("data/hg38_6mil_200bp_ID_UMAP50_100bp-golden-multimap-uniread.tab")
data100bp = fulldata100bp %>% select(5:8)
colnames(data100bp) = c("UMAP50", "gold", "multimap", "uniread")
rm(fulldata100bp)

map_bins <- seq(0, 1.01, by=0.01)
map_reads <- data100bp %>%
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
	summarise_at(vars(gold:uniread), mean)
plot_reads_by_umap <- reads_by_umap %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_reads_by_umap, aes(x = UMAPbin, y = value, color = variable)) + geom_line() +
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 25), expand = c(0,0)) +
	labs(x = "Average UMAP Score", y = "Average Read Depth")

#Keep my workspace clean!
rm(reads_by_umap, plot_reads_by_umap)

#Sorting each column of map_gold to generate a percentile plot (from gold_sub, with only ~6000 points so my GPU survives)
gold_sort <- map_gold %>%
	select(gold:uniread) %>%
	mutate_each(list(~sort(.)))
gold_sub <- gold_sort %>%
	mutate(ptile = row_number()/length(gold)*100) %>%
	filter(row_number() %% 1000 == 0)
plot_gold_sub <- gold_sub %>%
	gather(variable, value, -ptile)
ggplot(data = plot_gold_sub, aes(x = ptile, y = value, color = variable)) + geom_line() + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 50), expand = c(0,0)) + 
	labs(x = "Percentile", y = "Read Depth")

#Keep my workspace clean!
rm(excess_sub, plot_excess_ptile)

#Define the color palette that I'll be using so I don't have to keep dealing with it.
scl_pal <- wes_palette("Zissou1", type="continuous")

#Generate several qq plots of different samples
plt_gold = gold_sort %>%
	filter(row_number() %% ceiling(length(gold)/100) == 0) %>%
	mutate(ptile = 1:99)

#QQ plot of Map 1 scored vs Iteration 1 Scored
ggplot(data = plt_gold, aes(x = uniread, y = multimap, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "Uniread Scored", y = "MultiMap Scored") 
ggsave(filename = "QQ-M1S-I1S_100bp.eps", width = 3, height = 3, units = "in")

#QQ plot of Map 1 scored vs Gold Scored
ggplot(data = plt_gold, aes(x = uniread, y = gold, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "Uniread Scored", y = "Gold Standard") 
ggsave(filename = "QQ-M1S-Gold_100bp.eps", width = 3, height = 3, units = "in")

#QQ plot of MultiMap scored vs Gold Scored
ggplot(data = plt_gold, aes(x = multimap, y = gold, color = ptile)) + geom_point() + geom_line() + 
	theme_classic() + geom_abline(slope = 1) +
	theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
	scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
	scale_x_continuous(limits = c(0, 40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = seq(0, 70, 10)) + 
	labs(x = "MultiMap Scored", y = "Gold Standard") 
ggsave(filename = "QQ-I1S-Gold_100bp.eps", width = 3, height = 3, units = "in")

#Keep my workspace clean!
rm(plt_gold)

#Generate plots of mean absolute error by UMAP
mae_by_UMAP = map_gold %>%
	group_by(UMAPbin) %>%
	transmute_at(vars(multimap:uniread), ~abs(.-gold)) %>%
	summarize_all(mean)
plot_mae_umap <- mae_by_UMAP %>%
	gather(variable, value, -UMAPbin)
ggplot(data = plot_mae_umap, aes(x = UMAPbin, y = value, color = variable)) + theme_classic() + geom_line() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	labs(xaxis = "Average UMAP Score", yaxis = "Excess Read Depth")
rm(mae_by_UMAP, plot_mae_umap)