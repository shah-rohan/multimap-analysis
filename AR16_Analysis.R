library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)

AR16 = read.table("AR16_K562_windows_H3K4me1-2-3-Input_iter1-map1.tab")
colnames(AR16) = c("chr", "start", "stop", "me1I", "me1M", "me2I", "me2M", "me3I", "me3M", "II", "IM")
AR16_cals = AR16 %>% 
	filter(row_number() <= 136)
AR16 = AR16 %>% 
	filter(row_number() > 136)

hg_UMAP = read.table("hg38_UMAP50_windows.tab")
colnames(hg_UMAP) = c("chr", "start", "stop", "UMAP50")

#Check the genome-wide binning of UMAP50
UMAPbins = seq(0, 1.01, by = 0.01)
hg_UMAP_binned = hg_UMAP %>% 
	mutate(UMAPbin = cut(UMAP50, UMAPbins, include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01))) %>% 
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin]) 
UMAP_bins = hg_UMAP_binned %>%
	count(UMAPbin)
ggplot(data = UMAP_bins, aes(x = UMAPbin, y = n)) + geom_line() +
	theme_classic() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(limits = c(1e3, 1e8), expand = c(0,0), breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)) +
	labs(x = "UMAP Score", y = "Number of Windows")

#Binning density plot -- more rigorous, I suppose.
hg_dens = density(hg_UMAP_binned$UMAP50, from = 0, to = 1, bw = "nrd", n = 200)
hg_dens2 = as.data.frame(dens$x) %>% bind_cols(as.data.frame(dens$y))
colnames(hg_dens2) = c("x", "y")
ggplot(hg_dens2, aes(x = x, y = y)) + geom_line() + scale_y_log10()

#Create UMAP proportion plot
UMAP_prop = UMAP_bins %>%
	mutate(n = n/sum(n))
ggplot(data = UMAP_prop, aes(x = UMAPbin, y = n)) + geom_line() +
	theme_classic() +
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(limits = c(1e-4, 1), expand = c(0,0), breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
	labs(x = "UMAP Score", y = "Proportion of Windows")

#Keep my workspace clean!
rm(UMAP_bins, UMAP_prop, UMAPbins, hg_UMAP)

#Percentile plot of read depth across the windows for the Input, iter-1 vs. map-1
AR16_Input = AR16 %>%
	select(II:IM)
AR16_Input_sorted = AR16_Input %>%
	mutate_each(list(~sort(.)))
AR16_in_sort_sub = AR16_Input_sorted %>%
	mutate(ptile = row_number()/length(II)*100) %>%
	filter(row_number() %% 2500 == 0)
plt_AR16_in_sort_sub = AR16_in_sort_sub %>%
	gather(variable, value, -ptile)
ggplot(data = plt_AR16_in_sort_sub, aes(x = ptile, y = value, color = variable)) + geom_line() + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 80), expand = c(0,0)) + 
	labs(x = "Percentile", y = "Read Depth")

rm(AR16_in_sort_sub, plt_AR16_in_sort_sub)

#Read depth by UMAP score
AR16_UMAP_Input = hg_UMAP_binned %>% select(UMAP50, UMAPbin) %>% bind_cols(AR16_Input)
AR16_Input_UMAPbinned = AR16_UMAP_Input %>%
	group_by(UMAPbin) %>%
	summarize_at(vars(II:IM), list(~mean(.)))
AR16_Input_UB_ci = AR16_UMAP_Input %>%
	group_by(UMAPbin) %>%
	summarize_at(vars(II:IM), list(~ 1.96*sd(.)/sqrt(n())))
AR16_Input_shade = bind_cols(AR16_Input_UMAPbinned - AR16_Input_UB_ci, AR16_Input_UMAPbinned + AR16_Input_UB_ci) %>%
	transmute(UMAPbin = UMAPbin1/2, II_min = II, IM_min = IM, II_max = II1, IM_max = IM1)
AR16_reads_by_umap = AR16_Input_UMAPbinned %>%
	gather(variable, value, -UMAPbin)
ggplot(data = AR16_reads_by_umap, aes(x = UMAPbin, y = value, color = variable)) + geom_line() +
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
	scale_y_continuous(limits = c(0, 25), expand = c(0,0)) +
	labs(x = "Average UMAP Score", y = "Average Read Depth")

rm(AR16_Input_UMAPbinned, AR16_Input_UB_ci, AR16_Input_shade, AR16_reads_by_umap)

AR16_Input_excess  = AR16_UMAP_Input %>%
	transmute(UMAP50, UMAPbin, excess = II-IM)
AR16_excess_sorted = AR16_Input_excess %>%
	mutate_each(list(~sort(.)))
