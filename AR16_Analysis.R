library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)

scl_pal <- wes_palette("Zissou1", type="continuous")

#Read in the files for the datasets
AR7 = read.table("data/AR7_mESC_windows_H3K4me3-S1-S2-Input-S1-S2_iter1-map1.tab")
colnames(AR7) = c("chr", "start", "stop", "K4R1I", "K4R1M", "K4R2I", "K4R2M", "IR1I", "IR1M", "IR2I", "IR2M")
AR7_cals = AR7 %>%
	filter(row_number() <= 11)
AR7 = AR7 %>%
	filter(row_number() > 11)

AR8 = read.table("data/AR8_S2_windows_H3K27me3-Input_iter1-map1.tab")
colnames(AR8) = c("chr", "start", "stop", "K27I", "K27M", "II", "IM")
AR8_cals = AR8 %>%
	filter(row_number() <= 100)
AR8 = AR8 %>%
	filter(row_number() > 100)

AR9 = read.table("data/AR9_mESC_windows_H3K4me3-H3K9me3-H3K27me3-Input_iter1-map1.tab")
colnames(AR9) = c("chr", "start", "stop", "K4I", "K4M", "K9I", "K9M", "K27I", "K27M", "II", "IM")
AR9_cals = AR9 %>%
	filter(row_number() <= 100)
AR9 = AR9 %>%
	filter(row_number() > 100)

AR16 = read.table("data/AR16_K562_windows_H3K4me1-2-3-Input_iter1-map1.tab")
colnames(AR16) = c("chr", "start", "stop", "me1I", "me1M", "me2I", "me2M", "me3I", "me3M", "II", "IM")
AR16_cals = AR16 %>%
	filter(row_number() <= 136)
AR16 = AR16 %>% 
	filter(row_number() > 136)

AR17 = read.table("data/AR17_K562_windows_counts_H3K27me3-H3K9me3-Input_iter-map.tab")
colnames(AR17) = c("chr", "start", "stop", "K27I", "K27M", "K9I", "K9M", "II", "IM")
AR17_cals = AR17 %>%
	filter(row_number() <= 136)
AR17 = AR17 %>%
	filter(row_number() > 136)

#Load in the genome-wide UMAP files
hg_UMAP = read.table("UMAPS/hg38_UMAP50_windows.tab")
colnames(hg_UMAP) = c("chr", "start", "stop", "UMAP50")
hg_UMAP = hg_UMAP %>%
	mutate(UMAPbin = cut(UMAP50, seq(0, 1.01, by = 0.01), include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01))) %>%
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin])

dm_UMAP = read.table("UMAPS/dm3_UMAP50_windows.tab")
colnames(dm_UMAP) = c("chr", "start", "stop", "UMAP50")
dm_UMAP = dm_UMAP %>%
	mutate(UMAPbin = cut(UMAP50, seq(0, 1.01, by = 0.01), include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01))) %>%
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin])

mm_UMAP = read.table("UMAPS/mm10_UMAP50_windows.tab")
colnames(mm_UMAP) = c("chr", "start", "stop", "UMAP50")
mm_UMAP = mm_UMAP %>%
	mutate(UMAPbin = cut(UMAP50, seq(0, 1.01, by = 0.01), include.lowest = TRUE, right = FALSE, labels = seq(0, 1, by=0.01))) %>%
	mutate(UMAPbin = as.numeric(levels(UMAPbin))[UMAPbin])

#UMAP score density plot -- more rigorous, I suppose.
hg_dens = ggplot(hg_UMAP, aes(x = UMAP50)) + geom_density(bw = "nrd") + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(expand = c(0,0)) + 
	labs(x = "UMAP Score", y = "Density", title = "hg38")

dm_dens = ggplot(dm_UMAP, aes(x = UMAP50)) + geom_density(bw = "nrd0") + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(expand = c(0,0)) + 
	labs(x = "UMAP Score", y = "Density", title = "dm3")
	
mm_dens = ggplot(mm_UMAP, aes(x = UMAP50)) + geom_density(bw = "nrd") + 
	theme_classic() + 
	scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) +
	scale_y_log10(expand = c(0,0)) + 
	labs(x = "UMAP Score", y = "Density", title = "mm10")

print(hg_dens)
print(dm_dens)
print(mm_dens)
rm(hg_dens, dm_dens, mm_dens)

#Create a function to iterate over all the input analyses so I don't have to :)
input_analysis <- function(iter1, map1, umapbins, prefix){
	Input = bind_cols(as.data.frame(iter1), as.data.frame(map1), as.data.frame(umapbins))
	colnames(Input) = c("II", "IM", "UMAP")
	
	#Filenames
	percfile = paste(prefix, "_Reads_Percentile.tab", sep = "")
	umapfile = paste(prefix, "_Reads_by_UMAP.tab", sep = "")
	excessfile = paste(prefix, "_Excess_Percentile.tab", sep = "")
	excessumapfile = paste(prefix, "_Excess_by_UMAP.tab", sep = "")
	
	Input_sorted = Input %>%
		select(II, IM) %>%
		mutate_each(list(~sort(.)))
	Input_sorted_sub = Input_sorted %>%
		mutate(ptile = row_number()/length(II)*100) %>%
		filter(row_number() %% ceiling(length(II)/6000) == 0)
	Input_perc_plt = Input_sorted_sub %>%
		gather(variable, value, -ptile)
	Perc_plt <- ggplot(data = Input_perc_plt, aes(x = ptile, y = value, color = variable)) + geom_line() + 
		theme_classic() + 
		scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) + 
		scale_y_continuous(limits = c(0, 80), expand = c(0,0)) + 
		labs(x = "Percentile", y = "Read Depth", title = prefix)
	print(Perc_plt)
	write.table(Input_sorted_sub, percfile, sep = "\t", quote = F, row.names = F, col.names = T)
	
	Input_UMAP_binned = Input %>%
		group_by(UMAP) %>%
		summarize_at(vars(II:IM), list(~median(.)))
	Input_UMAP_plt = Input_UMAP_binned %>%
		gather(variable, value, -UMAP)
	UMAP_plt <- ggplot(data = Input_UMAP_plt, aes(x = UMAP, y = value, color = variable)) + geom_line() +
		theme_classic() + 
		scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
		scale_y_continuous(expand = c(0,0)) +
		labs(x = "UMAP Score", y = "Median Read Depth", title = prefix)
	print(UMAP_plt)
	write.table(Input_UMAP_binned, umapfile, sep = "\t", quote = F, row.names = F, col.names = T)
	
	Input_excess = Input %>%
		transmute(excess = II-IM, UMAP)
	
	Input_excess_sorted = Input_excess %>%
		select(excess) %>%
		mutate_each(list(~sort(.)))
	Input_excess_sorted_sub = Input_excess_sorted %>%
		mutate(ptile = row_number()/length(excess)*100) %>%
		filter(row_number() %% ceiling(length(excess)/6000) == 0)
	Input_excess_perc_plt = Input_excess_sorted_sub %>%
		gather(variable, value, -ptile)
	Excess_perc_plt <- ggplot(data = Input_excess_perc_plt, aes(x = ptile, y = value, color = variable)) + geom_line() + 
		theme_classic() + 
		scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) + 
		scale_y_continuous(limits = c(0, 80), expand = c(0,0)) + 
		labs(x = "Percentile", y = "Excess Read Depth", title = prefix)
	print(Excess_perc_plt)
	write.table(Input_excess_sorted_sub, excessfile, sep = "\t", quote = F, row.names = F, col.names = T)
	
	Input_excess_UMAP_binned = Input_excess %>%
		group_by(UMAP) %>%
		summarize_at(vars(excess), list(~median(.)))
	Input_excess_UMAP_plt = Input_excess_UMAP_binned %>%
		gather(variable, value, -UMAP)
	Excess_UMAP_plt <- ggplot(data = Input_excess_UMAP_plt, aes(x = UMAP, y = value, color = variable)) + geom_line() +
		theme_classic() + 
		scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0,0)) + 
		scale_y_continuous(expand = c(0,0)) +
		labs(x = "UMAP Score", y = "Median Excess Read Depth", title = prefix)
	print(Excess_UMAP_plt)
	write.table(Input_excess_UMAP_binned, excessumapfile, sep = "\t", quote = F, row.names = F, col.names = T)
	
	#QQ plot
	Input_qq = Input_sorted %>%
		filter(row_number() %% ceiling(length(II)/100) == 0) %>%
		mutate(ptile = 1:99)
	qq_plt = ggplot(data = Input_qq, aes(x = IM, y = II, color = ptile)) + geom_point() + geom_line() + 
		theme_classic() + geom_abline(slope = 1) +
		theme(text = element_text(size = 10), legend.position = c(0.05, 0.95), legend.justification = c("left", "top")) + 
		scale_color_gradientn(colors = scl_pal, name = "Percentile", limits = c(0, 100), position = "left") + 
		scale_x_continuous(limits = c(0, 60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
		scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = seq(0, 70, 10)) + 
		labs(x = "Uniread Scored", y = "Multimap Scored")
	print(qq_plt)
}

input_analysis(AR7$IR1I, AR7$IR1M, mm_UMAP$UMAPbin, "processed/AR7_Input_R1")
input_analysis(AR7$IR2I, AR7$IR2M, mm_UMAP$UMAPbin, "processed/AR7_Input_R2")
input_analysis(AR8$II, AR8$IM, dm_UMAP$UMAPbin, "processed/AR8_Input")
input_analysis(AR9$II, AR9$IM, mm_UMAP$UMAPbin, "processed/AR9_Input")
input_analysis(AR16$II, AR16$IM, hg_UMAP$UMAPbin, "processed/AR16_Input")
input_analysis(AR17$II, AR17$IM, hg_UMAP$UMAPbin, "processed/AR17_Input")

ratio_analysis <- function(input1, input2, reads_input1, reads_input2, filename) {
	Input_ratio = bind_cols(as.data.frame(input1), as.data.frame(input2))
	colnames(Input_ratio) = c("In1", "In2")
	Input_diff = Input_ratio %>%
		filter(In1>0 & In2>0) %>%
		mutate_each(list(~log10(.))) %>%
		transmute(dif = In2-In1 - log10(reads_input2/reads_input1)) %>%
		arrange(dif) %>%
		mutate(ptile = row_number()/length(dif)*100) %>%
		filter(row_number() %% ceiling(length(dif)/6000) == 0)
	mae = Input_diff %>%
		select(dif) %>%
		abs %>%
		colMeans
	print(mae)
	write.table(Input_diff, filename, sep = "\t", quote = F, row.names = F, col.names = T)
}

ratio_analysis(AR7$IR1I, AR7$IR2I, 396109479, 387757427, "processed/AR7_Input_R2-AR7_Input_R1_multimap_logratio.tab")
ratio_analysis(AR7$IR1M, AR7$IR2M, 311090692, 304127899, "processed/AR7_Input_R2-AR7_Input_R1_uniread_logratio.tab")
ratio_analysis(AR7$IR1I, AR9$II, 396109479, 620463606, "processed/AR9_Input-AR7_Input_R1_multimap_logratio.tab")
ratio_analysis(AR7$IR1M, AR9$IM, 311090692, 488503092, "processed/AR9_Input-AR7_Input_R1_uniread_logratio.tab")
ratio_analysis(AR7$IR2I, AR9$II, 387757427, 620463606, "processed/AR9_Input-AR7_Input_R2_multimap_logratio.tab")
ratio_analysis(AR7$IR2M, AR9$IM, 304127899, 488503092, "processed/AR9_Input-AR7_Input_R2_uniread_logratio.tab")
ratio_analysis(AR16$II, AR17$II, 342591891, 305008807, "processed/AR17_Input-AR16_Input_multimap_logratio.tab")
ratio_analysis(AR16$IM, AR17$IM, 285996344, 256373920, "processed/AR17_Input-AR16_Input_uniread_logratio.tab")

#To this point, the calibrants were ignored. Now, I need to review them and normalize with them.