library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)

source("ICeChIP_Analysis_Functions.R")

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

#Input analysis to compare different inputs to each other -- uniread vs. multimap
input_analysis(AR7$IR1I, AR7$IR1M, mm_UMAP$UMAPbin, "processed/AR7_Input_R1")
input_analysis(AR7$IR2I, AR7$IR2M, mm_UMAP$UMAPbin, "processed/AR7_Input_R2")
input_analysis(AR8$II, AR8$IM, dm_UMAP$UMAPbin, "processed/AR8_Input")
input_analysis(AR9$II, AR9$IM, mm_UMAP$UMAPbin, "processed/AR9_Input")
input_analysis(AR16$II, AR16$IM, hg_UMAP$UMAPbin, "processed/AR16_Input")
input_analysis(AR17$II, AR17$IM, hg_UMAP$UMAPbin, "processed/AR17_Input")

#Ratio analysis to compare biological replicates to each other.
ratio_analysis(AR7$IR1I, AR7$IR2I, 396109479, 387757427, "processed/AR7_Input_R2-AR7_Input_R1_multimap_logratio.tab")
ratio_analysis(AR7$IR1M, AR7$IR2M, 311090692, 304127899, "processed/AR7_Input_R2-AR7_Input_R1_uniread_logratio.tab")
ratio_analysis(AR7$IR1I, AR9$II, 396109479, 620463606, "processed/AR9_Input-AR7_Input_R1_multimap_logratio.tab")
ratio_analysis(AR7$IR1M, AR9$IM, 311090692, 488503092, "processed/AR9_Input-AR7_Input_R1_uniread_logratio.tab")
ratio_analysis(AR7$IR2I, AR9$II, 387757427, 620463606, "processed/AR9_Input-AR7_Input_R2_multimap_logratio.tab")
ratio_analysis(AR7$IR2M, AR9$IM, 304127899, 488503092, "processed/AR9_Input-AR7_Input_R2_uniread_logratio.tab")
ratio_analysis(AR16$II, AR17$II, 342591891, 305008807, "processed/AR17_Input-AR16_Input_multimap_logratio.tab")
ratio_analysis(AR16$IM, AR17$IM, 285996344, 256373920, "processed/AR17_Input-AR16_Input_uniread_logratio.tab")

#To this point, the calibrants were ignored. Now, I need to review them and normalize with them.
#First load in the calibration tables. AR7 and AR8 are easy, only one calibrant type
AR7_cal_table = AR7_cals %>%
	transmute(chr, mark = "H3K4me3")
AR8_cal_table = AR8_cals %>%
	transmute(chr, mark = "H3K27me3")
#Now load in AR9, AR16, and AR17 calibration tables. AR16 and AR17 are same.
AR9_cal_table = read.table("AR9_Calibration_Table.tab")
colnames(AR9_cal_table) = c("chr", "mark")

AR16_AR17_cal_table = read.table("AR16-AR17_Calibration_Table.tab")
colnames(AR16_AR17_cal_table) = c("chr", "mark")

#Getting the calibrant target ratios; it's easier to just manually input.
calibrant_analysis(AR7_cals$chr, AR7_cals$K4R1I, AR7_cals$IR1I, AR7_cal_table, "H3K4me3") # 20.0543
calibrant_analysis(AR7_cals$chr, AR7_cals$K4R1M, AR7_cals$IR1M, AR7_cal_table, "H3K4me3") # 19.87643
calibrant_analysis(AR7_cals$chr, AR7_cals$K4R2I, AR7_cals$IR2I, AR7_cal_table, "H3K4me3") # 18.99393
calibrant_analysis(AR7_cals$chr, AR7_cals$K4R2M, AR7_cals$IR2M, AR7_cal_table, "H3K4me3") # 18.9507
calibrant_analysis(AR8_cals$chr, AR8_cals$K27I, AR8_cals$II, AR8_cal_table, "H3K27me3") # 0.8789253
calibrant_analysis(AR8_cals$chr, AR8_cals$K27M, AR8_cals$IM, AR8_cal_table, "H3K27me3") # 0.8769298
calibrant_analysis(AR9_cals$chr, AR9_cals$K4I, AR9_cals$II, AR9_cal_table, "H3K4me3") # 28.31182
calibrant_analysis(AR9_cals$chr, AR9_cals$K4M, AR9_cals$IM, AR9_cal_table, "H3K4me3") # 27.66315
calibrant_analysis(AR9_cals$chr, AR9_cals$K9I, AR9_cals$II, AR9_cal_table, "H3K9me3") # 1.260581
calibrant_analysis(AR9_cals$chr, AR9_cals$K9M, AR9_cals$IM, AR9_cal_table, "H3K9me3") # 1.344
calibrant_analysis(AR9_cals$chr, AR9_cals$K27I, AR9_cals$II, AR9_cal_table, "H3K27me3") # 0.6774254
calibrant_analysis(AR9_cals$chr, AR9_cals$K27M, AR9_cals$IM, AR9_cal_table, "H3K27me3") # 0.6782338
calibrant_analysis(AR16_cals$chr, AR16_cals$me1I, AR16_cals$II, AR16_AR17_cal_table, "H3K4me1") # 4.843402
calibrant_analysis(AR16_cals$chr, AR16_cals$me1M, AR16_cals$IM, AR16_AR17_cal_table, "H3K4me1") # 4.343452
calibrant_analysis(AR16_cals$chr, AR16_cals$me2I, AR16_cals$II, AR16_AR17_cal_table, "H3K4me2") # 3.75365
calibrant_analysis(AR16_cals$chr, AR16_cals$me2M, AR16_cals$IM, AR16_AR17_cal_table, "H3K4me2") # 3.975822
calibrant_analysis(AR16_cals$chr, AR16_cals$me3I, AR16_cals$II, AR16_AR17_cal_table, "H3K4me3") # 31.09672
calibrant_analysis(AR16_cals$chr, AR16_cals$me3M, AR16_cals$IM, AR16_AR17_cal_table, "H3K4me3") # 32.43194
calibrant_analysis(AR17_cals$chr, AR17_cals$K9I, AR17_cals$II, AR16_AR17_cal_table, "H3K9me3") # 2.226229
calibrant_analysis(AR17_cals$chr, AR17_cals$K9M, AR17_cals$IM, AR16_AR17_cal_table, "H3K9me3") # 2.447601
calibrant_analysis(AR17_cals$chr, AR17_cals$K27I, AR17_cals$II, AR16_AR17_cal_table, "H3K27me3") # 1.73182
calibrant_analysis(AR17_cals$chr, AR17_cals$K27M, AR17_cals$IM, AR16_AR17_cal_table, "H3K27me3") # 1.824468

AR7_targ = c(20.0543, 19.87643, 18.99393, 18.9507)
AR8_targ = c(0.8789253, 0.8769298)
AR9_targ = c(28.31182, 27.66315, 1.260581, 1.344, 0.6774254, 0.6782338)
AR16_targ = c(4.843402, 4.343452, 3.75365, 3.975822, 31.09672, 32.43194)
AR17_targ = c(1.73182, 1.824468, 2.226229, 2.447601)

#Now compute HMD across windows
#AR7: This is incredibly ugly and a bit stupid but it works, and I'll take it
AR7_HMD_temp = (AR7[, 4:7]/AR7[, 8:11]) %>%
	mutate_all(~ ifelse(is.finite(.), ., 0))
AR7_HMD = bind_cols(AR7[,1:3], AR7_HMD_temp*100/AR7_targ)
rm(AR7_HMD_temp)

AR8_HMD = AR8 %>%
	mutate_at(vars(ends_with("I")), list(~ ./II)) %>%
	mutate_at(vars(ends_with("M")), list(~ ./IM)) %>%
	select(-II, -IM) %>%
	mutate_at(vars(-chr, -start, -stop), list(~ ifelse(is.finite(.), ., 0)))
AR8_HMD[, 4:5] = AR8_HMD[, 4:5]*100/AR8_targ

AR9_HMD = AR9 %>%
	mutate_at(vars(ends_with("I")), list(~ ./II)) %>%
	mutate_at(vars(ends_with("M")), list(~ ./IM)) %>%
	select(-II, -IM) %>%
	mutate_at(vars(-chr, -start, -stop), list(~ ifelse(is.finite(.), ., 0)))
AR9_HMD[, 4:9] = AR9_HMD[, 4:9]*100/AR9_targ

AR16_HMD = AR16 %>%
	mutate_at(vars(ends_with("I")), list(~ ./II)) %>%
	mutate_at(vars(ends_with("M")), list(~ ./IM)) %>%
	select(-II, -IM) %>%
	mutate_at(vars(-chr, -start, -stop), list(~ ifelse(is.finite(.), ., 0)))
AR16_HMD[, 4:9] = AR16_HMD[, 4:9]*100/AR16_targ

AR17_HMD = AR17 %>%
	mutate_at(vars(ends_with("I")), list(~ ./II)) %>%
	mutate_at(vars(ends_with("M")), list(~ ./IM)) %>%
	select(-II, -IM) %>%
	mutate_at(vars(-chr, -start, -stop), list(~ ifelse(is.finite(.), ., 0)))
AR17_HMD[, 4:7] = AR17_HMD[, 4:7]*100/AR17_targ