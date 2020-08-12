library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)
library(clipr)

scl_pal <- wes_palette("Zissou1", type="continuous")

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

calibrant_analysis <- function(cals_chrs, cals_IP, cals_In, cal_table, mark){
	cals_total = as.data.frame(cals_chrs) %>%
		bind_cols(as.data.frame(cals_IP), as.data.frame(cals_In)) %>%
		inner_join(cal_table, by = c("cals_chrs" = "chr"))
	colnames(cals_total) = c("chr", "IP", "Input", "Mark")
	cals_by_mark = cals_total %>%
		group_by(Mark) %>%
		summarize_at(vars(IP:Input), list(~sum(.)))
	write_clip(cals_by_mark %>%
		mutate(ratio = IP/Input) %>%
		mutate(spec = 100*ratio/ratio[Mark == mark]))
	cal_mark = cals_by_mark %>%
		filter(Mark == mark)
	return(cal_mark$IP/cal_mark$Input)
}