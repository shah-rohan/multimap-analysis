library(data.table)

for (i in 2:50){
startupstr = paste(c("Now starting map ", i), collapse = "")
print(startupstr)

strname = paste(c("hg38_6mil_200bp_30xCov_vf_k51_multiread_map-", i, ".bed"), collapse = "")
stroutname = paste(c("hg38_6mil_200bp_30xCov_vf_k51_multiread_randomread_map-", i, ".bed"), collapse = "")

data = fread(strname, sep = "\t", header = FALSE)
subdata = data[seq(0, dim(data)[1]/i-1, by = i) + as.numeric(sample.int(i, dim(data)[1]/i, replace = TRUE)),]

write.table(subdata, file = stroutname, quote = F, col.names = F, row.names = F, sep = "\t")
}
