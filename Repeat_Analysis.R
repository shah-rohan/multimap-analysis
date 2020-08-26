library(tidyverse)

#Since I have a header, for clarity: column names are Gene and a bunch of histone mods with either the suffix "I" or "M".
repdata = read.table("data/K562_Repeats_m100-to-p100_HistoneMods_iter1-map1.tab", header = T)

#Set to max 100 because above 100 is physically impossible for the data
rep100 = repdata %>%
	mutate_at(vars(-Gene), list(~ifelse(.>100, 100, .)))
line100 = repdata %>% 
	filter(grepl("LINE", Gene))
sine100 = repdata %>% 
	filter(grepl("SINE", Gene))
simprep100 = repdata %>% 
	filter(grepl("Simple", Gene))

#Principal component analysis using the columns with suffix "I"
rep100.pcaI = prcomp((rep100 %>% select(ends_with("I"))), scale = T, center = T)
line100.pcaI = prcomp((rep100 %>% select(ends_with("I"))), scale = T, center = T)
sine100.pcaI = prcomp((rep100 %>% select(ends_with("I"))), scale = T, center = T)
simprep100.pcaI = prcomp((rep100 %>% select(ends_with("I"))), scale = T, center = T)

#Optimizing the k-means clustering function -- using the "elbow" method
clustopt <- function(x)
{
	w <- vector()
	for (i in 1:10)
	{
		tws <- kmeans(x, centers = i, iter.max = 100, nstart = 20, algorithm = "Lloyd")$tot.withinss
		w <- c(w, tws)
	}
	w = as.data.frame(w)
	print(ggplot(data = w, aes(x = 1:10, y = w)) + geom_line())
	return(w)
}

#Running optimization on the clustering... this will take some time.
rep100.pcaI.wss <- clustopt(rep100.pcaI$x)
line100.pcaI.wss <- clustopt(line100.pcaI$x)
sine100.pcaI.wss <- clustopt(sine100.pcaI$x)
simprep100.pcaI.wss <- clustopt(simprep100.pcaI$x)

#For all groups, optimal number of clusters at 4 clusters.
rep100.pcaI.4k = kmeans(rep100.pcaI$x, centers = 4, iter.max = 100, nstart = 20, algorithm = "Lloyd")
line100.pcaI.4k = kmeans(line100.pcaI$x, centers = 4, iter.max = 100, nstart = 20, algorithm = "Lloyd")
sine100.pcaI.4k = kmeans(sine100.pcaI$x, centers = 4, iter.max = 100, nstart = 20, algorithm = "Lloyd")
simprep100.pcaI.4k = kmeans(simprep100.pcaI$x, centers = 4, iter.max = 100, nstart = 20, algorithm = "Lloyd")

#Get the average value of each of the histone modification columns for each cluster
rep100.pcaI.4k.avg = rep100 %>% 
	bind_cols(cluster = rep100.pcaI.4k$cluster) %>%
	group_by(cluster) %>% 
	summarize_at(vars(starts_with("K")), mean)
line100.pcaI.4k.avg = line100 %>% 
	bind_cols(cluster = line100.pcaI.4k$cluster) %>%
	group_by(cluster) %>% 
	summarize_at(vars(starts_with("K")), mean)
sine100.pcaI.4k.avg = sine100 %>% 
	bind_cols(cluster = sine100.pcaI.4k$cluster) %>%
	group_by(cluster) %>% 
	summarize_at(vars(starts_with("K")), mean)
simprep100.pcaI.4k.avg = simprep100 %>% 
	bind_cols(cluster = simprep100.pcaI.4k$cluster) %>%
	group_by(cluster) %>% 
	summarize_at(vars(starts_with("K")), mean)

#I'll use this one eventually... basically a post-hoc test for chi-squares.
pairwise.chisq.test <- function(x, met = "bonferroni"){
	R = nrow(x)
	C = ncol(x)
	
	ans = matrix(nrow = R, ncol = C)
	
	for (i in 1:R){
		for (j in 1:C){
			temp <- matrix(nrow = 2, ncol = 2)
			temp[1,1] = sum(x[i,j])
			temp[2,1] = sum(x[-i,j])
			temp[1,2] = sum(x[i,-j])
			temp[2,2] = sum(x[-i, -j])
			
			ans[i,j] = p.adjust(chisq.test(temp)$p.value, method = met, n = R*C)
		}
	}
	
	ans
}
