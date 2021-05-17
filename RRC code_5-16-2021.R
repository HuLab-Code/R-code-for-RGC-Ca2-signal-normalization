
# Code for normalization

# First, read in all the data and named as 'rgcraw'
# Then, set the first column as time
# Because imageJ automatically export the frame number as 'Time', so we need to transform it first
# We got 283 frames in 60s, so we need to divivde by 4.7 to get the correct time coloumn
time <- round(rgcraw$Time*1/4.7, digit = 1) 
# delect the useless coloumns with the colnames as 'Time' and 'Stim'
rgcraw <- rgcraw[, -c(grep("Time", colnames(rgcraw)), grep("Stim", colnames(rgcraw)))]
# Third, calculate the df/f and add the time coloumn
# Use the mean value of top 30 frames as a baseline
rgccolmean30 <- colMeans(rgcraw[1:30,])
rgcsubmean <- sweep(rgcraw, 2, rgccolmean30, '-')
rgcnormalization <- sweep(rgcsubmean, 2, rgccolmean30, '/')
rgc <- cbind(time, rgcnormalization)
# finally, as the normalized data as normalized_rgc.csv
write.csv(rgc,"normalized_rgc.csv")

# ------------------------------------------------------
# ------------------------------------------------------
# code for different RGC subtype visualization 

# First, use Hierarchical Clustering (hclust, method = "ward.D2"), to cluster the RGCs, and use cutree to exrtact the cell reaction information of different clusters.
# The number of k is depend on how many cluster you can get.
# We can firstly set a k value based on the Hierarchical Clustering result, 
# and estimate if the k value is suitable, based on the downstream visualization of RGC reaction waveform.
# It need to adjust several times, and sometime we need to perform the sub cluster untill the cluster shows a single waveform.
# Then, combine the similar clusters based on the waveform character

rgc_t <- t(rgc[,-1])
rgc_hc <- hclust(dist(rgc_t), method="ward.D2")
plot(hc)
rect.hclust(rgc_hc, k=9, border="green") 
rgc_hc_k <- cutree(rgc_hc, k=9)

clusterN <- data.frame(rgc_hc_k)
clusterN$cell <- rownames(clusterN)
names(clusterN) <- c("number", "cell")
cluster1 <- clusterN[clusterN$number=="1",2]
cluster2 <- clusterN[clusterN$number=="2",2]
cluster3 <- clusterN[clusterN$number=="3",2]
cluster4 <- clusterN[clusterN$number=="4",2]
cluster5 <- clusterN[clusterN$number=="5",2]
cluster6 <- clusterN[clusterN$number=="6",2]

react_1 <- cbind(time, rgc_reaction_n[,cluster1_c])
react_2 <- cbind(time, rgc_reaction_n[,cluster2_c])
react_3 <- cbind(time, rgc_reaction_n[,cluster3_c])
react_4 <- cbind(time, rgc_reaction_n[,cluster4_c])
react_5 <- cbind(time, rgc_reaction_n[,cluster5_c])
react_6 <- cbind(time, rgc_reaction_n[,cluster6_c])

# use ggplot to visualize the RGC 
# take cluster1 as an example
tiff(file="1m.tif", res=150, units='in', width=10, height=10)
react_1m <- melt(react_1, id="time")
rgc_on <- ggplot(data=react_1m,
                 mapping=aes(x=time, y=value, cluster=variable, colour=variable)) 
rgc_on + geom_line(alpha=0.5, size=1) + 
  scale_y_continuous(limits=c(-1, 3), breaks=seq(-1, 3, 1), labels = percent) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_base() +
  theme(legend.position="none") + 
  labs(x ="Time/s", y ="∆F/F0") +
  theme(axis.text.x=element_text(family="Helvetica", color="black", size=30),
        axis.text.y=element_text(family="Helvetica", color="black", size=30),
        axis.title.x=element_text(family="Helvetica", size=30), 
        axis.title.y=element_text(family="Helvetica", size=30))
dev.off()

# combine the cluster, which belong to the same functional type
# group_on_sustained1 <- cbind( )
# group_on_sustained2 <- cbind( )
# group_on_sustained3 <- cbind( )
# group_on_transient <- cbind( )
# group_off_sustained <- cbind( )
# group_off_transient <- cbind( )
# group_on_off <- cbind( )
# group_suppression1 <- cbind( )
# group_suppression2 <- cbind( )
# no_response <- cbind( )

group_on_sustained1 <- cbind(react_1, react_2[,-1], react_3[,-1])

# visualize the cluster of on_transient
group_on_transient_mean <- cbind(group_on_transient[,1], rowMeans(group_on_transient[,-1]), group_on_transient[,-1]) 
names(group_on_transient_mean)[1:2] <- c("time", "group_on_transient_mean")
group_on_transient_meanm <- melt(group_on_transient_mean, id="time")
group_on_transientg <- data.frame(group_on_transient_meanm %>% filter(variable=="group_on_transient_mean"))
group_on_transientc <- group_on_transient_meanm[!group_on_transient_meanm$variable=="group_on_transient_mean",] #use to draw all grey lines
ggplot(group_on_transient_meanm, aes(x = time, y = value, group = variable)) + 
  geom_line(data = group_on_transientc, alpha = 0.5, size = 1, color="lightgrey") +
  geom_line(data = group_on_transientg, alpha = 1, size = 2, color="red") +
  scale_y_continuous(limits=c(-1, 3), breaks=seq(-1, 3, 1), labels = percent) + 
  scale_x_continuous(expand = c(0, 0)) +
  theme_base() +
  theme(legend.position="none") + 
  labs(x = "Time/s", y = "dF/F0") +
  theme(axis.text.x = element_text(family="Helvetica", color="black", size=30),
        axis.text.y = element_text(family="Helvetica", color="black", size=30),
        axis.title.x=element_text(family="Helvetica", size=30), 
        axis.title.y=element_text(family="Helvetica", size=30))
ggsave("group_on_transient.tiff", dpi=300, units='in', width=15, height=10 )





