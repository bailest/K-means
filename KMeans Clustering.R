#after data is loaded and DGEList is defined on 20, the script will run until 83
#From there, you will input the amount of clusters on 88 and it will run from there
#insert filepath
setwd("~/")
library(edgeR)
library(SummarizedExperiment)
#Sample data set to test the script
#load(url("http://duffel.rail.bio/recount/SRP049355/rse_gene.Rdata"))
#counts <- assays(rse_gene)$counts

counts <- read.delim(".txt", row.names=1, header = TRUE)

y <- as.matrix((counts))
#Set the following to match your experimental setup. Like-numbers indicates they are 
#samples within the same group
#y <- DGEList(counts = y, group=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12))
#y <- DGEList(counts = y, group=c(1,2,3,4,5,6,7,8,9,10,11,12))
#y <- DGEList(counts = y, group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,8,8,8,9,9,10,10,11,11,11,12,12,12))
#y <- DGEList(counts = y, group=c(1,1,1,2,2,2,3,3,3,4,4,4))
y <-DGEList(counts = y, group=c(1,2,3,4))
#normalize
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)

#Data filtering
#Data 
#filtering based on mean expression and variance (dropping low counts and low variance)
z_var <- apply(z, 1, var)
z_mean <- apply(z, 1, mean)
plot(log2(z_mean), log2(z_var), pch='.')
abline(h=log2(50), col='red')
abline(v=log2(50), col='red')
text(x=13,y=23, labels="variance > 50 &\n mean > 50", col='red')
#plot shows the distribution of genes in the set based on their mean expression and variance
#We only want to look at the genes on the top right quadrant, next pulls those out
z <- z[which(z_var > 50 & z_mean > 50), 1:4]

#Scale the data because we want to identify clusters of genes that share similar 
#expression profiles rather than similar expression levels
scaledata <- t(scale(t(z)))

#clustering samples to identify outliers
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
TreeC = as.dendrogram(hc, method="average")

Sample_Cluster <- plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")


#choosing number of clusters
#There are multiple methods to choose the value of K
#sum of squared error (SSE)
#SSE is defined as the sum of the squared distance between each member of a cluster 
#and its cluster centroid. Basically, we repeatedly test this distance assigned as we 
#increase cluster amount, this distance (ie.similarity to centroid) decreases as we 
#increase the number of clusters made but eventually the distance does not significantly
#change between 2 cluster amounts. This is the point we want to choose for cluster amount

#in the following plot, we choose the elbow on the curve as the amount of clusters
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i,iter.max = 20)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#Average silhouette width
#silhouette value describes how similar a gene is to its own cluster (cohesion) 
#compared to other clusters (separation). A high value indicates that the 
#gene is well placed. So if the average of all of these silhouettes is high 
#then the number of clusters is good.
library(cluster)
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}
# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)



#performing K means clustering, change value of "centers" to match your amount of clusters
set.seed(20)
kClust <- kmeans(scaledata, centers=3, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster

#calculating centroids
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)

#plotting the centroids to see how they behave
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

#plot
position <- c("TP12_R_Mock", "TP24_R_Mock", "TP72_R_Mock", "TP12_R_TSWV",
              "TP24_R_TSWV", "TP72_R_TSWV","TP12_S_Mock","TP24_S_Mock",
              "TP72_S_Mock","TP12_S_TSWV","TP24_S_TSWV","TP72_S_TSWV")

p <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  scale_x_discrete(guide = guide_axis(n.dodge=3) ) +
      #             ,limits=position)+
  xlab("Group") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster Expression",color = "Cluster")
p

#calculate the correlation between the centroids to see thier similarity
#ideally, they are all dissimilar to one another
cor(kClustcentroids)

#Subset the cores molten dataframe so we can plot the core
#Note: if you are missing any of these cores, you don't need to comment them out,
#they will just produce an empty vector
core1 <- Kmolten[Kmolten$cluster=="1",]
core2 <- Kmolten[Kmolten$cluster=="2",]
core3 <- Kmolten[Kmolten$cluster=="3",]
core4 <- Kmolten[Kmolten$cluster=="4",]
core5 <- Kmolten[Kmolten$cluster=="5",]
#core6 <- Kmolten[Kmolten$cluster=="6",]

#get cluster 1
K1 <- (scaledata[kClusters==1,])
#calculate the correlation with the core
corscore1 <- function(x){cor(x,core1$value)}
score <- apply(K1, 1, corscore1)
#get the data frame into long format for plotting
K1molten <- melt(K1)
colnames(K1molten) <- c('gene','sample','value')
#add the score
K1molten <- merge(K1molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K1molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K1molten$order_factor <- 1:length(K1molten$gene)
#order the dataframe by score
K1molten <- K1molten[order(K1molten$score),]
#set the order by setting the factors
K1molten$order_factor <- factor(K1molten$order_factor , levels = K1molten$order_factor)

#get cluster 2
K2 <- (scaledata[kClusters==2,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)

#get cluster 3
K3 <- (scaledata[kClusters==3,])
#calculate the correlation with the core
corscore3 <- function(x){cor(x,core3$value)}
score <- apply(K3, 1, corscore)
#get the data frame into long format for plotting
K3molten <- melt(K3)
colnames(K3molten) <- c('gene','sample','value')
#add the score
K3molten <- merge(K3molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K3molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K3molten$order_factor <- 1:length(K3molten$gene)
#order the dataframe by score
K3molten <- K3molten[order(K3molten$score),]
#set the order by setting the factors
K3molten$order_factor <- factor(K3molten$order_factor , levels = K3molten$order_factor)

#get cluster 4
K4 <- (scaledata[kClusters==4,])
#calculate the correlation with the core
corscore4 <- function(x){cor(x,core4$value)}
score <- apply(K4, 1, corscore)
#get the data frame into long format for plotting
K4molten <- melt(K4)
colnames(K4molten) <- c('gene','sample','value')
#add the score
K4molten <- merge(K4molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K4molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K4molten$order_factor <- 1:length(K4molten$gene)
#order the dataframe by score
K4molten <- K4molten[order(K4molten$score),]
#set the order by setting the factors
K4molten$order_factor <- factor(K4molten$order_factor , levels = K3molten$order_factor)

#get cluster 5
K5 <- (scaledata[kClusters==5,])
#calculate the correlation with the core
corscore5 <- function(x){cor(x,core5$value)}
score <- apply(K5, 1, corscore)
#get the data frame into long format for plotting
K5molten <- melt(K5)
colnames(K5molten) <- c('gene','sample','value')
#add the score
K5molten <- merge(K5molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K5molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K5molten$order_factor <- 1:length(K5molten$gene)
#order the dataframe by score
K5molten <- K5molten[order(K5molten$score),]
#set the order by setting the factors
K5molten$order_factor <- factor(K5molten$order_factor , levels = K3molten$order_factor)

#Cluster 1 plot
p1 <- ggplot(K1molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  #this adds the core 
  geom_line(data=core1, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster 1 Expression",color = "Score")
p1

#Cluster 2 plot
p2 <- ggplot(K2molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  #this adds the core 
  geom_line(data=core2, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster 2 Expression",color = "Score")
p2

#Cluster 3 plot
p3 <- ggplot(K3molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  #this adds the core 
  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster 3 Expression",color = "Score")
p3

#Cluster 4 plot
p4 <- ggplot(K4molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  #this adds the core 
  geom_line(data=core4, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster 4 Expression",color = "Score")
p4

#Cluster 5 plot
p5 <- ggplot(K5molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  #this adds the core 
  geom_line(data=core5, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  theme_classic() +
  labs(title= "Cluster 5 Expression",color = "Score")
p5


#Export cluster lists into text files
write.table(K1molten, "k1_.txt",sep="\t",row.names=FALSE)
write.table(K2molten, "k2_.txt",sep="\t",row.names=FALSE)
write.table(K3molten, "k3_.txt",sep="\t",row.names=FALSE)
write.table(K4molten, "k4_.txt",sep="\t",row.names=FALSE)
write.table(K5molten, "k5_.txt",sep="\t",row.names=FALSE)

#Export plots
png(filename="allclusters.png", width=400, height=250, bg="white")
p
dev.off()

png(filename="k1.png", width=400, height=250, bg="white")
p1
dev.off()

png(filename="k2.png", width=400, height=250, bg="white")
p2
dev.off()

png(filename="k3.png", width=400, height=250, bg="white")
p3
dev.off()

png(filename="k4.png", width=400, height=250, bg="white")
p4
dev.off()

png(filename="k5.png", width=400, height=250, bg="white")
p5
dev.off()

png(filename="Kmeans_sample_clusters.png", width=900, height=1000, bg="white")
Sample_Cluster
dev.off()

#playing with subsetting based on score to see the genes that are most like the centroid
#If you create a small amount of clusters, likely many genes will have low scores
#relative to the centroid. In some case, you may want to ignore them and see how it looks
subset <- subset(K1molten, score > 0.9)

sub <- ggplot(subset, aes(x=sample,y=value)) + 
  geom_line(aes(color= "Red2", group=gene)) +
  scale_x_discrete(guide = guide_axis(n.dodge=3))+
  theme(legend.position="None") +
  #this adds the core 
  geom_line(data=core1, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE)
sub

png(filename="0.9cutoffKmean_example.png", width=400, height=250, bg="white")
sub
dev.off()
