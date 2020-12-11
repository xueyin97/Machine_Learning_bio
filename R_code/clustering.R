# hierarchical clustering
lung<-read.table(file.path("lung-cancer-sig-data.txt"),sep = "\t",header=T)
lung.fold<-as.data.frame(t(lung[,2:122]))
lung.fold<-lung.fold[-1,]
colnames(lung.fold)<-lung[,2]
library(e1071)
library(fpc)
lung.fold<-as.data.frame(apply(lung.fold,2,as.numeric))
d<-dist(lung.fold,method="euclidean")
fit<-hclust(d,method="ward.D")
plot(fit)
groups<-cutree(fit,k=2)
rect.hclust(fit,k=2,border="red")

# heatmap
a<-as.matrix(lung.fold)
row.names(a)<-colnames(lung)[3:122]
colnames(a)<-c(colnames(a)[1:34],"",colnames(a)[36:47])
b<-heatmap(a,distfun = dist,hclustfun = hclust)
col<-vector()
lab<-vector()
for(i in 1:120)
{if(b$rowInd[i]<=60)
{col[i]<-"#FF0000"
lab[i]<-"Cancer"
}
  else
  {col[i]<-"#0000FF"
  lab[i]<-"Ctl"}
}
heatmap(a,distfun = dist,hclustfun = hclust,RowSideColors = col,cexRow = 0.7,labRow = lab)

# K-means clustering
## Use 2 to be the initial number of clusters to do the k-Means cluster
km<-kmeans(a,2)
km.cluster<-as.data.frame(km$cluster)
km.centers<-as.data.frame(km$centers)
write.table(km.cluster,file = "cluster2.txt")
write.table(km.centers,file ="centers2.txt" )

par(mfrow=c(1,2))
plot(1:47,km.centers[1,],pch=16,main="Profile of the first cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[1,],pch=16)
plot(1:47,km.centers[2,],pch=16,main="Profile of the second cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[2,],pch=16)

pca.fit<-prcomp(a)
summary(pca.fit)
colnames(km.cluster)<-"km_cluster"
shape=km.cluster$km_cluster
col=c(rep("red",60),rep("green",60))
par(mfrow=c(1,1))
plot(pca.fit$x[,1],pca.fit$x[,2],col=col,pch=shape,xlab = "First Component",ylab = "Second Component")
legend(-19,15,legend = c("Cancer","CTL"),col = c("red","green"),cex = 0.8,lwd = 1, lty = 1)
legend(-19,11.5,legend=c("Cluster1","Cluster2"),pch = c(1,2,1,2))

## Use 4 to be the initial number of clusters to do the cluster
km<-kmeans(a,4)
km.cluster<-as.data.frame(km$cluster)
km.centers<-as.data.frame(km$centers)
write.table(km.cluster,file = "cluster4.txt")
write.table(km.centers,file ="centers4.txt" )

par(mfrow=c(2,2))
plot(1:47,km.centers[1,],pch=16,main="Profile of the first cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[1,],pch=16)
plot(1:47,km.centers[2,],pch=16,main="Profile of the second cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[2,],pch=16)
plot(1:47,km.centers[3,],pch=16,main="Profile of the third cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[3,],pch=16)
plot(1:47,km.centers[4,],pch=16,main="Profile of the fourth cluster's centroid",xlab = "Genes",ylab = "Expression Values")
lines(1:47,km.centers[4,],pch=16)

pca.fit<-prcomp(a)
summary(pca.fit)
colnames(km.cluster)<-"km_cluster"
shape=km.cluster$km_cluster
col=c(rep("red",60),rep("green",60))
par(mfrow=c(1,1))
plot(pca.fit$x[,1],pca.fit$x[,2],col=col,pch=shape,xlab = "First Component",ylab = "Second Component")
legend(-19,15,legend = c("Cancer","CTL"),col = c("red","green"),cex = 0.8,
       lwd = 1, lty = 1)
legend(-19,11.5,legend=c("Cluster1","Cluster2","Cluster3","Cluster4"),pch = c(1,2,3,4))


# My K-median clustering
manhat <- function(centroid, dataMatrix) {
  distance <- matrix(NA, nrow=dim(centroid)[1], ncol=dim(dataMatrix)[1])
  for(i in 1:nrow(dataMatrix)) {
    distance[,i] <- rowSums(abs(t(t(centroid)-dataMatrix[i,])))
  }
  distance
}

K_median <- function(x, centers, n.iter) {
  last.clusters<-vector()
  last.centers<-matrix(0,nrow=dim(centers)[1],ncol=dim(centers)[2])
  flag<-0
  i=1
  while(sum(i >= n.iter ,flag==1)<1) {
    distsToCenters <- manhat(centers,x)
    clusters <- apply(distsToCenters, 2, which.min)
    centers <- apply(x, 2, tapply, clusters, median)
    if(all(centers==last.centers)&all(clusters==last.clusters)){
      flag<-1}
    last.clusters<-clusters
    last.centers<-centers
    i=i+1
  }
  list(kclusters=clusters, kcenters=centers)
}
kcenters<-a[c(10,50,100),]
K_median(a,kcenters,5)
