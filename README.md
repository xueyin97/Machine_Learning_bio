# Machine Learning on the microarray dataset of Lung Cancer
## Python_code

## R_code
### Clustering.R: 
1. Do hierarchical clustering on samples, and cut the dendrogram to make samples into areasonable number of clusters
2. Generate a bi-clustering heatmap
3. Use K-means clustering to do clustering on the significant gene list by samples
4. Write R function to implement the k-median algorithm (Use median instead of mean to calculate the cluster centers)
– Use Manhattan instead of Euclidian to calculate distances
– Take the 3 arguments
• The data matrix
• The initial centers
• The number of iterations (Pre-defined by user, so that it can still finish even if the search does not converge)
– Return 2 types of Result Values
• Cluster Assignment for each Data Record
• Centers of each Cluster
