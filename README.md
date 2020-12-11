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
### RF_funuction.R:
1. Do 5-fold cross-validation on Random Forest using the R “randomForest” package, and each time:
– Use 80% dataset as training data to train the model
– leave 20% as testing data, and use it to test the trained model
– Balance the numbers of Cancer and Control samples in both training and testing subsets
– Output the confusion matrix, the averaged precision, recall, and F-1 score for both Cancer and Control classes
2. Write R function to do Random Forest (1/3 features, 100 trees) using the R package “rpart”
