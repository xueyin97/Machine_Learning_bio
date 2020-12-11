# Machine Learning on the microarray dataset of Lung Cancer
## Python_code
### significant_genes_selection.py:
Select significant genes by Row t-test for Every Gene
1. Get the raw p-value and Benjamini-Hochberg False Positive Rate (BH-FDR) Adjusted p-value, Insert them into the first two columns of the original data
2. Use the BH-FDR Adjusted p-value <0.05 as cutoff to take significant gene list, and export the data into one txt file.
### ML_models.py:
The First 60 Columns are Cancer Samples and the Second 60 Columns are Control Samples.
1. Run 5-fold Cross-Validation using Following Method, and Report Average F1 Scores
– Logistic Regression
– Support Vector Machine
– Random Forest
2. Use 100% data to build an Ensemble Learning Model (soft voting) using Logistic Regression, Support Vector Machine, Random Forest
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
### SVM_function.R:
Write R Function to do Step-down (Backward) Feature Selection using SVM:
Step-down algorithm:
Let F = {all features}
While not reduced to desired number of features
For each feature f ∈ F:
Estimate model’s performance on feature set F - f (using cross-validation)
Remove the f that leads to the best performed model with (F-f)
