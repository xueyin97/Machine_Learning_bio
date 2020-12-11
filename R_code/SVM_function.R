library(glmnet)
library(e1071)

lung <- read.table(file.path("lung-cancer-sig-data.txt"),sep = "\t",header=T)
label <- read.table(file.path("LungCancer.SampleType.txt"),sep = "\t",header=T)

svm.data <- as.data.frame(t(lung[,3:122]))
colnames(svm.data)<-lung[,1]

##My function to do stepdown feature selection using Support Vector Machine(SVM)
mySVM <- function(X,y,N){
  X$type <- y
  n.data <- nrow(X)
  k <- 5
  n.groups <- n.data%/%k
  shu <- runif(n.data)
  rk <- rank(shu)
  fold <-(rk-1)%/%n.groups+1
  fold <- as.factor(fold)
  
  while(dim(X)[2]-1>N){
    accuracy_score <- vector()
    for (f in 1:sum(dim(X)[2]-1)){
      X.data <- X[,-f]
      TP_cancer <- vector()
      FP_cancer <- vector()
      TN_cancer <- vector()
      FN_cancer <- vector()
      for(i in 1:5){
        svm.model <- svm(as.factor(type)~.,data = X.data[fold!=i,],cost=100,gamma=1)
        svm.pred <- predict(svm.model,X.data[fold==i,-dim(X.data)[2]],)
        conf <- table(X.data[fold==i,dim(X.data)[2]],svm.pred)
        
        TP_cancer[i] <- conf[1,1]
        FP_cancer[i] <- conf[1,2]
        TN_cancer[i] <- conf[2,2]
        FN_cancer[i] <- conf[2,1]
        }
      accuracy_score[f] <- sum(TP_cancer+TN_cancer/(TP_cancer+TN_cancer+FP_cancer+FN_cancer))/length(TP_cancer)
    }
    X <- X[,-which.max(accuracy_score)]
  }
  return.model <- svm(as.factor(type)~.,data=X,cost=100,gamma=1)
  return.list <- list("final model"=return.model,"selected features"=names(X[,-dim(X)[2]]))
  return(return.list)
}

#test mySVM function(selecting 46/45/44 features)
mySVM(svm.data,label$Type,46)
mySVM(svm.data,label$Type,45)
mySVM(svm.data,label$Type,44)
