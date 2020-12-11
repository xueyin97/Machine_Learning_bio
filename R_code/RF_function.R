library(randomForest)
lung <- read.table(file.path("lung-cancer-sig-data.txt"),sep = "\t",header=T)
type <- read.table(file.path("LungCancer.SampleType.txt"),sep = "\t",header=T)

lung.data <- data.frame(t(lung[,3:122]))
lung.data$type <- type$Type

##1
n.data <- nrow(lung.data)
k <- 5
n.groups <- n.data%/%k
shu <- runif(n.data)
rk <- rank(shu)
fold <-(rk-1)%/%n.groups+1
fold <- as.factor(fold)


TP_cancer <- vector()
FP_cancer <- vector()
TN_cancer <- vector()
FN_cancer <- vector()

TP_normal <- vector()
FP_normal <- vector()
TN_normal <- vector()
FN_normal <- vector()


for(i in 1:5){
  rf.train <- randomForest(as.factor(type)~.,data = lung.data[fold!=i,])
  rf.pred <- predict(rf.train,lung.data[fold==i,-48])
  conf <- table(lung.data[fold==i,48],rf.pred)
  print(conf)
  
  TP_cancer[i] <- conf[1,1]
  FP_cancer[i] <- conf[1,2]
  TN_cancer[i] <- conf[2,2]
  FN_cancer[i] <- conf[2,1]
  
  TP_normal[i] <- conf[2,2]
  FP_normal[i] <- conf[2,1]
  TN_normal[i] <- conf[1,1]
  FN_normal[i] <- conf[1,2]
}

cancer_precision <- sum(TP_cancer/(TP_cancer+FP_cancer))/length(TP_cancer)
cancer_recall <- sum(TP_cancer/(TP_cancer+FN_cancer))/length(TP_cancer)
cancer_F1 <- 2*cancer_precision*cancer_recall/(cancer_precision+cancer_recall)

normal_precision <- sum(TP_normal/(TP_normal+FP_normal))/length(TP_normal)
normal_recall <- sum(TP_normal/(TP_normal+FN_normal))/length(TP_normal)
normal_F1 <- 2*cancer_precision*cancer_recall/(cancer_precision+cancer_recall)

cancer_eval <- data.frame('precision'=cancer_precision,'recall'=cancer_recall,'F1'=cancer_F1,row.names = "Cancer")
normal_eval <- data.frame('precision'=normal_precision,'recall'=normal_recall,'F1'=normal_F1,row.names = "Normal")

cancer_eval
normal_eval

##2
library(rpart)
test.cancer <- sample(1:60,12)
test.normal <- sample(61:120,12)
test.data <- lung.data[c(test.cancer,test.normal),]
train.data <- lung.data[-c(test.cancer,test.normal),]
my.rf <- function(train,test){
  dt.decision <- data.frame(matrix(ncol = 0, nrow = dim(test)[1]))
  rownames(dt.decision) <- rownames(test.data)
  for(i in 1:100){
    boots <- sample(1:dim(train)[1],dim(train)[1],replace=TRUE)
    features<- sample(1:(dim(train)[2]-1),round(1/3*(dim(train)[2]-1)))
    boots.sample <- train[boots,c(features,48)]
    dt<- rpart(type~.,method="class",data=boots.sample)
    dt.pred <- data.frame(predict(dt,test[,-dim(test)[2]]))
    decision <- ifelse(dt.pred$Cancer>0.05,"Cancer","Normal")
    dt.decision <- cbind(dt.decision,decision)
  }
  ##majority voting
  test.predict <- apply(dt.decision,1,function(x) names(which.max(table(x))))
  conf <- table(test[,48],test.predict)
  return(conf)
}
my.rf(train.data,test.data)
