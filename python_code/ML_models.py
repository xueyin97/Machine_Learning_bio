import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import LinearSVC
from sklearn import model_selection
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn import metrics

#Run 5-fold Cross-Validation using Following Method, and Report Average F1 Scores
def CV(data,label):
 logit = LogisticRegression(solver = "sag", dual = False, max_iter=10000)
 L1SVC = LinearSVC(loss = 'squared_hinge', penalty = 'l1', dual = False, max_iter=10000)
 randomForest = RandomForestClassifier(n_estimators = 100)
 classifiers = [logit, L1SVC, randomForest]
 classifier_names = ['logit', 'L1SVC', 'randomForest']
 kf = StratifiedKFold(n_splits=5)
 i_count = 1
 for train_index, test_index in kf.split(data, label):
  X_train = [data.iloc[i,:] for i in train_index]
  X_test = [data.iloc[i,:] for i in test_index]
  y_train, y_test = label[train_index], label[test_index]
  resultF1 = {}
  for j, clf in enumerate(classifiers):	
   clf.fit(X_train, y_train)
   pred = clf.predict(X_test)
   score = metrics.f1_score(y_test, pred, average='macro')
   resultF1[classifier_names[j]] = score
  print("Average F1 Scores")
  print(resultF1)
  i_count+=1

#Use 100% data to build an Ensemble Learning Model (soft voting) using
def emsemble(data,label):
 #emsemble learning model
 kfold = model_selection.KFold(n_splits=5)
 estimators = []
 model1 = LogisticRegression(solver = "sag", dual = False, max_iter=10000)
 estimators.append(('Logistic', model1))
 model2 = LinearSVC(loss = 'squared_hinge', penalty = 'l1', dual = False, max_iter=10000)
 estimators.append(('SVM', model2))
 model3 = RandomForestClassifier(n_estimators = 100)
 estimators.append(('RandomForest', model3))
 ensemble = VotingClassifier(estimators)
 results = model_selection.cross_val_score(ensemble, data, label, cv=kfold)
 print("Ensemble Scores") 
 print(results.mean())


#run the functions for Lung_sig_result_py.txt
path='/Users/'
genes1 = pd.read_csv(path+'Lung_sig_result_py.txt', index_col = 0, sep = "\t")
genes1.drop(genes1.columns[[0,1,2]], axis=1, inplace=True)
data1 = genes1.transpose()
labels1 = np.zeros(120);
labels1[0:60]=1
labels1[60:120]=0
CV(data=data1,label=labels1)
emsemble(data=data1,label=labels1)

#run the functions for LungCancer.Preprocessed.txt
genes2 = pd.read_csv(path+'LungCancer.Preprocessed.txt', index_col = 0, sep = "\t")
data2 = genes2.transpose()
labels2 = np.zeros(120);
labels2[0:60]=1
labels2[60:120]=0
CV(data=data2,label=labels2)
emsemble(data=data2,label=labels2)
