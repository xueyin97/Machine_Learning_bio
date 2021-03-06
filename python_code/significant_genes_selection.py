import pandas as pd
from scipy import stats
path='/Users/' #path to the data file
predata=pd.read_csv(path+'LungCancer.Preprocessed.txt',sep='\t',index_col='ProbeID')
Cancer=predata.iloc[:,:60]
Control=predata.iloc[:,60:]
t_result=stats.ttest_ind(Cancer,Control,axis=1)
t_result.pvalue

from statsmodels.stats import multitest
p_adj=multitest.multipletests(t_result.pvalue, method='fdr_bh')[1]
result = predata
result.insert(0, 'BH_p-value', p_adj)
result.insert(0,'p-value', t_result.pvalue)
sig_result=result[result['p-value']<1e-10]
sig_result.to_csv(path+'hw4_result.txt', sep='\t')
