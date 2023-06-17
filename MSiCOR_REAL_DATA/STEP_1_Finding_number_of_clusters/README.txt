
 This section is to find the optimal number of clusters. Data is fit for
 different possible number of clusters, eventually finding BIC for each possible
 number of clusters.

1) Run 'Real_data_analysis_INITIAL_ESTIMATE.m' to find initial estimate of the
  parameters (set desired value of K, i.e., number of clusters).

2) Run 'Real_data_analysis_FINAL_ESTIMATE.m' to find final estimate of the
  parameters along with corresponding BIC values ((set desired value of K,
  i.e., number of clusters))

3) BIC values for different models can be found in the following generated files
  (i) For K = 3 : 'MS_REAL_DATA_AIC_BIC_corrected_num_clus_3.csv' 
 (ii) For K = 4 : 'MS_REAL_DATA_AIC_BIC_corrected_num_clus_4.csv' 
 and so on . . .

4) These BIC values are noted in Table 6 of the MSiCOR paper/draft.

