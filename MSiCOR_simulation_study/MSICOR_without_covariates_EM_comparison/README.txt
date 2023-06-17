In this section we compare MSiCOR with EM algorithm for clustering
MMM without covariates. Please follow the following instructions
to get the outputs for the scenarios n=500 and 1000.

Here we describe how to generate the outputs for n = 500 case.
Similar steps can be followed to get the corresponding outputs
for n = 1000 scenario.

1) Go to 'Sparse_model_500' folder.

2) Run 'DATA_Generate_500_sparse.m' to generate data for 50 
 replications.

3) Run 'MSICOR.m' to generate summarized outputs for MSiCOR
 method.

4) Run 'seqHMM_EM_all_reps.R' to generate summarized outputs
 for EM algorithm, using 'seqHMM' R packages.

5) 'AAA_MSICOR_output.csv' gives the output values of TVD 
 (se of TVD), MR (se of MR) for MSiCOR, which is noted
 down in Table 5 of MSiCOR paper/draft.

6) 'AAA_seqHMM_output.csv' gives the output values of TVD 
 (se of TVD), MR (se of MR) for EM algorithm, which is noted
 down in Table 5 of MSiCOR paper/draft.

7) Run 'Box_plot.R' to generate the box-plot shown at Figure 2
 of the MSiCOR paper/draft. [This boxplot function is available
 only for n = 500 scenario].
 