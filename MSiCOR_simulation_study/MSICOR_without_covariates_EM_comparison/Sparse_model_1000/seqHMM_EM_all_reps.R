rm(list = ls())
#install.packages("permutations")
library(permutations)
library(seqHMM)
setwd("D:/MSiCOR_JCGS/MSiCOR_simulation_study/MSICOR_without_covariates_EM_comparison/Sparse_model_1000")
#set.seed(4)

K <- 3  # Number of clusters
M <- 10 # number of states
ZZZ_Theta_TRUE <-
  read.table(
    file = paste("AAA_Theta_TRUE.csv", sep = ""),
    header = FALSE,
    sep = ","
  )
num_rep <- 50

TVD_EM_array <- array(NA, num_rep)
TPR_array <- array(NA, num_rep)

#######################################################################
######## Computation Begins ###########################################
#######################################################################

for (rept in 1:num_rep)
{
  #set.seed(rept+2)
  print(rept)
  DATA_sequences <-
    read.table(
      file = paste0("ZZZ_DATA_sequences_rep_", rept, ".csv", sep = ""),
      header = FALSE,
      sep = ","
    )
  ZZZ_TRUE_membership <-
    read.table(
      file = paste("ZZZ_TRUE_membership_rep_", rept, ".csv", sep = ""),
      header = FALSE,
      sep = ","
    )
  
  unique_patient_nums <- unique(DATA_sequences[, 1])
  num_unique_patients <- length(unique_patient_nums)
  
  seq_mat <- matrix(NA, num_unique_patients, 30)
  med_seq_lengths <- rep(NA, num_unique_patients)
  
  for (i in 1:num_unique_patients)
  {
    med_seq <- DATA_sequences[which(DATA_sequences[, 1] == i), 2]
    seq_mat[i, 1:length(med_seq)] <- med_seq
    med_seq_lengths[i] <- length(med_seq)
  }
  
  max_med_seq_length <- max(med_seq_lengths)
  seq_mat_cut <- seq_mat[, 1:max_med_seq_length]
  t1 <- seqdef(seq_mat_cut)
  
  
  
  #######################################################################
  ######## RUNNING EM Algorithm #########################################
  #######################################################################
  
  transition_probs <- list()
  initial_probs <- list()
  
  for (i in 1:K)
  {
    random_mat <- matrix(NA, M, M)
    random_init_prob <- array(NA, M)
    for (jjjj in 1:M)
    {
      random_init_prob[jjjj] <- runif(1)
      for (kkkk in 1:M)
      {
        random_mat[jjjj, kkkk] <- runif(1)
      }
      random_mat[jjjj,] <-
        random_mat[jjjj,] / sum(random_mat[jjjj,])
    }
    random_init_prob <- random_init_prob / sum(random_init_prob)
    
    
    
    transition_probs[[i]] <- random_mat
    initial_probs[[i]] <- as.vector(random_init_prob)
  }
  
  mmm <-
    build_mmm(t1,
              K,
              transition_probs = transition_probs,
              initial_probs = initial_probs)
  
  fit_EM <- fit_model(mmm, global_step = FALSE, local_step = FALSE)
  
  EST_mat <-
    rbind(
      fit_EM[[1]][4]$initial_probs$`Cluster 1`,
      fit_EM[[1]][2]$transition_probs$`Cluster 1`,
      fit_EM[[1]][4]$initial_probs$`Cluster 2`,
      fit_EM[[1]][2]$transition_probs$`Cluster 2`,
      fit_EM[[1]][4]$initial_probs$`Cluster 3`,
      fit_EM[[1]][2]$transition_probs$`Cluster 3`
    )
  PROP_mat <- summary(fit_EM$model)$prior_cluster_probabilities[1,]
  
  
  #######################################################################
  ######## Calculating Total Variation Distance #########################
  #######################################################################
  
  EST_mat_list <- list()
  TRUE_mat_list <- list()
  
  for (k in 1:K)
  {
    EST_mat_list[[k]] <-
      EST_mat[((k - 1) * (M + 1) + 1):(k * (M + 1)),]
    TRUE_mat_list[[k]] <-
      ZZZ_Theta_TRUE[((k - 1) * (M + 1) + 1):(k * (M + 1)),]
  }
  
  
  perm_sequences <- matrix(NA, 6, 3)
  
  perm_sequences[1,] <- c(1, 2, 3)
  perm_sequences[2,] <- c(1, 3, 2)
  perm_sequences[3,] <- c(2, 1, 3)
  perm_sequences[4,] <- c(2, 3, 1)
  perm_sequences[5,] <- c(3, 1, 2)
  perm_sequences[6,] <- c(3, 2, 1)
  
  TVDs <- rep(NA, 6)
  
  for (i in 1:6)
  {
    TVDs[i] <-
      sum(abs(TRUE_mat_list[[1]] - EST_mat_list[[perm_sequences[i, 1]]])) +
      sum(abs(TRUE_mat_list[[2]] - EST_mat_list[[perm_sequences[i, 2]]])) +
      sum(abs(TRUE_mat_list[[3]] - EST_mat_list[[perm_sequences[i, 3]]]))
  }
  
  right_order <- perm_sequences[which.min(TVDs),]
  TVD_EM <- TVDs[which.min(TVDs)]
  
  
  EST_mat_ORDERED <-
    rbind(EST_mat_list[[right_order[1]]], EST_mat_list[[right_order[2]]],
          EST_mat_list[[right_order[3]]])
  
  PROP_mat_ORDERED <- PROP_mat[perm_sequences[which.min(TVDs),]]
  
  sum(abs(ZZZ_Theta_TRUE - EST_mat_ORDERED))
  
  # 
  # write.table(
  #   EST_mat_ORDERED,
  #   paste0("initial_mat_using_seqHMM_rep_",rept,".csv"),
  #   row.names = FALSE,
  #   sep = ",",
  #   col.names = FALSE
  # )
  # 
  # write.table(
  #   PROP_mat_ORDERED,
  #   paste0("initial_prop_using_seqHMM_rep_",rept,".csv"),
  #   row.names = FALSE,
  #   sep = ",",
  #   col.names = FALSE
  # )
  # 
  
  #######################################################################
  ######## Finding estimated cluster for each med-sequence ##############
  #######################################################################
  
  product <- matrix(NA, num_unique_patients, K)
  
  for (k in 1:K)
  {
    est_mat_this_clus <-
      EST_mat_ORDERED[((k - 1) * (M + 1) + 1):(k * (M + 1)),]
    initial_prob <- est_mat_this_clus[1,]
    transition_mat <- est_mat_this_clus[2:(M + 1),]
    prior_prop <- PROP_mat_ORDERED[k]
    
    for (j in 1:num_unique_patients)
    {
      num_non_zeros <- length(which(is.na(seq_mat_cut[j,]) == FALSE))
      element_position_mat <- matrix(NA, num_non_zeros - 1, 2)
      product[j, k] <- prior_prop * initial_prob[seq_mat_cut[j, 1]]
      for (kk in 1:(num_non_zeros - 1))
      {
        element_position_mat[kk,] <-
          c(seq_mat_cut[j, kk], seq_mat_cut[j, kk + 1])
        product[j, k] <-
          product[j, k] * transition_mat[element_position_mat[kk, 1], element_position_mat[kk, 2]]
      }
    }
  }
  
  max_index <- rep(NA, num_unique_patients)
  
  for (ii in 1:num_unique_patients)
  {
    max_index[ii] <- which.max(product[ii,])
  }
  
  TPR <-
    length(which(ZZZ_TRUE_membership - max_index == 0)) / num_unique_patients
  
  True_props <- rep(NA, K)
  Est_props <- rep(NA, K)
  for (k in 1:K)
  {
    True_props[k] <-
      length(which(ZZZ_TRUE_membership == k)) / num_unique_patients
    Est_props[k] <-
      length(which(max_index == k)) / num_unique_patients
  }
  TVD_EM_array[rept] <- TVD_EM
  TPR_array[rept] <- TPR
  
}


std <- function(x)
  sd(x) / sqrt(length(x))

output <- matrix(NA, 2, 2)
output[1, 1] <- mean(TVD_EM_array)
output[1, 2] <- std(TVD_EM_array)
output[2, 1] <- mean(TPR_array)
output[2, 2] <- std(TPR_array)


TVD_TPR <- cbind(TVD_EM_array,TPR_array)

write.table(
  TVD_TPR,
  paste0("AAA_seqHMM_TVD_TPR.csv"),
  row.names = FALSE,
  sep = ",",
  col.names = FALSE
)

write.table(
  output,
  paste0("AAA_seqHMM_output.csv"),
  row.names = FALSE,
  sep = ",",
  col.names = FALSE
)
