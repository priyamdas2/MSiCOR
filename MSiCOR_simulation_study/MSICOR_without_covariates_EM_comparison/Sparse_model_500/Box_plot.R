rm(list = ls())
#install.packages("permutations")
library(permutations)
library(seqHMM)
setwd("D:/MSiCOR_JCGS/MSiCOR_simulation_study/MSICOR_without_covariates_EM_comparison/Sparse_model_500")


seqHMM_TVD_TPR <-
  read.table(
    paste0("AAA_seqHMM_TVD_TPR.csv"),
    header = FALSE,
    sep = ","
  )

MSiCOR_TVD_TPR <-
  read.table(
    paste0("AAA_MSICOR_all_TVD_TPR.csv"),
    header = FALSE,
    sep = ","
  )

# boxplot(seqHMM_TVD_TPR[,1], MSiCOR_TVD_TPR[,1],
#         main = "Multiple boxplots for comparision",
#         at = c(1,2),
#         names = c("EM", "MSiCOR"),
#         #col = c("orange","red"),
#         border = "brown",
#         horizontal = FALSE,
#         notch = TRUE
# )

boxplot(log(MSiCOR_TVD_TPR[,1]),log(seqHMM_TVD_TPR[,1]), names=c("MSiCOR", "EM"), #main="Car Milage Data", 
        #xlab="Methods",
        ylab ="log (TVD)")

max(seqHMM_TVD_TPR[,1])
1 - min(seqHMM_TVD_TPR[,2])

max(MSiCOR_TVD_TPR[,1])
1 - min(MSiCOR_TVD_TPR[,2])

# abline(h=possibles_1[indices[1]],lwd=1, lty=2)
# abline(h=possibles_2[indices[2]],lwd=1, lty=2)