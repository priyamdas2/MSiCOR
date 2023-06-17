# To draw example sequence heatmaps
# install.packages("ComplexHeatmap")
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(DEoptim)
library(optimr)
library(expm)
library(Rfast)
library(purrr)
library(readxl)
library(gplots)
#library(ComplexHeatmap)

setwd("D:/MSiCOR_JCGS/MSiCOR_REAL_DATA/STEP_2_REAL_DATA_result")

# 1 = 121191 = rituximab + ocrelizumab
# 2 = 214582 = glatiramer acetate
# 3 = NA = Interferon-beta
# 4 = 1373478 = dimethyl fumarate
# 5 = 354770 = natalizumab
# 6 = 1012892 = fingolimod
# 7 = 1310520 = teriflunomide
# 8 = 92299 = cyclophosphamide
# 9 = 7005 = mitoxantrone
# 10 = 117055 = alemtuzumab

MAT <- as.matrix(read.table("MS_REAL_DATA_EXAMPLE_SEQ_num_clus_3_num_1_100_by_10.csv", header = FALSE, sep=","))

cluster_membership <- read.table("Cluster_membership_ordered.csv", header = FALSE, sep=",")

length(which(cluster_membership[,2] == 1))
length(which(cluster_membership[,2] == 2))
length(which(cluster_membership[,2] == 3))

length(which(cluster_membership[,2] == 1))/dim(cluster_membership)[1]
length(which(cluster_membership[,2] == 2))/dim(cluster_membership)[1]
length(which(cluster_membership[,2] == 3))/dim(cluster_membership)[1]


par(mar = c(0, 0, 2, 2))
heatmap.2(MAT, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',offsetRow=3,key=FALSE,
          col= c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"),
          labels=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","NA"), 
          xlab = "CODES", ylab = "Patients",labRow = FALSE, main = "K=3, 1st CLUSTER, 100 patients")


legend("left", title = "Drug",legend=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","blank"), 
       fill=c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"), cex=1.4, box.lty=1)





MAT <- as.matrix(read.table("MS_REAL_DATA_EXAMPLE_SEQ_num_clus_3_num_2_100_by_10.csv", header = FALSE, sep=","))


par(mar = c(0, 0, 2, 2))
heatmap.2(MAT, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',offsetRow=3,key=FALSE,
          col= c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"),
          labels=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","NA"), 
          xlab = "CODES", ylab = "Patients",labRow = FALSE, main = "K=3, 2nd CLUSTER, 100 patients")


legend("left", title = "Drug",legend=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","blank"), 
       fill=c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"), cex=1.4, box.lty=1)




MAT <- as.matrix(read.table("MS_REAL_DATA_EXAMPLE_SEQ_num_clus_3_num_3_100_by_10.csv", header = FALSE, sep=","))


par(mar = c(0, 0, 2, 2))
heatmap.2(MAT, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',offsetRow=3,key=FALSE,
          col= c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"),
          labels=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","NA"), 
          xlab = "CODES", ylab = "Patients",labRow = FALSE, main = "K=3, 3rd CLUSTER, 100 patients")


legend("left", title = "Drug",legend=c("ritu + ocre","glat. ac.","Inf-beta","dim. fum.","natal.","fing.", "terifmd.", "cyc-phosp.","mitox.","alemtu.","blank"), 
       fill=c("red","brown","yellow","green","blue","blueviolet","black","aquamarine","azure4","darkorange", "white"), cex=1.4, box.lty=1)

