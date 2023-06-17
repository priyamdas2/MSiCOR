clear all
K = 3;
num_simulations = 20;
shift = 10; % put random values >= 0, <= 80

TPR_MAT = zeros(num_simulations,1);
PROP_MAT = zeros(3,1,num_simulations);
TRUE_PROP_MAT = zeros(3,1,num_simulations);
BETA_MAT = zeros(3,4,num_simulations);

for i = 1:num_simulations
    
    
    filename = ['ZZZ_TPR_clus_',num2str(K),'_seed_',num2str(i+shift),'.csv'];
    TPR = csvread(filename);
    TPR_MAT(i,1) = TPR;
    
    filename = ['ZZZ_EST_proportions_',num2str(K),'_seed_',num2str(i+shift),'.csv'];
    PROP = csvread(filename);
    PROP_MAT(:,:,i) = PROP;
    
    filename = ['ZZZ_TRUE_proportions_',num2str(K),'_seed_',num2str(i+shift),'.csv'];
    TRUE_PROP = csvread(filename);
    TRUE_PROP_MAT(:,:,i) = TRUE_PROP;
    
    filename = ['ZZZ_BETA_EST_num_clus_',num2str(K),'_seed_',num2str(i+shift),'.csv'];
    BETA = csvread(filename);
    BETA_MAT(:,:,i) = BETA;
    
end

TPR_mean_se = [mean(TPR_MAT),std(TPR_MAT)/sqrt(num_simulations)];

PROP_MAT_mean_se = zeros(K,2);
TRUE_PROP_MAT_mean_se = zeros(K,2);
BETA_MAT_mean_se = zeros(K,8);

for ii = 1:K
    elements = PROP_MAT(ii,1,1:num_simulations);
    PROP_MAT_mean_se(ii,1) = mean(elements);
    PROP_MAT_mean_se(ii,2) = std(elements)/sqrt(num_simulations);
    
    elements_2 = TRUE_PROP_MAT(ii,1,1:num_simulations);
    TRUE_PROP_MAT_mean_se(ii,1) = mean(elements_2);
    TRUE_PROP_MAT_mean_se(ii,2) = std(elements_2)/sqrt(num_simulations);
    
    for jj = 1:4
        elements = BETA_MAT(ii,jj,1:num_simulations);
        BETA_MAT_mean_se(ii,jj) = mean(elements);
        BETA_MAT_mean_se(ii,jj+4) = std(elements)/sqrt(num_simulations);
    end
end
 
 filename = ['ZZZ_AAA_TPR_num_clus_',num2str(K),'_numiters_',num2str(num_simulations),'.csv'];
 csvwrite(filename,TPR_mean_se);
 
 filename = ['ZZZ_AAA_BETA_EST_num_clus_',num2str(K),'_numiters_',num2str(num_simulations),'.csv'];
 csvwrite(filename,BETA_MAT_mean_se);
 
 filename = ['ZZZ_AAA_EST_proportions_',num2str(K),'_numiters_',num2str(num_simulations),'.csv'];
 csvwrite(filename, PROP_MAT_mean_se);

 filename = ['ZZZ_AAA_TRUE_proportions_',num2str(K),'_numiters_',num2str(num_simulations),'.csv'];
 csvwrite(filename, TRUE_PROP_MAT_mean_se);
