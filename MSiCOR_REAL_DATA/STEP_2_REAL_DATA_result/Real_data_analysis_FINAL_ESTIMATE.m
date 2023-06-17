%clear all

filename = ['ZZZ_Theta_initial_value.csv'];
theta_start = csvread(filename);

filename = ['ZZZ_Beta_initial_value.csv'];
beta_start = csvread(filename);

filename = ['MSiCOR_DATA_after_2006_10_codes_have_covariates.csv'];
MAT = csvread(filename);

unique_patient_ids = unique(MAT(:,1));

filename = ['Covariate_values_ordered.csv'];
X = csvread(filename);

M = size(theta_start,2);
K = size(theta_start,1)/(M+1);
p = size(beta_start,2);
N = length(unique(MAT(:,1)));
n = (M+1)*K;


sequences_temp = cell(1,N);
patient_num = 1;
sequences_temp{1} = MAT(1,2);
seq_lengths = zeros(N,1);
for i = 2:size(MAT,1)
    if(MAT(i,1) == MAT(i-1,1))
        sequences_temp{patient_num} = [sequences_temp{patient_num} MAT(i,2)];
    else
        patient_num = patient_num + 1;
        sequences_temp{patient_num} = MAT(i,2);
    end
end
sequences = sequences_temp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZATION INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = theta_start;
beta_vec = beta_2_beta_vec(beta_start);

fun = @(beta_vec) value_at_BETA_vec(beta_vec,X,theta,sequences);

value_at_beta_vec(beta_vec,X,theta,sequences)

tic;

no_iter = 5000;
lh_array = zeros(no_iter,1);
tol_fun = 10^(-1);                 % minimum improvement to keep epsilon same after iteration
tol_fun_2 = 10^(-1);               % minimum improvement to start another run
parameter_cut_off = 10^(-20);      % hard-thresholding of each coordinate value
epsilon_cut_off = 10^(-3);         % smallest step-size
epsilon_decreasing_factor_1 = 2;   % step-decay rate (default is 2)
epsilon_decreasing_factor_2 = 2;
max_runs = 10;                     % Number of runs





like_corresponding_run = zeros(max_runs,1);

P = (K-1)*p; % Total number of parameters 

for ii = 1:max_runs
    old_theta = theta;
    epsilon = 1;
    if(ii == 1)
        epsilon_decreasing_factor = epsilon_decreasing_factor_1;
    else
        epsilon_decreasing_factor = epsilon_decreasing_factor_2;
    end
    
    array_of_values = zeros(no_iter,1);
    
    for i=1:no_iter
        count_1 = 0;
        count_2 = 0;
        check_each_best_prop = zeros(1,K);
        check_each_best = zeros(n,M);
        corresponding_likelihood = zeros(n,1);
        
        current_lh_1 = value_at_beta_vec(beta_vec,X,theta,sequences);
        if(i == 1)
            num = current_lh_1;
        end
        
       %%%%%%%%%%%%%%%%%%%% UPDATING  beta_vec %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for f = 1:P
            total_lh_beta = zeros(2*P,1);
            matrix_update_at_h_beta = zeros(2*P, P);
            
            for paral_no = 1:(2*P)   % CAN BE RUN IN PARALLEL
                % pos
                if(mod(paral_no,2) == 1)
                    positive_change_loc_1 = (paral_no + 1)/2;
                    possibility_pos_1 = beta_vec;
                    possibility_pos_1(positive_change_loc_1) = beta_vec(positive_change_loc_1) + epsilon;
                    total_lh_beta(paral_no) = value_at_beta_vec(possibility_pos_1,X,theta,sequences);
                    matrix_update_at_h_beta(paral_no,:) = possibility_pos_1;
                else
                    % neg
                    negative_change_loc_1 = paral_no/2;
                    possibility_neg_1 = beta_vec;
                    possibility_neg_1(negative_change_loc_1) = beta_vec(negative_change_loc_1) - epsilon;
                    total_lh_beta(paral_no) = value_at_beta_vec(possibility_neg_1,X,theta,sequences);
                    matrix_update_at_h_beta(paral_no,:) = possibility_neg_1; 
                end
            end
        end
        
        [num_beta, idx_beta] = min(total_lh_beta);
        
        
       
       %%%%%%%%%%%%%%%%%%%% UPDATING TRANSITION MATRICES %%%%%%%%%%%%%%%%%%
        for j = 1:n
            if(min(ge(theta(j,:),0)) == 0)
                stop('error')
            end
            if(sum(theta(j,:))<0.99 || sum(theta(j,:))>1.01)
                stop('error')
            end
            [ii,i,j,log(current_lh_1)/log(10), epsilon]
            %beta_vec'
            %current_lh_1
            total_lh = zeros(2*M,1);
            matrix_update_at_h = zeros(2*M,M);
            
            for parallel_no = 1:(2*M)   % CAN BE RUN IN PARALLEL
                % pos
                if(mod(parallel_no,2) == 1)
                    positive_change_loc = (parallel_no + 1)/2;
                    possibility_pos = theta(j,:);
                    temp_possibility_pos = theta(j,:);
                    temp_possibility_pos(positive_change_loc) = 0; % To find all significant positions except the positive_change_loc-th
                    significant_positions = find(gt(temp_possibility_pos, parameter_cut_off*ones(1,M)));
                    if(isempty(significant_positions) == 1)
                        possibility_pos = theta(j,:);
                    else
                        possibility_pos(positive_change_loc) = theta(j,positive_change_loc) + epsilon;
                        possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon/length(significant_positions);
                        epsilon_temp_pos = epsilon;
                        
                        if(min(ge(possibility_pos,0)) == 0 && epsilon_temp_pos > epsilon_cut_off)
                            epsilon_temp_pos = epsilon_temp_pos/epsilon_decreasing_factor;
                            possibility_pos = theta(j,:);
                            possibility_pos(positive_change_loc) = theta(j,positive_change_loc) + epsilon_temp_pos;
                            possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon_temp_pos/length(significant_positions);
                        end
                    end
                    if(sum(possibility_pos) <0.99 || sum(possibility_pos)>1.01)  %% Checking if sum is 1
                        stop('error')
                    end
                    
                    if(min(ge(possibility_pos,0)) == 0 || isequal(possibility_pos,theta(j,:)) == 1)
                        possibility_pos = theta(j,:);
                        total_lh(parallel_no) = current_lh_1;
                    else
                        proxy_coef_theta_pos = theta;
                        proxy_coef_theta_pos(j,:) = possibility_pos;
                        total_lh(parallel_no) = value_at_beta_vec(beta_vec,X,proxy_coef_theta_pos,sequences);
                    end
                    matrix_update_at_h(parallel_no,:) = possibility_pos;
                else
                    % neg
                    negative_change_loc = parallel_no/2;
                    possibility_neg = theta(j,:);
                    temp_possibility_neg = theta(j,:);
                    temp_possibility_neg(negative_change_loc) = 0;
                    significant_positions = find(gt(temp_possibility_neg, parameter_cut_off*ones(1,M)));
                    if(isempty(significant_positions) == 1)
                        possibility_neg = theta(j,:);
                    else
                        possibility_neg(negative_change_loc) = theta(j,negative_change_loc) - epsilon;
                        possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon/length(significant_positions);
                        epsilon_temp_neg = epsilon;
                        
                        if(min(ge(possibility_neg,0)) == 0 && epsilon_temp_neg > epsilon_cut_off)
                            epsilon_temp_neg = epsilon_temp_neg/epsilon_decreasing_factor;
                            possibility_neg = theta(j,:);
                            possibility_neg(negative_change_loc) = theta(j,negative_change_loc) - epsilon_temp_neg;
                            possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon_temp_neg/length(significant_positions);
                        end
                    end
                    if(sum(possibility_neg) <0.99 || sum(possibility_neg)>1.01)  %% Checking if sum is 1
                        stop('error')
                    end
                    if(min(ge(possibility_neg,0)) == 0 || isequal(possibility_neg,theta(j,:)) == 1)
                        possibility_neg = theta(j,:);
                        total_lh(parallel_no) = current_lh_1;
                    else
                        proxy_coef_theta_neg = theta;
                        proxy_coef_theta_neg(j,:) = possibility_neg;
                        total_lh(parallel_no) = value_at_beta_vec(beta_vec,X,proxy_coef_theta_neg,sequences);
                    end
                    matrix_update_at_h(parallel_no,:) = possibility_neg;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [Min_here,I] = min(total_lh);
            check_each_best(j,:) = theta(j,:);
            corresponding_likelihood(j) = current_lh_1;
            
            if(Min_here < current_lh_1)
                count_1 = count_1+1;
                check_each_best(j,:) = matrix_update_at_h(I,:);
                corresponding_likelihood(j) = Min_here;
            end
            
        end
        
        [num, idx] = min(corresponding_likelihood(:));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simultaneous Best move selection

        if(min(num_beta,num) < current_lh_1)
            if(num < num_beta)
                parameter_1 = check_each_best(idx,:);
                sparsity_positions = lt(parameter_1,parameter_cut_off*ones(1,M));
                garbage = sum(parameter_1(sparsity_positions));
                if(garbage > 0)
                    parameter_1(sparsity_positions) = 0;
                    rich_positions = ge(parameter_1,parameter_cut_off*ones(1,M));
                    parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
                end
                theta(idx,:) = parameter_1;
            else
                beta_vec = matrix_update_at_h_beta(idx_beta,:)';
            end
        end
        array_of_values(i) = value_at_beta_vec(beta_vec,X,theta,sequences);
        
        % Epsilon change
        if(i > 1)
            if(abs(array_of_values(i) - array_of_values(i-1)) < tol_fun)
                if(epsilon > epsilon_decreasing_factor*epsilon_cut_off)
                    epsilon = epsilon/epsilon_decreasing_factor;
                else
                    break
                end
            end
        end
        
    end
    
    like_corresponding_run(ii)= current_lh_1;
    if(ii >= 2)
        if(abs(like_corresponding_run(ii) - like_corresponding_run(ii-1)) < tol_fun_2)
            break
        end
    end
end




beta_start
beta_end = beta_vec_2_beta(beta_vec,K)

filename = ['MS_REAL_DATA_EXAMPLE_SEQ_BETA_num_clus_',num2str(K),'.csv'];
    csvwrite(filename,round(beta_end,2))

%%%%%%%%%% Finding Predicted Clusters USING COVARIATES ONLY %%%%%%%%%%%%%%%
% Likelihoods_only_covs = zeros(N,K);
% cluster_index_PREDICTED = zeros(N,1);
% M = length(unique(MAT(:,2)));
% for i = 1:N
%     for k = 1:K    % K = num of clusters
%         theta_temp =  theta(((M+1)*(k-1)+1):(k*(M+1)),:);
%         initial_dist_temp = theta_temp(1,:);
%         matrice_temp = theta_temp(2:(M+1),:);
%         m = M * (sequences{i}(2:end) - 1) + sequences{i}(1:(end - 1));
%         temp = initial_dist_temp(sequences{i}(1)) * prod(matrice_temp(m));
%         Likelihoods_only_covs(i,k) = temp;
%     end
%     [Max,I] = max(Likelihoods_only_covs(i,:));
%     cluster_index_PREDICTED(i) = I;
% end


membership_prob_EST = zeros(N,K);
%cluster_index_EST = zeros(N,1);

for i = 1:N
    membership_factors = zeros(1,K);
    for k = 1:K
        membership_factors(k) = exp(X(i,:)*beta_end(k,:)');
    end
    membership_prob_EST(i,:) = membership_factors/sum(membership_factors);
    [MM, II] = max(membership_prob_EST(i,:));
    %cluster_index_EST(i) = II;
end
% proportions_EST = zeros(K,1);
% for k = 1:K
%     proportions_EST(k) = sum(cluster_index_EST == k)/N;
% end
% proportions_EST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Finding Predicted Clusters USING SEQUENCES ONLY %%%%%%%%%%%%%%%%
Likelihoods_only_seqs = zeros(N,K);
max_index = zeros(N,1);
M = length(unique(MAT(:,2)));
for i = 1:N
    for k = 1:K    % K = num of clusters
        theta_temp =  theta(((M+1)*(k-1)+1):(k*(M+1)),:);
        initial_dist_temp = theta_temp(1,:);
        matrice_temp = theta_temp(2:(M+1),:);
        m = M * (sequences{i}(2:end) - 1) + sequences{i}(1:(end - 1));
        temp = initial_dist_temp(sequences{i}(1)) * prod(matrice_temp(m));
        Likelihoods_only_seqs(i,k) = temp;
    end
    Likelihood_covs_and_seqs = membership_prob_EST(i,:).*Likelihoods_only_seqs(i,:);
    [Max,I] = max(Likelihood_covs_and_seqs);
    max_index(i) = I;
end


MS_patient_membership = [unique_patient_ids,max_index];

filename = ['Cluster_membership_ordered.csv'];
csvwrite(filename,MS_patient_membership)

proportions_EST_from_cov_and_seq = zeros(K,1);
for k = 1:K
    proportions_EST_from_cov_and_seq(k) = sum(max_index == k)/N;
end
proportions_EST_from_cov_and_seq

filename = ['MS_REAL_DATA_EXAMPLE_SEQ_PROPORTIONS_num_clus_',num2str(K),'.csv'];
    csvwrite(filename,proportions_EST_from_cov_and_seq)
%%%%% For plots of example patients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_patient_each_clus = 100;
min_length = 10;

Heat_matrix_cluster = cell(K,1);

Cluster_index = cell(K,1);
for k = 1:K
    Cluster_index{k} = -1;
    for i = 1:N
        if(max_index(i) == k) % Or use cluster_index_EST
            Cluster_index{k} = [Cluster_index{k},i];
        end
    end
    Cluster_index{k}(1) = [];
end

for ii = 1:N
    sequences{ii} = [sequences{ii},repelem(11,min_length)];
end

for k = 1:K
    Heat_matrix_cluster{k} = zeros(Num_patient_each_clus,min_length);
    count = 0;
    j = 1;
    while(count < Num_patient_each_clus)
        index_temp = Cluster_index{k}(j);
        if(length(sequences{index_temp}) >= min_length)
            count = count+1;
            Heat_matrix_cluster{k}(count,:) = sequences{index_temp}(1:min_length);
        end
        j = j+1;
    end
end


for k = 1:K
    filename = ['MS_REAL_DATA_EXAMPLE_SEQ_num_clus_',num2str(K),'_num_',num2str(k),'_',...
        num2str(Num_patient_each_clus),'_by_',num2str(min_length),'.csv'];
    csvwrite(filename,Heat_matrix_cluster{k})
end


proportions_EST_from_cov_and_seq

