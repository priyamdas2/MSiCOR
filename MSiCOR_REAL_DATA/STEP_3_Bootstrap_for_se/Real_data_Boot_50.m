clear all

filename = ['ZZZ_Theta_initial_value.csv'];
theta_start_temp = csvread(filename);

filename = ['ZZZ_Beta_initial_value.csv'];
beta_start_temp = csvread(filename);

filename = ['Covariate_values_ordered.csv'];
X = csvread(filename);

M = size(theta_start_temp,2);
K = size(theta_start_temp,1)/(M+1);

filename = ['MSiCOR_DATA_after_2006_10_codes_have_covariates.csv'];
MAT = csvread(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Changing serial numbers to 1:822 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unique_patient_nums = unique(MAT(:,1));
MAT2 = MAT;
for ii = 1:size(MAT,1)
    if(ii == 1)
        new_serial_num = 1;
        MAT2(ii,1) = new_serial_num; 
    else
        if(MAT(ii,1) == MAT(ii-1,1))
            MAT2(ii,1) = new_serial_num;
        else
            new_serial_num = new_serial_num + 1;
            MAT2(ii,1) = new_serial_num;
        end
    end
end
        



p = size(beta_start_temp,2);
N = length(unique(MAT2(:,1)));
n = (M+1)*K;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generating Bootsrap sample from MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_boots = 50;

ALL_BETA_EST = zeros(K*num_boots, 4);
ALL_EST_proportions = zeros(num_boots,3);

filename = ['Real_ALL_BETA_EST.csv'];
csvwrite(filename,ALL_BETA_EST);
filename = ['Real_ALL_EST_proportions.csv'];
csvwrite(filename,ALL_EST_proportions);


for boot_no = 1:num_boots
    %%%% Changing beta_start_temp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_start = zeros(K,p);
    
    for ii = 2:K
        for jj = 1:p
            beta_start(ii,jj) = beta_start_temp(ii,jj) + normrnd(0,1);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Changing theta_start_temp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_start = zeros((M+1)*K,M);
    
    for ii = 1:((M+1)*K)
        for jj = 1:M
            theta_start(ii,jj) = theta_start_temp(ii,jj) + abs(normrnd(0,0.01));
        end
        theta_start(ii,:) = theta_start(ii,:)/sum(theta_start(ii,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    tic;
    random_sample = randsample(N,N, true);
    arranged_random_sample = sort(random_sample);
    
    for i = 1:N
        row_nums = (arranged_random_sample(i) == MAT2(:,1));
        num_rows = sum(row_nums);
        if(i == 1)
            NEW_MAT = zeros(num_rows,2);
            NEW_MAT(:,1) = i;
            NEW_MAT(:,2) = MAT2(row_nums,2);
        else
            sub = zeros(num_rows,2);
            sub(:,1) = i;
            sub(:,2) = MAT2(row_nums,2);
            NEW_MAT = [NEW_MAT; sub];
        end
    end
    
    NEW_X = X(arranged_random_sample,:);
    
    sequences_temp = cell(1,N);
    patient_num = 1;
    sequences_temp{1} = NEW_MAT(1,2);
    seq_lengths = zeros(N,1);
    for i = 2:size(NEW_MAT,1)
        if(NEW_MAT(i,1) == NEW_MAT(i-1,1))
            sequences_temp{patient_num} = [sequences_temp{patient_num} NEW_MAT(i,2)];
        else
            patient_num = patient_num + 1;
            sequences_temp{patient_num} = NEW_MAT(i,2);
        end
    end
    sequences = sequences_temp;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIMIZATION INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = theta_start;
    beta_vec = beta_2_beta_vec(beta_start);
    
    fun = @(beta_vec) value_at_BETA_vec(beta_vec,NEW_X,theta,sequences);
    
    value_at_beta_vec(beta_vec,NEW_X,theta,sequences)
    
    tic;
    
    no_iter = 5000;
    lh_array = zeros(no_iter,1);
    tol_fun = 10^(1);                 % minimum improvement to keep epsilon same after iteration
    tol_fun_2 = 10^(1);               % minimum improvement to start another run
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
            
            current_lh_1 = value_at_beta_vec(beta_vec,NEW_X,theta,sequences);
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
                        total_lh_beta(paral_no) = value_at_beta_vec(possibility_pos_1,NEW_X,theta,sequences);
                        matrix_update_at_h_beta(paral_no,:) = possibility_pos_1;
                    else
                        % neg
                        negative_change_loc_1 = paral_no/2;
                        possibility_neg_1 = beta_vec;
                        possibility_neg_1(negative_change_loc_1) = beta_vec(negative_change_loc_1) - epsilon;
                        total_lh_beta(paral_no) = value_at_beta_vec(possibility_neg_1,NEW_X,theta,sequences);
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
                [boot_no, ii,i,j,log(current_lh_1)/log(10), epsilon]
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
                            total_lh(parallel_no) = value_at_beta_vec(beta_vec,NEW_X,proxy_coef_theta_pos,sequences);
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
                            total_lh(parallel_no) = value_at_beta_vec(beta_vec,NEW_X,proxy_coef_theta_neg,sequences);
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
            array_of_values(i) = value_at_beta_vec(beta_vec,NEW_X,theta,sequences);
            
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
    
    beta_start;
    beta_end = beta_vec_2_beta(beta_vec,K);
    
    membership_prob_EST = zeros(N,K);
    
    for i = 1:N
        membership_factors = zeros(1,K);
        for k = 1:K
            membership_factors(k) = exp(NEW_X(i,:)*beta_end(k,:)');
        end
        membership_prob_EST(i,:) = membership_factors/sum(membership_factors);
        [MM, II] = max(membership_prob_EST(i,:));
        %cluster_index_EST(i) = II;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%% Finding Predicted Clusters USING SEQUENCES ONLY %%%%%%%%%%%%%%%%
    Likelihoods_only_seqs = zeros(N,K);
    max_index = zeros(N,1);
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
    
    ESTIMATED_membership = max_index;
    proportions_EST_from_cov_and_seq = zeros(K,1);
    for k = 1:K
        proportions_EST_from_cov_and_seq(k) = sum(max_index == k)/N;
    end
    proportions_EST_from_cov_and_seq;
    

    filename = ['Real_ALL_BETA_EST.csv'];
    ALL_BETA_EST = csvread(filename);
    ALL_BETA_EST((3*(boot_no-1)+1):(3*(boot_no-1)+3),1:p) = beta_end;
    
    filename = ['Real_ALL_EST_proportions.csv'];
    ALL_EST_proportions = csvread(filename);
    ALL_EST_proportions(boot_no, 1:K) = proportions_EST_from_cov_and_seq';
    
    
    filename = ['Real_ALL_BETA_EST.csv'];
    csvwrite(filename,ALL_BETA_EST);
    filename = ['Real_ALL_EST_proportions.csv'];
    csvwrite(filename,ALL_EST_proportions);
    
    if(boot_no == 1)
        time_taken = toc;
        filename = ['Real_ZZZ_TOC.csv']; % time to compute 1 bootstrap
        csvwrite(filename,time_taken);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 6;
filename = ['Real_ALL_BETA_EST.csv'];
ALL_BETA_EST = csvread(filename);
se_ALL_BETA_EST = zeros(3,p);
for jj = 1:3
    for kk = 1:p
        for i = 1:num_boots
            mat_corresponding = ALL_BETA_EST((3*(i-1)+1):(3*(i-1)+3),1:p);
            if(i == 1)
                elements = mat_corresponding(jj,kk);
            else
                elements = [elements,mat_corresponding(jj,kk)];
            end
        end
        se_ALL_BETA_EST(jj,kk) = std(elements)/sqrt(num_boots);
    end
end
filename = ['se_Real_ALL_BETA_EST.csv'];
csvwrite(filename,se_ALL_BETA_EST);

filename = ['Real_ALL_EST_proportions.csv'];
ALL_EST_proportions = csvread(filename);
se_ALL_EST_proportions = zeros(1,3);
for jj = 1:3
    for i = 1:num_boots
            row_corresponding = ALL_EST_proportions(i,:);
            if(i == 1)
                elements = row_corresponding(jj);
            else
                elements = [elements,row_corresponding(jj)];
            end
        end
        se_ALL_EST_proportions(jj) = std(elements)/sqrt(num_boots);
end
filename = ['se_Real_ALL_EST_proportions.csv'];
csvwrite(filename,se_ALL_EST_proportions);