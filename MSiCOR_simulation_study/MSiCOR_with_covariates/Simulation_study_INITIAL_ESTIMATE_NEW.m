clear all
num_simulations = 100;
rng(1)
% Consider K many markov chains, each having M many states, so there will
% be K many MxM transition matrices and one vector of length M (proportion
% vector). Here 'theta' denotes K many (M+1)x M matrices stacked one after
% another. So 'theta' is of ((M+1)K)x K diemnsional. 'theta_prop' denotes the
% proportion vetor of Kx1 diemnsion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10; % number of states
K = 3; % number of clusters
p = 4; % number of variables
N = 1000;
%%%%%% Generating Transition matrix+initial prob vectors %%%%%%%%%%%%%%%%%

Transitions = cell(K,1);
initial_DIST_TRUE = cell(K,1);

dirich_par = ones(1,M);

for k = 1:K
    Transitions{k} = drchrnd(dirich_par, M);
    initial_DIST_TRUE{k} = drchrnd(dirich_par, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Generating covaraites %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta_TRUE = zeros(K,p);
X = ones(N,p);

for i = 2:K
    Beta_TRUE(i,:) = round(normrnd(0,1,1,p),2);
end

filename = ['ZZZ_BETA_TRUE_num_clus_',num2str(K),'.csv'];
csvwrite(filename,round(Beta_TRUE,2))

for i = 1:N
    X(i,2:p) = normrnd(0,1,1,p-1);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Finding Membership %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

membership_prob_TRUE = zeros(N,K);
cluster_index_DETERMINISTIC = zeros(N,1);

for i = 1:N
    membership_factors = zeros(1,K);
    for k = 1:K
        membership_factors(k) = exp(X(i,:)*Beta_TRUE(k,:)');
    end
    membership_prob_TRUE(i,:) = membership_factors/sum(membership_factors);
    %     [MM, II] = max(membership_prob_TRUE(i,:));
    %     cluster_index(i) = II;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Building THETA biscuits: putting THETA vector for optimization %%%

for k = 1:K
    theta_biscuit = [initial_DIST_TRUE{k};Transitions{k}];
    if(k == 1)
        theta_TRUE = theta_biscuit;
    else
        theta_TRUE = [theta_TRUE; theta_biscuit];
    end
end






for zzzz = 1:num_simulations
    rng(zzzz)
    
    
    %%%%%%%%%% DATA GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    chain_length = unidrnd(7,N,1)+5;
    cumsum_chain_length = cumsum(chain_length);
    
    TRUE_membership = zeros(N,1);
    MAT = zeros(sum(chain_length),2);
    
    
    
    
    for i = 1:N
        % generate chain
        temp_seq = -1*ones(chain_length(i),1);
        pd = makedist('Multinomial','probabilities',membership_prob_TRUE(i,:));
        r = random(pd);
        pp = r;
        %     [aa,bb] = max(membership_prob_TRUE(i,:));
        %     pp = bb;
        TRUE_membership(i) = pp;
        initial_value_pppp = mnrnd(1,initial_DIST_TRUE{pp});
        initial_value = initial_value_pppp*[1:M]';
        temp_seq(1) = initial_value;
        
        
        for j = 2:chain_length(i)
            temp_prob_vec = Transitions{pp}(temp_seq(j-1),:);
            qqqq = mnrnd(1,temp_prob_vec);
            qq = qqqq*[1:M]';
            temp_seq(j) = qq;
        end
        
        
        if(i == 1)
            MAT(1:chain_length(i),1) = i;
            MAT(1:chain_length(i),2) = temp_seq;
        else
            MAT((cumsum_chain_length(i-1)+1):cumsum_chain_length(i),1) = i;
            MAT((cumsum_chain_length(i-1)+1):cumsum_chain_length(i),2) = temp_seq;
        end
        
    end
    
    % filename = ['ZZZ_DATA_sequences.csv'];
    % csvwrite(filename,MAT)
    
    
    % filename = ['ZZZ_TRUE_membership.csv'];
    % csvwrite(filename,TRUE_membership)
    
    proportions_TRUE = zeros(K,1);
    
    for k = 1:K
        proportions_TRUE(k) = sum(TRUE_membership == k)/N;
    end
    
    
    
    % filename = ['ZZZ_TRUE_proportions_',num2str(K),'.csv'];
    % csvwrite(filename,proportions_TRUE)
    
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
    objective_fun_identifiable_no_prop_CONTROLLED(theta_TRUE,sequences)
    
    %%%%%%%%%% Finding transition places %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = length(unique(MAT(:,2)));
    Indicator_matrix = zeros(M,M);
    for i = 1:N
        for j = 1:(length(sequences{i})-1)
            aa = sequences{i}(j);
            bb = sequences{i}(j+1);
            Indicator_matrix(aa,bb) = Indicator_matrix(aa,bb) + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    initial_states = zeros(N,1);
    for gg = 1:N % = size(sequences,2)
        initial_states(gg) = sequences{gg}(1);
    end
    
    length(initial_states(initial_states  == 1)) % 4871 in INF based
    tabulate(initial_states)
    tabulate(MAT(:,2))
    
    for j = 1:N
        seq_lengths(j) = size(sequences_temp{j},2);
    end
    
    M = length(unique(MAT(:,2)));
    n = (M+1)*K;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIMIZATION INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %starting_point_prop = proportions_TRUE;
    starting_point = drchrnd(dirich_par, (M+1)*K);
    %starting_point =  theta_TRUE;
    
    %theta_prop = starting_point_prop;
    theta = starting_point;
    
    tic;
    
    no_iter = 5000;
    lh_array = zeros(no_iter,1);
    tol_fun = 10^(-1);                 % minimum improvement to keep epsilon same after iteration
    tol_fun_2 = 10^(-1);               % minimum improvement to start another run
    parameter_cut_off = 10^(-20);      % hard-thresholding of each coordinate value
    epsilon_cut_off = 10^(-3);         % smallest step-size
    epsilon_decreasing_factor_1 = 2;   % step-decay rate (default is 2)
    epsilon_decreasing_factor_2 = 2;
    max_runs = 10;    % Number of runs
    
    like_corresponding_run = zeros(max_runs,1);
    
    
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
            
            current_lh_1 = objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences);
            if(i == 1)
                num = current_lh_1;
            end
            
            %%%%%%%%%%%%%%%%%%%% UPDATING TRANSITION MATRICES %%%%%%%%%%%%%%%%%%
            for j = 1:n
                if(min(ge(theta(j,:),0)) == 0 || min(ge(theta(j,:),0)) == 0 )
                    stop('error')
                end
                if(sum(theta(j,:))<0.99 || sum(theta(j,:))>1.01)
                    stop('error')
                end
                [ii,i,j,log(current_lh_1)/log(10), epsilon]
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
                            total_lh(parallel_no) = objective_fun_identifiable_no_prop_CONTROLLED(proxy_coef_theta_pos,sequences);
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
                            total_lh(parallel_no) = objective_fun_identifiable_no_prop_CONTROLLED(proxy_coef_theta_neg,sequences);
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
            
            [num,current_lh_1]
            
            if(num < current_lh_1)
                parameter_1 = check_each_best(idx,:);
                sparsity_positions = lt(parameter_1,parameter_cut_off*ones(1,M));
                garbage = sum(parameter_1(sparsity_positions));
                if(garbage > 0)
                    parameter_1(sparsity_positions) = 0;
                    rich_positions = ge(parameter_1,parameter_cut_off*ones(1,M));
                    parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
                end
                theta(idx,:) = parameter_1;
            end
            
            array_of_values(i) = objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences);
            
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
    
    das_algo_times = toc;
    das_algo_values = objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences);
    
    theta_TRUE_CELL = cell(K,1);
    theta_EST_CELL = cell(K,1);
    theta_START_CELL = cell(K,1);
    
    for k = 1:K
        theta_TRUE_CELL{k} = theta_TRUE(((k-1)*(M+1)+1):(k*(M+1)),:);
        theta_EST_CELL{k} = theta(((k-1)*(M+1)+1):(k*(M+1)),:);
        theta_START_CELL{k} = starting_point(((k-1)*(M+1)+1):(k*(M+1)),:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    objective_fun_identifiable_no_prop_CONTROLLED(theta_TRUE,sequences)
    objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)
    
    
    
    
    
    
    
    %%%%%%%% Finding which cluster got mapped to which one %%%%%%%%%%%%%%%%%%%%
    diff_with_clus = zeros(K,K);
    
    min_distance_index = zeros(K,1);
    for k = 1:K
        for h = 1:K
            diff_with_clus(k,h) = norm(theta_TRUE_CELL{k} - theta_EST_CELL{h}, 'fro');
        end
    end
    diff_with_clus
    
    for k = 1:K
        [MM,II] = min(diff_with_clus(k,:));
        min_distance_index(k) = II;
        diff_with_clus(:,II) = 99999;
    end
    
    %%%%%%%% Theta rearrange %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    THETA_CELL = cell(K,1);
    THETA_REARRANGED_CELL = cell(K,1);
    for k = 1:K
        THETA_CELL{k} = theta(((k-1)*(M+1)+1): (k*(M+1)),:);
    end
    
    for k = 1:K
        THETA_REARRANGED_CELL{k} = THETA_CELL{min_distance_index(k)};
        if(k == 1)
            theta_re = THETA_REARRANGED_CELL{k};
        else
            theta_re = [theta_re;THETA_REARRANGED_CELL{k}];
        end
    end
    
    theta_re - theta_TRUE
    
    
    %%%%%%%%%% Finding Predicted Clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Likelihoods = zeros(N,K);
    cluster_index_PREDICTED = zeros(N,1);
    M = length(unique(MAT(:,2)));
    for i = 1:N
        for k = 1:K    % K = num of clusters
            theta_temp =  theta_re(((M+1)*(k-1)+1):(k*(M+1)),:);
            initial_dist_temp = theta_temp(1,:);
            matrice_temp = theta_temp(2:(M+1),:);
            m = M * (sequences{i}(2:end) - 1) + sequences{i}(1:(end - 1));
            temp = initial_dist_temp(sequences{i}(1)) * prod(matrice_temp(m));
            Likelihoods(i,k) = temp;
        end
        [Max,I] = max(Likelihoods(i,:));
        cluster_index_PREDICTED(i) = I;
    end
    % %%%% Finding mapping of original cluster 1 and naming that cluster as %%%%%
    % %%%% cluster 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % [MM,II] = min(diff_with_clus(1,:));
    %
    % if(II ~= 1)
    %     cluster_index_PREDICTED(cluster_index_PREDICTED == II) = 999;
    %     cluster_index_PREDICTED(cluster_index_PREDICTED == 1) = II;
    %     cluster_index_PREDICTED(cluster_index_PREDICTED == 999) = 1;
    % end
    
    %%%%%% Finding initial estimate of Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fun = @(Y) value_at_beta_vec(Y,X,theta_re,sequences);
    
    beta = zeros(K,p);
    beta_vec = beta_2_beta_vec(beta);
    fun(beta_vec)
    
    x0 = beta_vec;
    lb = []; %-1*ones((K-1)*p,1);
    ub = []; %1*ones((K-1)*p,1);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    [x,fval] = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);
    
    beta_vec_new = x;
    beta_new = beta_vec_2_beta(beta_vec_new,K);
    
    value_at_beta_vec(beta_vec_new,X,theta_re,sequences)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    prop_MAT(X,beta_new);
    
    membership_prob_EST = zeros(N,K);
    cluster_index_EST = zeros(N,1);
    
    for i = 1:N
        membership_factors = zeros(1,K);
        for k = 1:K
            membership_factors(k) = exp(X(i,:)*beta_new(k,:)');
        end
        membership_prob_EST(i,:) = membership_factors/sum(membership_factors);
        [MM, II] = max(membership_prob_EST(i,:));
        cluster_index_EST(i) = II;
    end
    proportions_EST = zeros(K,1);
    for k = 1:K
        proportions_EST(k) = sum(cluster_index_EST == k)/N;
    end
    % proportions_EST
    % proportions_TRUE
    
    
    
    theta_start = theta_re;
    beta_start = beta_new;
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
                if(min(ge(theta(j,:),0)) == 0 || min(ge(theta(j,:),0)) == 0 )
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
    
    beta_end = beta_vec_2_beta(beta_vec,K)
    
    membership_prob_EST = zeros(N,K);
    
    for i = 1:N
        membership_factors = zeros(1,K);
        for k = 1:K
            membership_factors(k) = exp(X(i,:)*beta_end(k,:)');
        end
        membership_prob_EST(i,:) = membership_factors/sum(membership_factors);
        [MM, II] = max(membership_prob_EST(i,:));
        %cluster_index_EST(i) = II;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%% Finding Predicted Clusters USING SEQUENCES ONLY %%%%%%%%%%%%%%%%
    Likelihoods_only_seqs = zeros(N,K);
    max_index = zeros(N,1);
    %M = length(unique(MAT(:,2)));
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
    
    diff = TRUE_membership - ESTIMATED_membership;
    
    TPR = sum(diff == 0)/N;
    
    filename = ['ZZZ_TPR_clus_',num2str(K),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename,TPR)
    
    proportions_EST_from_cov_and_seq = zeros(K,1);
    for k = 1:K
        proportions_EST_from_cov_and_seq(k) = sum(max_index == k)/N;
    end
    
    
    filename = ['ZZZ_EST_proportions_',num2str(K),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename,proportions_EST_from_cov_and_seq)
    
    filename = ['ZZZ_BETA_EST_num_clus_',num2str(K),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename,round(beta_end,2))
    
    
end
