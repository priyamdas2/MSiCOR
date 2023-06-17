clear all
rng(1)
% Consider K many markov chains, each having M many states, so there will
% be K many MxM transition matrices and one vector of length M (proportion
% vector). Here 'theta' denotes K many (M+1)x M matrices stacked one after
% another. So 'theta' is of ((M+1)K)x K diemnsional. 'theta_prop' denotes the
% proportion vetor of Kx1 diemnsion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10; % number of states
K = 6; % number of clusters
% Covariate order : c(Intercept, "age_at_diag", "disease_duration", "gender","white", "black")
filename = ['Covariate_values_ordered.csv'];
X = csvread(filename);
p = size(X,2); % number of variables
N = size(X,1);
filename = ['MSiCOR_DATA_after_2006_10_codes_have_covariates.csv'];
MAT = csvread(filename);

N

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

Indicator_matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_states = zeros(N,1);
for gg = 1:N % = size(sequences,2)
    initial_states(gg) = sequences{gg}(1);
end

length(initial_states(initial_states  == 1)) 
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
dirich_par = ones(1,M);
starting_point = drchrnd(dirich_par, (M+1)*K);
%starting_point =  theta_TRUE;

%theta_prop = starting_point_prop;
theta = starting_point;


objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)
tic;

no_iter = 5000;
lh_array = zeros(no_iter,1);
tol_fun = 10^(-1);                 % minimum improvement to keep epsilon same after iteration
tol_fun_2 = 10^(-1);               % minimum improvement to start another run
parameter_cut_off = 10^(-20);      % hard-thresholding of each coordinate value
epsilon_cut_off = 10^(-2);         % smallest step-size
epsilon_decreasing_factor_1 = 2;   % step-decay rate (default is 2)
epsilon_decreasing_factor_2 = 2;
max_runs = 1;    % Number of runs

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


%%%%%%%%%%%%%%%%%%%%%%%%
objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)


%%%%%%%% Theta rearrange %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
THETA_CELL = cell(K,1);
THETA_REARRANGED_CELL = cell(K,1);
for k = 1:K
    THETA_CELL{k} = theta(((k-1)*(M+1)+1): (k*(M+1)),:);
end

for k = 1:K
    if(k == 1)
        theta_re = THETA_CELL{k};
    else
        theta_re = [theta_re;THETA_CELL{k}];
    end
end




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
proportions_EST


filename = ['ZZZ_Theta_initial_value_num_clus_',num2str(K),'.csv'];
csvwrite(filename,theta_re)

filename = ['ZZZ_Beta_initial_value_num_clus_',num2str(K),'.csv'];
csvwrite(filename,beta_new)

beta_new

% % filename = ['Final_NEWWW_algo_rastrigin_high_dim_', num2str(M-1),'_' ,num2str(n),'.xls'];
% % xlswrite(filename,[das_algo_values(1),das_algo_times(1)])