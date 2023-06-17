clear all
%%% CHECK EPSILON
% Consider K many markov chains, each having M many states, so there will
% be K many MxM transition matrices and one vector of length K (proportion
% vector). Here 'theta' denotes K many MxM matrices stacked one after
% another. So 'theta' is of (MK)x K diemnsional. 'theta_prop' denotes the
% proportion vetor of Kx1 diemnsion.

%%%%%%%%%%%%%%%%%%%%%%%%% READ DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_rep = 50;


TVD_array = zeros(num_rep,1);
TPR_array = zeros(num_rep,1);

parfor rep = 1:num_rep
    
    M = 10;
    K = 3; % number of mixture components
    n = (M+1)*K;
    
    filename = ['ZZZ_DATA_sequences_rep_',num2str(rep),'.csv'];
    MAT = csvread(filename);
    filename = ['ZZZ_TRUE_membership_rep_',num2str(rep),'.csv'];
    TRUE_membership = csvread(filename);
    THETA_TRUE = csvread('AAA_Theta_TRUE.csv');
    
    num_unique_patients = size(unique(MAT(:,1)),1);
    sequences_temp = cell(1,num_unique_patients);
    
    patient_num = 1;
    sequences_temp{1} = MAT(1,2);
    seq_lengths = zeros(num_unique_patients,1);  % check
    
    for i = 2:size(MAT,1)
        if(MAT(i,1) == MAT(i-1,1))
            sequences_temp{patient_num} = [sequences_temp{patient_num} MAT(i,2)];
        else
            patient_num = patient_num + 1;
            sequences_temp{patient_num} = MAT(i,2);
        end
    end
    
    sequences = sequences_temp;
    
    sequence_lengths = zeros(length(sequences),1);
    for i = 1:length(sequences)
        sequence_lengths(i) = length(sequences{i});
    end
    
    %%%%%%%%%% Finding transition places %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Indicator_matrix = zeros(M,M);
    for i = 1:num_unique_patients
        for j = 1:(length(sequences{i})-1)
            aa = sequences{i}(j);
            bb = sequences{i}(j+1);
            Indicator_matrix(aa,bb) = Indicator_matrix(aa,bb) + 1;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    initial_states = zeros(num_unique_patients,1);
    for gg = 1:num_unique_patients % = size(sequences,2)
        initial_states(gg) = sequences{gg}(1);
    end
    
    length(initial_states(initial_states  == 1)) % 4871 in INF based
    tabulate(initial_states)
    
    % Indicator_matrix
    % tabulate(MAT(:,2))
    
    sum(sum(Indicator_matrix))
    
    for j = 1:num_unique_patients
        seq_lengths(j) = size(sequences_temp{j},2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIMIZATION INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rng(rep)
    
    %starting_point = (1/M)*ones(n,M); % M by n matrix
    
    starting_point_temp = rand(n,M);
    for i = 1:n
        starting_point_temp(i,:) = starting_point_temp(i,:)/sum(starting_point_temp(i,:));
    end
    starting_point = starting_point_temp;
    theta = starting_point;
    objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)
    
    tic;
    
    no_iter = 5000;
    tol_fun = 10^(1);                 % minimum improvement to keep epsilon same after iteration
    tol_fun_2 = 10^(1);              % minimum improvement to start another run
    parameter_cut_off = 10^(-20);       % hard-thresholding of each coordinate value
    epsilon_cut_off = 5*10^(-3);        % smallest step-size
    epsilon_decreasing_factor_1 = 2; % step-decay rate (default is 2)
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
                if(min(ge(theta(j,:),0)) == 0)
                    stop('error')
                end
                if(sum(theta(j,:))<0.99 || sum(theta(j,:))>1.01)
                    stop('error')
                end
                [ii,i,j,log(current_lh_1), epsilon]
                %[num,current_lh_1]
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
            % Simultaneous Best move selection + sparsity control
            
            %[num,current_lh_1]
            
            if(num < current_lh_1)
                parameter_1 = check_each_best(idx,:);
                %             sparsity_positions = lt(parameter_1,parameter_cut_off*ones(1,M));
                %             garbage = sum(parameter_1(sparsity_positions));
                %             if(garbage > 0)
                %                 parameter_1(sparsity_positions) = 0;
                %                 rich_positions = ge(parameter_1,parameter_cut_off*ones(1,M));
                %                 parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
                %             end
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
    
    % objective_fun_identifiable_no_prop_CONTROLLED(starting_point,sequences)
    % objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)
    
    %%% Each patient belongs to which cluster
    
    Likelihoods = zeros(num_unique_patients,K);
    max_index = zeros(num_unique_patients,1);
    M = length(unique(MAT(:,2)));
    for i = 1:num_unique_patients
        for k = 1:K    % K = num of clusters
            theta_temp =  theta(((M+1)*(k-1)+1):(k*(M+1)),:);
            initial_dist_temp = theta_temp(1,:);
            matrice_temp = theta_temp(2:(M+1),:);
            m = M * (sequences{i}(2:end) - 1) + sequences{i}(1:(end - 1));
            temp = initial_dist_temp(sequences{i}(1)) * prod(matrice_temp(m));
            Likelihoods(i,k) = temp;
        end
        [Max,I] = max(Likelihoods(i,:));
        max_index(i) = I;
    end
    kkk = tabulate(max_index);
    PROPORTIONS = kkk(:,3)/100;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WARM STARTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    starting_point = theta;
    starting_proportions = PROPORTIONS;
    
    theta = starting_point;
    proportions = starting_proportions';
    
    % objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,proportions)
    % objective_fun_identifiable_no_prop_CONTROLLED(theta,sequences)
    
    tic;
    
    no_iter = 5000;
    tol_fun = 10^(-1);                 % minimum improvement to keep epsilon same after iteration
    tol_fun_2 = 10^(-1);              % minimum improvement to start another run
    parameter_cut_off = 10^(-3);       % hard-thresholding of each coordinate value
    epsilon_cut_off = 5*10^(-4);        % smallest step-size
    epsilon_decreasing_factor_1 = 2; % step-decay rate (default is 2)
    epsilon_decreasing_factor_2 = 2;
    max_runs = 10;    % Number of runs
    
    like_corresponding_run = zeros(max_runs,1);
    
    
    for ii = 1:max_runs
        old_theta = theta;
        old_proportions = proportions;
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
            
            current_lh_1 = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,proportions);
            if(i == 1)
                num = current_lh_1;
            end
            %%%%%%%%%%%%%%%%%%%% UPDATING TRANSITION MATRICES %%%%%%%%%%%%%%%%%%
            for j = 1:n
                if(min(ge(theta(j,:),0)) == 0)
                    stop('error')
                end
                if(sum(theta(j,:))<0.99 || sum(theta(j,:))>1.01)
                    stop('error')
                end
                [ii,i,j,log(current_lh_1), epsilon]
                proportions
                %[num,current_lh_1]
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
                            total_lh(parallel_no) = objective_fun_identifiable_WITH_prop_CONTROLLED(proxy_coef_theta_pos,sequences,proportions);
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
                            total_lh(parallel_no) = objective_fun_identifiable_WITH_prop_CONTROLLED(proxy_coef_theta_neg,sequences,proportions);
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
            % Simultaneous Best move selection + sparsity control
            
            %[num,current_lh_1]
            
            if(num < current_lh_1)
                parameter_1 = check_each_best(idx,:);
                sparsity_positions = lt(parameter_1,parameter_cut_off*ones(1,M));
                garbage = sum(parameter_1(sparsity_positions));
                if(garbage > 0)
                    parameter_1(sparsity_positions) = 0;
                    %                 rich_positions = ge(parameter_1,parameter_cut_off*ones(1,M));
                    %                 parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
                    parameter_1 = parameter_1/sum(parameter_1);
                end
                theta_copy = theta;
                theta_copy(idx,:) = parameter_1;
                sparse_likelihood = objective_fun_identifiable_WITH_prop_CONTROLLED(theta_copy,sequences,proportions);
                if(sparse_likelihood < num)
                    theta(idx,:) = parameter_1;
                else
                    theta(idx,:) = check_each_best(idx,:);
                end
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%% UPDATING PROPROTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            current_lh_2 = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,proportions);
            
            total_lh_proportions = zeros(2*K,1);
            matrix_update_at_h_proportions = zeros(2*K,K);
            
            for parallel_no = 1:(2*K)   % CAN BE RUN IN PARALLEL
                % pos
                if(mod(parallel_no,2) == 1)
                    positive_change_loc = (parallel_no + 1)/2;
                    possibility_pos = proportions;
                    temp_possibility_pos = proportions;
                    temp_possibility_pos(positive_change_loc) = 0; % To find all significant positions except the positive_change_loc-th
                    significant_positions = find(gt(temp_possibility_pos, parameter_cut_off*ones(1,K)));
                    if(isempty(significant_positions) == 1)
                        possibility_pos = proportions;
                    else
                        possibility_pos(positive_change_loc) = proportions(positive_change_loc) + epsilon;
                        possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon/length(significant_positions);
                        epsilon_temp_pos = epsilon;
                        
                        if(min(ge(possibility_pos,0)) == 0 && epsilon_temp_pos > epsilon_cut_off)
                            epsilon_temp_pos = epsilon_temp_pos/epsilon_decreasing_factor;
                            possibility_pos = proportions;
                            possibility_pos(positive_change_loc) = proportions(positive_change_loc) + epsilon_temp_pos;
                            possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon_temp_pos/length(significant_positions);
                        end
                    end
                    if(sum(possibility_pos) <0.99 || sum(possibility_pos)>1.01)  %% Checking if sum is 1
                        stop('error')
                    end
                    
                    if(min(ge(possibility_pos,0)) == 0 || isequal(possibility_pos,proportions) == 1)
                        possibility_pos = proportions;
                        total_lh_proportions(parallel_no) = current_lh_2;
                    else
                        total_lh_proportions(parallel_no) = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,possibility_pos);
                    end
                    matrix_update_at_h_proportions(parallel_no,:) = possibility_pos;
                else
                    % neg
                    negative_change_loc = parallel_no/2;
                    possibility_neg = proportions;
                    temp_possibility_neg = proportions;
                    temp_possibility_neg(negative_change_loc) = 0;
                    significant_positions = find(gt(temp_possibility_neg, parameter_cut_off*ones(1,K)));
                    if(isempty(significant_positions) == 1)
                        possibility_neg = proportions;
                    else
                        possibility_neg(negative_change_loc) = proportions(negative_change_loc) - epsilon;
                        possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon/length(significant_positions);
                        epsilon_temp_neg = epsilon;
                        
                        if(min(ge(possibility_neg,0)) == 0 && epsilon_temp_neg > epsilon_cut_off)
                            epsilon_temp_neg = epsilon_temp_neg/epsilon_decreasing_factor;
                            possibility_neg = proportions;
                            possibility_neg(negative_change_loc) = proportions(negative_change_loc) - epsilon_temp_neg;
                            possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon_temp_neg/length(significant_positions);
                        end
                    end
                    if(sum(possibility_neg) <0.99 || sum(possibility_neg)>1.01)  %% Checking if sum is 1
                        stop('error')
                    end
                    if(min(ge(possibility_neg,0)) == 0 || isequal(possibility_neg,proportions) == 1)
                        possibility_neg = proportions;
                        total_lh_proportions(parallel_no) = current_lh_2;
                    else
                        total_lh_proportions(parallel_no) = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,possibility_neg);
                    end
                    matrix_update_at_h_proportions(parallel_no,:) = possibility_neg;
                end
            end
            
            
            [Min_here_proportions,I_proportions] = min(total_lh_proportions);
            check_each_best_proportions = proportions;
            corresponding_likelihood_proportions = current_lh_2;
            
            if(Min_here_proportions < current_lh_2)
                count_1 = count_1+1;
                proportions = matrix_update_at_h_proportions(I_proportions,:);
            end
            
            
            
            
            array_of_values(i) = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences,proportions);
            
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
    
%     das_algo_times = toc;
%     das_algo_values = objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences, proportions);
    
    objective_fun_identifiable_WITH_prop_CONTROLLED(starting_point,sequences,starting_proportions)
    objective_fun_identifiable_WITH_prop_CONTROLLED(theta,sequences, proportions)
    
    
    %%%%%%%%%%%%%% Finding right order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    THETA_NOT_ORDERED = theta;
    
    
    
    permutations = [1:11,12:22,23:33;
        1:11,23:33,12:22;
        12:22,1:11,23:33;
        12:22,23:33,1:11;
        23:33,1:11,12:22;
        23:33,12:22,1:11];
    perm_Clus = [1,2,3;
        1,3,2;
        2,1,3;
        2,3,1;
        3,1,2;
        3,2,1];
    
    values = zeros(6,1);
    for test = 1:6
        values(test) = sum(sum(abs(THETA_TRUE - THETA_NOT_ORDERED(permutations(test,:),:))));
    end
    
    [aa,bb] = min(values);
    
    THETA_ORDERED = THETA_NOT_ORDERED(permutations(bb,:),:);
    proportions_ORDERED = proportions(perm_Clus(bb,:));
    
    %%% Each patient belongs to which cluster
    
    Likelihoods = zeros(num_unique_patients,K);
    max_index = zeros(num_unique_patients,1);
    M = length(unique(MAT(:,2)));
    for i = 1:num_unique_patients
        for k = 1:K    % K = num of clusters
            theta_temp =  THETA_ORDERED(((M+1)*(k-1)+1):(k*(M+1)),:);
            initial_dist_temp = theta_temp(1,:);
            matrice_temp = theta_temp(2:(M+1),:);
            m = M * (sequences{i}(2:end) - 1) + sequences{i}(1:(end - 1));
            temp = initial_dist_temp(sequences{i}(1)) * prod(matrice_temp(m));
            Likelihoods(i,k) = proportions_ORDERED(k)*temp;
        end
        [Max,I] = max(Likelihoods(i,:));
        max_index(i) = I;
    end
    kkk = tabulate(max_index);
    
    TVD_array(rep) = aa;
    TPR_array(rep) = sum((TRUE_membership - max_index) == 0)/num_unique_patients;
    
end

output = zeros(2,2);
output(1,1) = mean(TVD_array);
output(1,2) = std(TVD_array)/sqrt(length(TVD_array));
output(2,1) = mean(TPR_array);
output(2,2) = std(TPR_array)/sqrt(length(TPR_array));

das_algo_times = toc;
das_algo_times


filename = ['AAA_MSICOR_all_TVD_TPR.csv'];
csvwrite(filename,[TVD_array,TPR_array]);


filename = ['AAA_MSICOR_output.csv'];
csvwrite(filename,output)
% 
% filename = ['AAA_Theta_output', num2str(M-1),'_' ,num2str(n),'.xls'];
% xlswrite(filename,[das_algo_values(1),das_algo_times(1)])