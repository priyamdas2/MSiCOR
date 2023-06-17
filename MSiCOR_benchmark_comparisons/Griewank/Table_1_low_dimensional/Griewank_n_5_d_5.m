clear all

M = 6; % d = M-1
n = 5;

sim_reps = 100;

one_by_root_i_start = 1./sqrt([1:(M-1)]');
one_by_root_i = repmat(one_by_root_i_start,1,n);
lb_here = (-500)*ones(M-1,n);
ub_here = (500)*ones(M-1,n);


fun = @(x) sum((1/4000)*sum(((M-1)*x(1:(M-1),:).*(ub_here-lb_here)+lb_here).^2)) - sum(prod(cos(((M-1)*x(1:(M-1),:).*(ub_here-lb_here)+lb_here).*one_by_root_i))) + n;


one_by_root_i_ga = repmat(one_by_root_i_start',1,n);
lb_here_ga = (-500)*ones(1,n*(M-1));
ub_here_ga = (500)*ones(1,n*(M-1));

fun_ga = @(xxx) (1/4000)*sum(((M-1)*xxx.*(ub_here_ga-lb_here_ga)+lb_here_ga).^2) - prod(cos(((M-1)*xxx.*(ub_here_ga-lb_here_ga)+lb_here_ga).*(one_by_root_i_ga))) + 1;


solution = vertcat((1/(2*(M-1)))*ones(M-1,n),0.5*ones(1,n));
solution_ga = (1/(2*(M-1)))*ones(1,n*(M-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZATION INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aequal = zeros(n,M*n);

for i = 1:n
    for j = ((i-1)*M+1):(i*M)
        Aequal(i,j) = 1;
    end
end

A_here = zeros(n,(M-1)*n);

for i = 1:n
    for j = ((i-1)*(M-1)+1):(i*(M-1))
        A_here(i,j) = 1;
    end
end




das_algo_values = zeros(sim_reps,1);
das_algo_times = zeros(sim_reps,1);

fmincon_values_ip = zeros(sim_reps,1);
fmincon_times_ip = zeros(sim_reps,1);

fmincon_values_sqp = zeros(sim_reps,1);
fmincon_times_sqp = zeros(sim_reps,1);

ga_values = zeros(sim_reps,1);
ga_times = zeros(sim_reps,1);


for jjj = 1: sim_reps
    %%% DAS ALGO
    
    starting_point_trans = drchrnd(ones(1,M), n);
    starting_point = starting_point_trans';
    
    
    theta = starting_point;
    
    
    
    
    tic;
    
    no_iter = 5000;
    lh_array = zeros(no_iter,1);
    tol_fun = 10^(-6);
    tol_fun_2 = 10^(-6);
    parameter_cut_off = 10^(-6);
    epsilon_cut_off = 10^(-4);
    epsilon_decreasing_factor_1 = 1.01; % default is 2
    epsilon_decreasing_factor_2 = 1.01;
    max_runs = 200;
    
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
            check_each_best = zeros(M,n);
            corresponding_likelihood = zeros(n,1);
            
            current_lh_1 = fun(theta);
            
            for j = 1:n
                if(min(ge(theta(:,j),0)) == 0)
                    stop('error')
                end
                [jjj,ii,i,j,current_lh_1, epsilon,current_lh_1]
                
                total_lh_pos = zeros(1,M);
                total_lh_neg = zeros(1,M);
                
                matrix_pos_update_at_h = zeros(M,M);
                matrix_neg_update_at_h = zeros(M,M);
                
                % pos
                for positive_change_loc = 1:M
                    possibility_pos = theta(:,j);
                    temp_possibility_pos = theta(:,j);
                    temp_possibility_pos(positive_change_loc) = 0; % To find all significant positions except the positive_change_loc-th
                    significant_positions = find(gt(temp_possibility_pos, parameter_cut_off*ones(M,1)));
                    if(isempty(significant_positions) == 1)
                        possibility_pos = theta(:,j);
                    else
                        possibility_pos(positive_change_loc) = theta(positive_change_loc,j) + epsilon;
                        possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon/size(significant_positions,1);
                        epsilon_temp_pos = epsilon;
                        
                        if(min(ge(possibility_pos,0)) == 0 && epsilon_temp_pos > epsilon_cut_off)
                            epsilon_temp_pos = epsilon_temp_pos/epsilon_decreasing_factor;
                            possibility_pos = theta(:,j);
                            possibility_pos(positive_change_loc) = theta(positive_change_loc,j) + epsilon_temp_pos;
                            possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon_temp_pos/size(significant_positions,1);
                        end
                    end
                    
                    if(min(ge(possibility_pos,0)) == 0 || isequal(possibility_pos,theta(:,j)) == 1)
                        possibility_pos = theta(:,j);
                        total_lh_pos(positive_change_loc) = current_lh_1;
                    else
                        proxy_coef_theta_pos = theta;
                        proxy_coef_theta_pos(:,j) = possibility_pos;
                        total_lh_pos(positive_change_loc) = fun(proxy_coef_theta_pos);
                    end
                    
                    matrix_pos_update_at_h(:,positive_change_loc) = possibility_pos;
                end
                
                % neg
                for negative_change_loc = 1:M
                    possibility_neg = theta(:,j);
                    temp_possibility_neg = theta(:,j);
                    temp_possibility_neg(negative_change_loc) = 0;
                    significant_positions = find(gt(temp_possibility_neg, parameter_cut_off*ones(M,1)));
                    if(isempty(significant_positions) == 1)
                        possibility_neg = theta(:,j);
                    else
                        possibility_neg(negative_change_loc) = theta(negative_change_loc,j) - epsilon;
                        possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon/size(significant_positions,1);
                        epsilon_temp_neg = epsilon;
                        
                        if(min(ge(possibility_neg,0)) == 0 && epsilon_temp_neg > epsilon_cut_off)
                            epsilon_temp_neg = epsilon_temp_neg/epsilon_decreasing_factor;
                            possibility_neg = theta(:,j);
                            possibility_neg(negative_change_loc) = theta(negative_change_loc,j) - epsilon_temp_neg;
                            possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon_temp_neg/size(significant_positions,1);
                        end
                    end
                    
                    if(min(ge(possibility_neg,0)) == 0 || isequal(possibility_neg,theta(:,j)) == 1)
                        possibility_neg = theta(:,j);
                        total_lh_neg(negative_change_loc) = current_lh_1;
                    else
                        proxy_coef_theta_neg = theta;
                        proxy_coef_theta_neg(:,j) = possibility_neg;
                        total_lh_neg(negative_change_loc) = fun(proxy_coef_theta_neg);
                    end
                    matrix_neg_update_at_h(:,negative_change_loc) = possibility_neg;
                end
                
                
                % general
                
                [M_pos_1,I_pos_1] = min(total_lh_pos);
                [M_neg_1,I_neg_1] = min(total_lh_neg);
                
                check_each_best(:,j) = theta(:,j);
                corresponding_likelihood(j) = current_lh_1;
                
                if(min(M_pos_1,M_neg_1) < current_lh_1)
                    count_1 = count_1+1;
                    if(M_pos_1 < M_neg_1)
                        check_each_best(:,j) = matrix_pos_update_at_h(:,I_pos_1);
                        corresponding_likelihood(j) = M_pos_1;
                    else
                        check_each_best(:,j) = matrix_neg_update_at_h(:,I_neg_1);
                        corresponding_likelihood(j) = M_neg_1;
                    end
                end
                
            end
            
            
            [num, idx] = min(corresponding_likelihood(:));
            
            
            % Best move selection + sparsity control
            
            parameter_1 = check_each_best(:,idx);
            sparsity_positions = lt(parameter_1,parameter_cut_off*ones(M,1));
            garbage = sum(parameter_1(sparsity_positions));
            if(garbage > 0)
                parameter_1(sparsity_positions) = 0;
                rich_positions = ge(parameter_1,parameter_cut_off*ones(M,1));
                parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
            end
            theta(:,idx) = parameter_1;
            
            
            array_of_values(i) = current_lh_1;
            
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
    
    das_algo_times (jjj) = toc;
    das_algo_values (jjj) = fun(theta);
    
    tic;
    x0 = starting_point;
    A = [];
    b = [];
    Aeq = Aequal;
    beq = ones(n,1);
    lb = zeros(M,n);
    ub = ones(M,n);
    options = struct('MaxFunEvals', Inf,'MaxIter',Inf,'UseParallel',false,'Algorithm','interior-point');
    x_3 = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[], options);
    fmincon_times_ip(jjj) = toc;
    fmincon_values_ip(jjj) = fun(x_3);
    
    
    
    
    % fmincon (sqp)
    
    tic;
    x0 = starting_point;
    A = [];
    b = [];
    Aeq = Aequal;
    beq = ones(n,1);
    lb = zeros(M,n);
    ub = ones(M,n);
    options = struct('MaxFunEvals', Inf,'MaxIter',Inf,'UseParallel',false,'Algorithm','sqp');
    x_3 = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[], options);
    fmincon_times_sqp(jjj) = toc;
    fmincon_values_sqp(jjj) = fun(x_3);
    
    
    % GA
    
    
    
    tic;
    X00 = starting_point(1:(M-1),:);
    X0 = X00(:)';
    A = A_here;
    b = ones(n,1);
    Aeq = [];
    beq = [];
    lb = zeros((M-1)*n,1);
    ub = [];
    options = gaoptimset(options,'UseParallel',false,'InitialPopulation',X0);
    x_ga = ga(fun_ga,size(X0,2),A,b,Aeq,beq,lb,ub,[],options);
    ga_times(jjj) = toc;
    ga_values(jjj) = fun_ga(x_ga);

    
    
end


filename = ['Griewank_output_MSiCOR_IP_SQP_GA_d_', num2str(M-1),'_n_' ,num2str(n),'.xls'];
xlswrite(filename,[min(das_algo_values),std(das_algo_values),mean(das_algo_times),std(das_algo_times);...
    min(fmincon_values_ip),std(fmincon_values_ip),mean(fmincon_times_ip),std(fmincon_times_ip); 
    min(fmincon_values_sqp),std(fmincon_values_sqp),mean(fmincon_times_sqp),std(fmincon_times_sqp); ...
    min(ga_values),std(ga_values),mean(ga_times),std(ga_times)]);







