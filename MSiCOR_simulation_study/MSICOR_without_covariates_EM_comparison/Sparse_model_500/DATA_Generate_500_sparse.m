clear all
rng(1)
% Consider K many markov chains, each having M many states, so there will
% be K many MxM transition matrices and one vector of length M (proportion
% vector). Here 'theta' denotes K many (M+1)x M matrices stacked one after
% another. So 'theta' is of ((M+1)K)x K diemnsional. 'theta_prop' denotes the
% proportion vetor of Kx1 diemnsion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10; % number of states
K = 3; % number of clusters
N = 500;
num_rep = 50;
%%%%%% Generating Transition matrix+initial prob vectors %%%%%%%%%%%%%%%%%

Transitions = cell(K,1);
initial_DIST_TRUE = cell(K,1);

dirich_par = ones(1,M);

for k = 1:K
    random_mat = drchrnd(dirich_par, M);
    if (k == 1)
        num_zero = .3*M*M;
    end
    if (k == 2)
        num_zero = .5*M*M;
    end
    if(k == 3)
        num_zero = .7*M*M;
    end
    
    if(num_zero ~= 0)
        zero_locations = randsample(M*M,num_zero);
        random_mat(zero_locations) = 0;
    end
    
    for vv = 1:M
        random_mat(vv,:) = random_mat(vv,:)/sum(random_mat(vv,:));
    end
            
        
    Transitions{k} = random_mat;
    initial_DIST_TRUE{k} = drchrnd(dirich_par, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Finding Membership %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

membership_prob_TRUE = zeros(N,K);
cluster_index_DETERMINISTIC = zeros(N,1);

for i = 1:N
    membership_factors = zeros(1,K);
    for k = 1:K
        membership_factors(k) = unifrnd(0,1);
    end
    membership_prob_TRUE(i,:) = membership_factors/sum(membership_factors);
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

filename = ['AAA_Theta_TRUE.csv'];
csvwrite(filename,theta_TRUE)

%%%%%%%%%% DATA GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rep = 1:num_rep
    
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
    
    filename = ['ZZZ_DATA_sequences_rep_',num2str(rep),'.csv'];
    csvwrite(filename,MAT)
    
    
    filename = ['ZZZ_TRUE_membership_rep_',num2str(rep),'.csv'];
    csvwrite(filename,TRUE_membership)
    
    
    
end

% filename = ['ZZZ_TRUE_proportions_',num2str(K),'.csv'];
% csvwrite(filename,proportions_TRUE)
