function value = value_at_beta_vec(beta_vec,X,theta,sequences)
M = size(theta,2);
K = size(theta,1)/(M+1);
p = length(beta_vec)/(K-1);
N = size(X,1);

% Beta_mat = zeros(K,p);
% Beta_mat(1,:) = zeros(1,p);
% count = 1;
% for ii = 2:K
%     for jj = 1:p
%         Beta_mat(ii,jj) = beta_vec(count);
%         count = count + 1;
%     end
% end
Beta_mat = beta_vec_2_beta(beta_vec,K);
membership_prob = zeros(N,K);

for i = 1:N
    membership_factors = zeros(1,K);
    for k = 1:K
        membership_factors(k) = exp(X(i,:)*Beta_mat(k,:)');
    end
    membership_prob(i,:) = membership_factors/sum(membership_factors);
end

value = objective_fun_identifiable_WITH_prop_at_BETA(theta,sequences,membership_prob);

end

% BETA = Beta_TRUE



