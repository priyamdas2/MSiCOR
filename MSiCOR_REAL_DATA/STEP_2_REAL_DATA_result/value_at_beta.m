function value = value_at_beta(beta,X,theta,sequences)
M = size(theta,2);
K = size(theta,1)/(M+1);
N = size(X,1);
Beta_mat = beta;
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



