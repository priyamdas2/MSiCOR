function value = beta_vec_2_beta(beta_vec,K)
p = length(beta_vec)/(K-1);
beta = zeros(K,p);
count = 1;
for k = 2:K
    for j = 1:p
        beta(k,j) = beta_vec(count);
        count = count + 1;
    end
end
value = beta;
end
