function value = prop_MAT(X,beta)
N = size(X,1);
p = size(X,2);
K = size(beta,1);
propMAT = zeros(N,K);
exp_xBeta = exp(X*beta');
for i = 1:N
    for j = 1:K
        propMAT(i,j) = exp_xBeta(i,j)/sum(exp_xBeta(i,:));
    end
end

value = propMAT;
end