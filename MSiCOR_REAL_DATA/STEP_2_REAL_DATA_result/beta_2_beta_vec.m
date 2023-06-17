function value = beta_2_beta_vec(beta)
K = size(beta,1);
p = size(beta,2);

BETA_trans = beta';
BETA_vec_whole = BETA_trans(:);
BETA_vec = BETA_vec_whole;
BETA_vec(1:p) = [];
value = BETA_vec;
end