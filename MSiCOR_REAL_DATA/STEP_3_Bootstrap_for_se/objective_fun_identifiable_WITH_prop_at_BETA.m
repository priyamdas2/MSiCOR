% theta is K(M+1) X M matrix where in each (M+1) X M matrix, first row
% is the initial distribution vector, rest is the transsition matrix
function value = objective_fun_identifiable_WITH_prop_at_BETA(theta,sequences,membership_prob)
    M = size(theta,2);
    K = size(theta,1)/(M+1);
    
    initial_dists = cell(1,K);
    matrices = cell(1,K);
    
    for i = 1:K
        chunk = theta(((i-1)*(M+1)+1) : i*(M+1),:);
        initial_dists{i} = chunk(1,:);
        matrices{i} = chunk(2:end,:);
    end
    value = -estimate_likelihood_WITH_propMAT_CONTROLLED(matrices, initial_dists, sequences, membership_prob);
end

    


