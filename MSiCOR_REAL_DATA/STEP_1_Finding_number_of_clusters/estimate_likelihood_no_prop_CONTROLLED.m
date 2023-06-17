
function log_likelihood = estimate_likelihood_no_prop_CONTROLLED( matrices, initial_dists, sequences)

    K = size(sequences, 2); % gives number of sequences
    L = size(initial_dists, 2); % number of clusters
    n = length(initial_dists{1});
    
    cache = 0;

    
    for k = 1:K
        
        m = n * (sequences{k}(2:end) - 1) + sequences{k}(1:(end - 1));
        
        temps = zeros(L,1);
        for l = 1:L
            
            temps(l) = (1/L)*initial_dists{l}(sequences{k}(1)) * prod(matrices{l}(m));
            
        end
        
        temp = sum(temps);
        if(temp == 0)  % if sequence too large, may become zero
            temp =  10^(-323);
        end
        
        cache = cache + log(temp);
    end

    log_likelihood = cache ;



end