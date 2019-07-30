function [next_pop, next_fit, next_con, evals, prob] = CSA1(pop, fit, con, dis, copy, evals, prob, cht, algRand)
    [NP, D] = size(pop);
    c = max(dis(:));
    if NP == 1
        c = 1;
    end
    next_pop = zeros(NP, D);
    next_fit = zeros(NP, 1);
    next_con = zeros(NP, 1);
    alpha = 0.5 * c / sqrt(D);                  % alpha
    
    for i = 1:NP
        anti = pop(i, :);                       % current antibody
        pool = ones(copy+1, 1) * anti;          % clone pool
        pool_fit = ones(copy+1, 1) .* fit(i);   % clone fitnees
        pool_con = ones(copy+1, 1) .* con(i);
        if i == 1 % the seed uses Gaussian mutation
            pool = pool + [zeros(1, D); alpha .* randn(algRand, copy, D)];
        else
            if rand(algRand) < 0.5
                step = [zeros(1, D); 0.5 .* (pop(1, :) - anti); alpha .* randn(algRand, copy-1, D)];
                pool = pool + step;
            else
                pool = pool + [zeros(1, D); alpha .* randn(algRand, copy, D)];
            end
        end
        pool = boundary_check(pool, prob.lb, prob.ub);      % check the boundary
        
        if prob.rest > copy % selection
            [pool_fit(2:end), pool_con(2:end), prob, evals] = evaluate(prob, pool(2:end, :), evals, true); % evaluate the children
            sort_idx = sort_indis(pool_fit, pool_con, cht);
            best_idx = sort_idx(1);
            next_pop(i, :) = pool(best_idx, :);
            next_fit(i) = pool_fit(best_idx);
            next_con(i) = pool_con(best_idx);         
        else % prob.rest <= copy
            count = prob.rest;
            [pool_fit(2:1+count), pool_con(2:1+count), prob, evals] = evaluate(prob, pool(2:1+count, :), evals, true); 
            sort_idx = sort_indis(pool_fit(1:1+count), pool_con(1:1+count), cht);
            best_idx = sort_idx(1);
            next_pop(i, :) = pool(best_idx, :);
            next_fit(i) = pool_fit(best_idx);
            next_con(i) = pool_con(best_idx);
            
            next_pop = next_pop(1:i, :);
            next_fit = next_fit(1:i);
            next_con = next_con(1:i);
            break;
        end
    end
end