function [next_pop, next_fit, next_con, next_eta, evals, prob] = FEP(pop, fit, con, eta, tau, tau1, q, evals, prob, cht, algRand);
    [NP, D] = size(pop);

    count = min(NP, prob.rest);

    offspring = pop(1:count, :) + eta(1:count, :) .* randn(algRand, count, D);
    offspring = boundary_check(offspring, prob.lb, prob.ub);      % check the boundary
    
    offspring_eta = eta(1:count, :) .* exp(tau1 .* repmat(randn(algRand, count, 1), 1, D) + tau .* randn(algRand, count, D));

    [offspring_fit, offspring_con, prob, evals] = evaluate(prob, offspring, evals, true);
    
    indis_pop = [pop; offspring];
    indis_fit = [fit; offspring_fit];
    indis_con = [con; offspring_con];
    indis_eta = [eta; offspring_eta];
    
    q = min(q, count);
    cmp_idx = zeros(size(indis_pop, 1), q);
    for i = 1:size(indis_pop, 1)
        cmp_idx(i, :) = randperm(algRand, count, q);
    end
    
    win = zeros(size(indis_pop, 1), 1);
    for i = 1:q
        win = win + cmp_indis(indis_fit(cmp_idx(:, i)), indis_fit, indis_con(cmp_idx(:, i)), indis_con, cht);
    end
    [~, win_idx] = sort(win, 'descend');
    select_idx = win_idx(1:count);
    
    next_pop = indis_pop(select_idx, :);
    next_fit = indis_fit(select_idx, :);
    next_con = indis_con(select_idx, :);
    next_eta = indis_eta(select_idx, :);
end