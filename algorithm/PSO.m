function [next_pop, next_fit, next_con, pbest_pop, pbest_fit, pbest_con, v, evals, prob] = PSO(pop, fit, con, v, v_min, v_max, pbest_pop, pbest_fit, pbest_con, evals, prob, cht, algRand)
    [NP, D] = size(pop);
    
    count = min(NP, prob.rest);
    pop = pop(1:count, :);
    fit = fit(1:count, :);
    con = con(1:count, :);
    pbest_pop = pbest_pop(1:count, :);
    pbest_fit = pbest_fit(1:count);
    pbest_con = pbest_con(1:count);
    v = v(1:count, :);
    
    w = 0.729;
    c1 = 2.05 * w;
    c2 = 2.05 * w;
    v = w * v + c1 .* rand(algRand, count, D) .* (pbest_pop - pop) + c2 .* rand(algRand, count, D) .* (repmat(pbest_pop(1, :), count, 1) - pop);
    v = max(v, v_min);
    v = min(v, v_max);
    next_pop = pop + v;
    next_pop = boundary_check(next_pop, prob.lb, prob.ub);      % check the boundary
    v = next_pop - pop;                                         % reset the velocity
    [next_fit, next_con, prob, evals] = evaluate(prob, next_pop, evals, true);
    
%     select = next_fit > pbest_fit;
    select = cmp_indis(pbest_fit, next_fit, pbest_con, next_con, cht);
    pbest_pop = repmat(select, 1, D) .* next_pop + repmat(1-select, 1, D) .* pbest_pop;
    pbest_fit = select .* next_fit + (1-select) .* pbest_fit;
    pbest_con = select .* next_con + (1-select) .* pbest_con;
end