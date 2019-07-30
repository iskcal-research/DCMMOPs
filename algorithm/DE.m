function [next_pop, next_fit, next_con, evals, prob] = DE(pop, fit, con, evals, prob, cht, algRand)
    [NP, D] = size(pop);
    if NP < 4                       % if the size is less than 4, do nothing
        next_pop = pop;
        next_fit = fit;
        next_con = con;
        return;
    end
    count = min(NP, prob.rest);
    idx = zeros(count, 3); % get the index of DE/rand/1
    for i = 1:count
       idx(i, :) = randperm(algRand, NP-1, 3);
       idx(i, idx(i, :) >= i) = idx(i, idx(i, :) >= i) + 1;
    end
    offspring = pop(idx(:, 1), :) + 0.5 * (pop(idx(:, 2), :) - pop(idx(:, 3)));
    offspring = boundary_check(offspring, prob.lb, prob.ub);      % check the boundary
    
    cross = rand(algRand, count, D) < 0.9;
    for i = 1:count
        if sum(cross(i, :)) == 0
            cross(i, randi(algRand, D)) = true;
        end
    end
    offspring = cross .* offspring + (1 - cross) .* pop(1:count, :);
    [offspring_fit, offspring_con, prob, evals] = evaluate(prob, offspring, evals, true);
    
    select  = cmp_indis(fit(1:count), offspring_fit, con(1:count), offspring_con, cht);
    %offspring_fit > fit(1:count);
    next_pop = repmat(select, 1, D) .* offspring + repmat(1-select, 1, D) .* pop(1:count, :);
    next_fit = select .* offspring_fit + (1-select) .* fit(1:count);
    next_con = select .* offspring_con + (1-select) .* con(1:count);
end