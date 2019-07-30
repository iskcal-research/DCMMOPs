function [mem, pop, fit, con, prob, evals] = dyanmic_response(mem, pop, dy, prob, seeds, mem_maxsize, evals, cht, algRand)
    [NP, D] = size(pop);
    num_select = min(5, length(seeds));
    
    % the individuals to be put into temp memory
    temp_M = pop(seeds(1:num_select), :);
    
    % update the population
    [fit, con, prob, evals] = evaluate(prob, pop, evals, true);
    temp_newfit = fit(seeds(1:num_select));
    temp_newcon = con(seeds(1:num_select));
    sort_index = sort_indis(fit, con, cht);    
%     [fit, sort_index] = sort(fit, 'descend');
    pop = pop(sort_index, :);
    fit = fit(sort_index);
    con = con(sort_index);
    
    % update the memory's fitness and add the individuals to the population
    s_index = floor(NP/2);
    count = 0;
    if ~isempty(mem.pop)
        [mem.fit, mem.con, prob, evals] = evaluate(prob, mem.pop, evals, true);
        sort_idx = sort_indis(mem.fit, mem.con, cht);
        mem.pop = mem.pop(sort_idx, :);
        mem.fit = mem.fit(sort_idx);
        mem.con = mem.con(sort_idx);
        mem.dy = mem.dy(sort_idx);
        
%         [mem.fit, sort_index] = sort(mem.fit, 'descend');
%         mem.pop = mem.pop(sort_index, :);
%         mem.dy = mem.dy(sort_index);
        
        best_index = find(mem.dy == mem.dy(1));
        count = min(length(best_index), NP - s_index);
        pop(s_index+1:s_index+count, :) = mem.pop(best_index(1:count), :);
        fit(s_index+1:s_index+count) = mem.fit(best_index(1:count), :);
        con(s_index+1:s_index+count) = mem.con(best_index(1:count), :);
    end
    
    pop(s_index+count+1:NP, :) = rand(algRand, NP - s_index - count, D) * (prob.ub - prob.lb) + prob.lb;
    [fit(s_index+count+1:NP), con(s_index+count+1:NP), prob, evals] = evaluate(prob, pop(s_index+count+1:NP, :), evals, true);
        
    % update memory
    for i = num_select:-1:1
        mem_size = size(mem.pop, 1);
        if mem_size < mem_maxsize
            mem.pop = [mem.pop; temp_M(i, :)];
            mem.fit = [mem.fit; temp_newfit(i)];
            mem.con = [mem.con; temp_newcon(i)];
            mem.dy = [mem.dy; dy];
        else
            [~, mn] = min(pdist2(temp_M(i, :), mem.pop));
            mem.pop(mn, :) = temp_M(i, :);
            mem.fit(mn) = temp_newfit(i);
            mem.con(mn) = temp_newcon(i);
            mem.dy(mn) = dy;
        end
    end
end