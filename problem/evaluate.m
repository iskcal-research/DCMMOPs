function [fit, con, prob, evals] = evaluate(prob, pop, evals, is_evals)
    if is_evals && prob.rest == 0 % environment changes
        prob = change(prob);
        prob.rest = prob.freq;
    end
    dist = pdist2(pop, prob.X);
    fit = max(prob.h - prob.w .* dist, [], 2);
    
    cdist = pdist2(pop, prob.cX);
    con = prob.phi - max(prob.ch ./ (1 + prob.cw .* (cdist .^ 2)), [], 2);
    con = max(con, 0);
    
    if is_evals
        prob.rest = prob.rest - size(pop, 1);
    end
    evals = evals + size(pop, 1);
end