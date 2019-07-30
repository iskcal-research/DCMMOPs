function [opt, opt_fit, p_idx] = optimal(prob)
    pc = prob.X;
    
    [fit, con, ~, ~] = evaluate(prob, pc, 0, false);
    
    opt_fit = max(fit);
    p_idx = find((fit == opt_fit & (con == 0)));

    opt = pc(p_idx, :);
end