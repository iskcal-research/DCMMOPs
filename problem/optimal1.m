function [opt, opt_fit, p_idx] = optimal1(prob)
    pc = prob.X;
    rc = prob.cX;
%     r = (prob.ch - prob.phi) ./ prob.cw;
    r = sqrt((prob.ch ./ prob.phi - 1) ./ prob.cw);
    
    local = [];
    fea_p_idx = [];
    for i = 1:size(rc, 1)
        for j = 1:size(pc, 1)
            dist = pdist2(rc(i, :), pc(j, :));
            if dist <= r(i)
                local = [local; pc(j, :)];
            else
                cur_local = (1 - r(i) / dist) .* rc(i, :) + r(i) / dist .* pc(j, :); 
                local = [local; cur_local];
            end
            fea_p_idx = [fea_p_idx; j];
        end
    end
    
    [local, se_idx, ~] = unique(local, 'rows', 'stable');
    fea_p_idx = fea_p_idx(se_idx);
    
    [fit, con, ~, ~] = evaluate(prob, local, 0, false);
%     best_fit = max(fit(con == 0));
%     best_idx = (fit == best_fit) & (con == 0);
    best_fit = max(fit);
    best_idx = fit == best_fit;

    opt = local(best_idx, :);
    p_idx = fea_p_idx(best_idx);
    opt_fit = best_fit;
end