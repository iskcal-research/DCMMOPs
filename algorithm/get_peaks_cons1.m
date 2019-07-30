function [peak, all_p] = get_peaks_cons1(pop, fit, con, prob, epsilon)
    [opt, opt_fit, p_idx] = optimal(prob);
    
    dist = pdist2(pop, prob.X(p_idx, :));
    [~, peak_idx] = max(prob.h(p_idx) - prob.w(p_idx) .* dist, [], 2);
    
    select_idx = (abs(opt_fit - fit) <= epsilon) & (con == 0);
    
    peak_idx = peak_idx(select_idx);
    
    peak = length(unique(peak_idx));
    all_p = size(opt, 1);
end