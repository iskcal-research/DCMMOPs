function prob = change(prob)
    % update the position
    for i = 1:prob.num_peaks
        accept = false;
        while ~accept
            v = rand(prob.proRand, 1, prob.D) - 0.5;
            len = sqrt(sum(v.^2, 2));
            if len ~= 0
                len = prob.s / len;
            end
            v = v * len;
            prob.X(i,:) = prob.X(i, :) + v;
            prob.X(i,:) = boundary_check(prob.X(i,:), prob.lb, prob.ub);
            
            % check the distance
            if i == 1 
                accept = true;
            else
                dis = pdist2(prob.X(i, :), prob.X(1:i-1, :));
                if all(dis(:) > prob.dpeaks) % accept the position
                    accept = true;
                end
            end
        end
    end
    % update the height
    prob.height = prob.height + randn(prob.proRand, 1, 1) * prob.h_s;
    prob.height = max(prob.height, prob.h_min+prob.lh_min);
    prob.height = min(prob.height, prob.h_max);
    
    prob.h = ones(1, prob.num_peaks) * prob.height - [zeros(1, prob.G), rand(prob.proRand, 1, prob.L) * (prob.lh_max-prob.lh_min) + prob.lh_min];
    prob.h = max(prob.h, prob.h_min);
    prob.h = min(prob.h, prob.h_max);
    
    % update the width
    prob.w = prob.w + randn(prob.proRand, 1, prob.num_peaks) * prob.w_s;
    prob.w = max(prob.w, prob.w_min);
    prob.w = min(prob.w, prob.w_max);
    
    % update the constraints
    prob.ch = prob.ch + randn(prob.proRand, 1, prob.cnum_peaks) * prob.ch_s;
    prob.ch = max(prob.ch, prob.ch_min);
    prob.ch = min(prob.ch, prob.ch_max);
    prob.phi = min(prob.ch)/prob.cv;
    prob.cw = prob.cw + randn(prob.proRand, 1, prob.cnum_peaks) * prob.cw_s;
    prob.cw = max(prob.cw, prob.cw_min);
    prob.cw = min(prob.cw, prob.cw_max);
    prob.cr = sqrt((prob.ch ./ prob.phi - 1) ./ prob.cw); 
    
    % generate the constraint peaks location
    select_peaks = randperm(prob.proRand, prob.G, prob.P);
    select_lb = prob.X(select_peaks, :) - repmat(prob.cr(1:prob.P)', 1, prob.D) ./ sqrt(prob.D);
    select_ub = prob.X(select_peaks, :) + repmat(prob.cr(1:prob.P)', 1, prob.D) ./ sqrt(prob.D);
    
    prob.cX = rand(prob.proRand, prob.P, prob.D) .* (select_ub - select_lb) + select_lb;
    prob.cX = [prob.cX; rand(prob.proRand, prob.R, prob.D) .* (prob.ub - prob.lb) + prob.lb];
end
