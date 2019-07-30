function prob = initialize_problem(D, G, L, P, R)
    if P > G
        R = R + (P - G);
        P = G;
    end

    prob.proRand = RandStream('mt19937ar','Seed', 0);    % rand stream
    prob.D = D;                                          % dimension (D)
    prob.G = G;                                          % the No. of the global peaks (G)
    prob.L = L;                                          % the No. of the local peaks (L)
    prob.num_peaks = prob.G + prob.L;                    % the total No. of peaks (G+L)
    prob.s = 0.1;                                        % shift length (s)
    prob.lb = -5;                                        % lower boundary (lb)
    prob.ub = 5;                                         % upper boundary (ub)
    prob.num_change = 60;                                % the number of the environment change (Dy)
    prob.freq = 3000 * D;                                % the maximum evoluations of an environment (freq)
    prob.maxFEs = prob.num_change * prob.freq;           % the total evolutions (MaxFEs)
    prob.dpeaks = 0.1;                                   % the minimum distance between the peaks (dpeaks)
    prob.height = 50;                                    % initial height (height)
    prob.lh_min = 5;                                     % minimal height difference (dif_{min})
    prob.lh_max = 20;                                    % maximal height difference (dif_{max})
    prob.h_min = 30;                                     % minimal height (h_{min})
    prob.h_max = 70;                                     % maximal height (h_{max})
    prob.h_s = 7.0;                                      % change severity of height (h_s)
    prob.w_min = 1;                                      % minimal width (w_{min})
    prob.w_max = 12;                                     % maximal width (w_{max})
    prob.w_s = 1.0;                                      % change severity of width (w_s)
    
    prob.P = P;                                          % the No. of first type of feasible regions (P)
    prob.R = R;                                          % the No. of second type of feasible regions (R)
    prob.cnum_peaks = prob.P + prob.R;                   % the total No. of feasible regions (P+R)
    prob.cheight = 50;                                   % initial constraint height (height^c)
    prob.ch_min = 30;                                    % minimal constraint height (h_{min}^c)
    prob.ch_max = 70;                                    % maximal constraint height (h_{max}^c)
	prob.ch_s = 7.0;                                     % change severity of constraint height (h_s^c)
    prob.cw_min = 1;                                     % minimal constraint width (w_{min}^c)
    prob.cw_max = 12;                                    % maximal constraint width (w_{max}^c)
    prob.cw_s = 1.0;                                     % change severity of constraint width (w_s^c)
    prob.cv = 6;                                         % constraint control value (cv)
    
    prob.rest = prob.freq;                               % the rest evaluations in the current envirnment
    
    % initialize the position, height and width
    prob.X = zeros(prob.num_peaks, prob.D);              % positions of the fitness peaks
    for i = 1:prob.num_peaks
       accept = false;
       while ~accept
           prob.X(i, :) = rand(prob.proRand, 1, prob.D) * (prob.ub - prob.lb) + prob.lb;
           if i == 1
               accept = true;
           else
               dis = pdist2(prob.X(i, :), prob.X(1:i-1, :));
               if all(dis(:) > prob.dpeaks)
                   accept = true;
               end
           end
       end
    end
    prob.h = ones(1, prob.num_peaks) * prob.height - [zeros(1, prob.G), rand(prob.proRand, 1, prob.L) * (prob.lh_max-prob.lh_min) + prob.lh_min];
    prob.w = rand(prob.proRand, 1, prob.num_peaks) * (prob.w_max-prob.w_min)+prob.w_min;
    
    prob.ch = ones(1, prob.cnum_peaks) * prob.cheight;
    prob.cw = rand(prob.proRand, 1, prob.cnum_peaks) * (prob.cw_max-prob.cw_min)+prob.cw_min;
    prob.phi = min(prob.ch) / prob.cv;                  % the parameter \delta(t) in the constraint landscape
    prob.cr = sqrt((prob.ch ./ prob.phi - 1) ./ prob.cw);  % the radius of feasible regions
    
    % generate the constraint peaks location
    select_peaks = randperm(prob.proRand, G, P);
    select_lb = prob.X(select_peaks, :) - repmat(prob.cr(1:prob.P)', 1, prob.D) ./ sqrt(prob.D);
    select_ub = prob.X(select_peaks, :) + repmat(prob.cr(1:prob.P)', 1, prob.D) ./ sqrt(prob.D);
    
    prob.cX = rand(prob.proRand, prob.P, prob.D) .* (select_ub - select_lb) + select_lb;
    prob.cX = [prob.cX; rand(prob.proRand, prob.R, prob.D) .* (prob.ub - prob.lb) + prob.lb];
    
end