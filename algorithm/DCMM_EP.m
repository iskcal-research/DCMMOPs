function [peaks, all_ps] = DCMM_EP(D, G, L, P, R, cmp, run)
    prob = initialize_problem(D, G, L, P, R);    % initialize the information of the benchmark
    algRand = RandStream('mt19937ar','Seed', run); % the algorithm rand stream
    
    peaks = zeros(5, 60);
    all_ps = zeros(5, 60);
    
    evals = 0;
    NP = 100;                                   % population size
    phi = 2.0;                                  % the weight of NBC
    p_max = 10;                                 % the maximum size of species
    min_delete = 10;                            % the minimun delete number
    
    tau = 1/(sqrt(2*sqrt(D)));
    tau1 = 1/(sqrt(2*D));
    q = 10;
    
	maxG = prob.freq / NP;
    curG = 0;
    
    % memory
    mem.pop = [];                               % the memory population
    mem.dy = [];                                % the environment index
    mem.fit = [];                               % the memory fitness
    mem.con = [];                               % the memory constraints

    mem_maxsize = 50;                           % maximum size of the memory
    
    dy = 1;                                     % environment index
%     gen = 1;                                    % generation
        
    % initialize the population
    pop = rand(algRand, NP, D) * (prob.ub - prob.lb) + prob.lb;
    [fit, con, prob, evals] = evaluate(prob, pop, evals, true);
    eta = ones(size(pop)) .* 3;
    
    % calculate initial epsilon
    cht = struct();
    cht.cmp = cmp;
    
    % ec
    cht.ec.p = 0.5;
    sort_con = sort(con, 'ascend');
    theta = round(0.9 * NP);
    cht.ec.epsilon_0 = sort_con(theta);
    cht.ec.epsilon_g = 0;
    cht.ec.cp = -(log(cht.ec.epsilon_0)+5) / log(1-cht.ec.p);
    
    cht.pf.r = 10;
    cht.apf.r = 10;
    
    % sr
    cht.sr.algRand = algRand;
    cht.sr.pf = 0.45;
    
    while evals < prob.maxFEs
        
        % update the constraint handle technology
        if cmp == 2
            if curG / maxG < cht.ec.p
                cht.ec.epsilon_g = cht.ec.epsilon_0 * (1 - curG/maxG) ^ cht.ec.cp;
            else
                cht.ec.epsilon_g = 0;
            end
        elseif cmp == 4
            cht.apf.v = 1 - exp(-10*curG / maxG);
        end
        
        % sort the population from the best to the worst
        idx = sort_indis(fit, con, cht);
        fit = fit(idx, :);
        con = con(idx, :);
        pop = pop(idx, :);
        eta = eta(idx, :);

        % the distance between each pair of individuals
        matdis = pdist2(pop, pop);
        
        % NBC
        species = NBC(matdis, phi);
        
        % delete individuals in each species
        [species, dlist] = delete_redundancy(species, p_max, min_delete);

        evo_num = min(length(dlist), prob.rest);
        pop(dlist(1:evo_num), :) = rand(algRand, evo_num, D) .* (prob.ub - prob.lb) + prob.lb;
        [fit(dlist(1:evo_num)), con(dlist(1:evo_num)), prob, evals] = evaluate(prob, pop(dlist(1:evo_num), :), evals, true);
        eta(dlist(1:evo_num), :) = 3;
        
        if prob.rest > 0
            for i = 1:length(species)                
                sub_pop = pop(species(i).idx, :);                       
                sub_fit = fit(species(i).idx);
                sub_con = con(species(i).idx);
                sub_eta = eta(species(i).idx, :);
                
                [next_pop, next_fit, next_con, next_eta, evals, prob] = FEP(sub_pop, sub_fit, sub_con, sub_eta, tau, tau1, q, evals, prob, cht, algRand);
                pop(species(i).idx(1:length(next_fit)), :) = next_pop;
                fit(species(i).idx(1:length(next_fit))) = next_fit;
                con(species(i).idx(1:length(next_fit))) = next_con;
                eta(species(i).idx(1:length(next_fit)), :) = next_eta;
                
                if prob.rest == 0
                    break;
                end
            end
        end
        
        curG = curG + 1;
        
        % the environment changes
        if prob.rest == 0
            for i = 1:5
                [peak, all_p] = get_peaks_cons1(pop, fit, con, prob, 10^(-i));
                peaks(i, dy) = peak;
                all_ps(i, dy) = all_p;
            end
            [mem, pop, fit, con, prob, evals] = dyanmic_response(mem, pop, dy, prob, cat(2, species.seed), mem_maxsize, evals, cht, algRand);
            dy = dy + 1;
            curG = 0;
            eta = ones(size(pop)) .* 3;
            
            sort_con = sort(con, 'ascend');
            cht.ec.epsilon_0 = sort_con(theta);
            cht.ec.epsilon_g = 0;
            cht.ec.cp = -(log(cht.ec.epsilon_0)+5) / log(1-cht.ec.p);
        end
    end
end