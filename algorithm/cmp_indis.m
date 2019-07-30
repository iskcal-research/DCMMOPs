function select = cmp_indis(pop_fit, off_fit, pop_con, off_con, cht)
    cmp_fit = off_fit > pop_fit;
	cmp_con = off_con < pop_con;

    if cht.cmp == 1
        cmp_fit_index = (pop_con == 0 & off_con == 0) | (pop_con == off_con);
        cmp_con_index = ~cmp_fit_index;
        select = cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con;
    elseif cht.cmp == 2
        cmp_fit_index = (pop_con <= cht.ec.epsilon_g & off_con <= cht.ec.epsilon_g) | (pop_con == off_con);
        cmp_con_index = ~cmp_fit_index;
        select = cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con;
    elseif cht.cmp == 3
        off_fit_con = off_fit - off_con .^2 .* cht.pf.r;
        pop_fit_con = pop_fit - pop_con .^2 .* cht.pf.r;
        select = off_fit_con > pop_fit_con;
    elseif cht.cmp == 4
        % Varying Fitness Functions in Genetic Algorithms: Studying the Rate of Increase of the Dynamic Penalty Terms
        
        off_fit_con = off_fit - off_con .^2 .* cht.apf.v.* cht.apf.r;
        pop_fit_con = pop_fit - pop_con .^2 .* cht.apf.v.* cht.apf.r;
        select = off_fit_con > pop_fit_con;     
    elseif cht.cmp == 5
        cmp_fit_index = pop_con == 0 & off_con == 0;
        cmp_con_count = sum(~cmp_fit_index);
        cmp_fit_index(~cmp_fit_index) = rand(cht.sr.algRand, cmp_con_count, 1) < cht.sr.pf;
        cmp_con_index = ~cmp_fit_index;
        select = cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con;
    else
        select = cmp_fit;
    end
end
