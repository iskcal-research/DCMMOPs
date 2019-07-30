function index = sort_indis(fit, con, cht)

    if cht.cmp == 1
        data = [con -fit];
        [~, index] = sortrows(data);
    elseif cht.cmp == 2
        new_con = con;
        new_con(new_con <= cht.ec.epsilon_g) = 0;
        data = [new_con -fit];
        [~, index] = sortrows(data);
    elseif cht.cmp == 3
        fit_con = fit - con .^2 .* cht.pf.r;
        [~, index] = sort(fit_con, 'descend');
    elseif cht.cmp == 4
        fit_con = fit - con .^2 .* cht.apf.v.* cht.apf.r;
        [~, index] = sort(fit_con, 'descend');        
    elseif cht.cmp == 5
        index = stochastic_sort(fit, con, cht.sr);
    else
        index = 1:length(fit);
    end
end