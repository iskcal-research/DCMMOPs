function [species, dlist] = delete_redundancy(species, p_max, min_delete)
    cur_d = 0;  % current delete number
    dlist = []; % the index of delete
    % the size of the species can not exceed 'p_max'
    for i = 1:length(species)
       if species(i).len > p_max 
           cur_d = cur_d + species(i).len - p_max;
           dlist = [dlist; species(i).idx(p_max+1:end)];
           species(i).idx = species(i).idx(1:p_max);
           species(i).len = p_max;
       end
    end
    % ensure minimum number of deletions 
    i = length(species);
    while cur_d < min_delete
        cur_d = cur_d + 1;
        dlist = [dlist; species(i).idx(end)];
        species(i).idx(end) = [];
        species(i).len = species(i).len - 1;
        if species(i).len == 0
            species(i) = [];
        end
        i = i - 1;
        if i == 0
            i = length(species);
        end
    end
end