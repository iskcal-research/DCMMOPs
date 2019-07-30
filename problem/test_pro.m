function test_pro(prob)
    x = [prob.lb:(prob.ub-prob.lb)/100:prob.ub];
    y = [prob.lb:(prob.ub-prob.lb)/100:prob.ub];
    u = [];
    for i = 1:length(x)
        for j = 1:length(y)
            u = [u; x(i), y(j)];
        end
    end
    
    [fit, con, ~,~] = evaluate(prob, u, 1, false);
    feasible = u(con == 0, :);
    plot(feasible(:, 1), feasible(:, 2), 'ro');
    hold on;
    
    plot(prob.X(1:prob.G, 1), prob.X(1:prob.G, 2), 'bo', 'MarkerFaceColor', 'b');
    plot(prob.X(prob.G+1:prob.num_peaks, 1), prob.X(prob.G+1:prob.num_peaks, 2), 'go', 'MarkerFaceColor', 'g');
    
    opt = optimal(prob);
    plot(opt(:, 1), opt(:, 2), 'cs', 'MarkerFaceColor', 'c');
    hold off;
end