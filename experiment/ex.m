function ex()
    max_run = 30;
    
    DArr = [2 5];
    GArr = [2 5 10];
    PArr = [5 10];
    cmpArr = 1:5;
    algArr = 1:4;
    for cmp = cmpArr
        for alg = algArr
            for D = DArr
                for G = GArr
                    for P = PArr
%                         labindex = 4;
                        disp(['D:',num2str(D), ' | G:', num2str(G), ' | P:', num2str(P), ' | alg:', num2str(alg), ' | cmp:', num2str(cmp)]);
                        delete(gcp('nocreate'));
                        parpool('local',max_run);
                        spmd(max_run)
                            switch(alg)
                                case 1
                                    [peaks, all_p] = DCMM_CSA(D, G, 2, P, 2, cmp, labindex);
                                case 2
                                    [peaks, all_p] = DCMM_DE(D, G, 2, P, 2, cmp, labindex);
                                case 3
                                    [peaks, all_p] = DCMM_PSO(D, G, 2, P, 2, cmp, labindex);
                                case 4
                                    [peaks, all_p] = DCMM_EP(D, G, 2, P, 2, cmp, labindex);
                            end
                        end
%                         result = cat(1, result{1:end});

                        peaks = cat(1, peaks{1:end});
                        all_p = cat(1, all_p{1:end});

                        acc = zeros(5, 1);
                        for i = 1:5
                            e_peaks = peaks(i:5:end, :);
                            e_all = all_p(i:5:end, :);
                            acc(i) = sum(e_peaks(:)) / sum(e_all(:));
                        end
                        dlmwrite(sprintf('./result/data/A%d+D%d+G%d+P%d+C%d', alg, D, G, P, cmp), acc);
                        dlmwrite(sprintf('./result/Details/A%d+D%d+G%d+P%d+C%d', alg, D, G, P, cmp), [peaks; all_p]);
                    end
                end
            end
        end
    end
end
