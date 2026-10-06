
function generate_csv_summaries(base_folder, latest_folder, cfg)
    combined_rows = {};
    for e = 1:cfg.num_experiments
        exp_folder = fullfile(base_folder, sprintf('exp_%d', e));
        runs_folder = fullfile(exp_folder, 'runs');
        if ~exist(runs_folder, 'dir'); continue; end
        
        for ni = 1:numel(cfg.n_range)
            n = cfg.n_range(ni);
            for t = 1:numel(cfg.type_names)
                files = dir(fullfile(runs_folder, sprintf('n%d_type%d_run*.mat', n, t)));
                ee = []; re = []; ue = []; ok_count = 0; fail_count = 0;
                
                for f = 1:numel(files)
                    d = load(fullfile(runs_folder, files(f).name));
                    if d.failed
                        fail_count = fail_count + 1;
                    else
                        ok_count = ok_count + 1;
                        ee = [ee; d.errs];
                        re = [re; d.rel_errs];
                        ue = [ue; d.unit_errs];
                    end
                end
                
                if isempty(ee); ee = 0; end
                if isempty(re); re = 0; end
                if isempty(ue); ue = 0; end
                
                combined_rows(end+1, :) = { ...
                    e, cfg.type_names{t}, n, 2^n, ok_count, fail_count, ...
                    std(ee), median(ee), std(re), median(re), std(ue), median(ue) };
            end
        end
    end
    
    summary_table = cell2table(combined_rows, 'VariableNames', ...
        {'Experiment','Type','n','N','NumOK','NumFailed', ...
         'StdAbsErr','MedianAbsErr', 'StdRelErr','MedianRelErr', 'StdUnitErr','MedianUnitErr'});
     
    writetable(summary_table, fullfile(base_folder, 'combined_summary.csv'));
    writetable(summary_table, fullfile(latest_folder, 'combined_summary.csv'));
end
