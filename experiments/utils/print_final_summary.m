function print_final_summary(base_folder,cfg)

fprintf('\n');
fprintf('===============================================================\n');
fprintf('                    FINAL BENCHMARK SUMMARY\n');
fprintf('===============================================================\n');

n_types = numel(cfg.type_names);

for t = 1:n_types

    fprintf('\nGenerator: %s\n',cfg.type_names{t});
    fprintf('---------------------------------------------------------------------------------------------------------------\n');
    fprintf('%6s %6s %6s %14s %14s %14s %14s %14s\n', ...
        'N','OK','Fail','MeanErr','StdErr','MedianErr','MeanRelErr','MeanUnitErr');
    fprintf('---------------------------------------------------------------------------------------------------------------\n');

    for ni = 1:numel(cfg.n_range)

        n = cfg.n_range(ni);
        N = 2^n;

        absErr = [];
        relErr = [];
        unitErr = [];

        ok = 0;
        fail = 0;

        fprintf('\nN = %d\n',N);
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('%4s %4s %14s %14s %18s %8s\n',...
    'Exp','Run','Abs Error','Rel Error','Unitarity','Status');
fprintf('-----------------------------------------------------------------------------------------------\n');

for e = 1:cfg.num_experiments

    runs_folder = fullfile(base_folder,...
        sprintf('exp_%d',e),'runs');

    files = dir(fullfile(runs_folder,...
        sprintf('n%d_type%d_run*.mat',n,t)));

    for f = 1:numel(files)

        d = load(fullfile(runs_folder,files(f).name));

        if d.failed
            fprintf('%4d %4d %14s %14s %18s %8s\n',...
                e,d.r_save,'-','-','-','FAILED');
        else
            fprintf('%4d %4d %14.3e %14.3e %18.3e %8s\n',...
                e,...
                d.r_save,...
                d.errs,...
                d.rel_errs,...
                d.unit_errs,...
                'OK');
        end

    end
end

    fprintf('---------------------------------------------------------------------------------------------------------------\n');

end

fprintf('\n');

end

end
