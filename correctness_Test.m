%RUN_ALL_SEMISEPARABLE_BENCHMARKS_V2
% Benchmark script with:
%   - progress bar
%   - 4 independent full repetitions of the experiment (no averaging
%     across repetitions -- each is saved and plotted separately)
%   - IMMEDIATE saving after every single run (so killing the script
%     mid-way still leaves all completed runs on disk) -- toggle with
%     cfg.enable_saving
%   - saves u, v, S, S_tilde for every run (the full circuit unitary U
%     is NOT saved -- only the scalar unitarity error derived from it)
%   - per-run diagnostics: absolute error, relative error, and
%     (optionally) the unitarity error ||U'*U - I||_2 -- toggle with
%     cfg.enable_unitarity_check
%   - summary statistics (std/median -- NO mean) exported to CSV
%   - robust to per-run circuit-build failures (recorded, not fatal)
%   - overall live error plot, updated and saved after every single
%     run, one legend entry per generator type; the unitarity-error
%     plot is only created/updated if cfg.enable_unitarity_check is true
%
% Computes semiseparable block-encoding errors.

clear;
clc;

%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.n_range               = 2:5;
cfg.num_runs               = 1;                    % runs per (experiment, n, type)
cfg.num_experiments        = 4;                     % full independent repetitions
cfg.type_names              = {'Random','Tridiagonal','Exponential'};
cfg.save_formats            = {'png','eps','fig'};  % remove entries to save disk/time
cfg.use_parfor              = false;                % set true if Parallel Computing Toolbox is available
cfg.enable_saving           = true;                 % master switch: save anything to disk (mat/csv/figures)?
cfg.enable_unitarity_check  = false;                 % compute ||U'*U - I||_2 and its plot?

n_types = numel(cfg.type_names);
colors  = lines(n_types);
markers = {'o','s','d','^','v','p','h'};    % extend automatically if n_types > 3
markers = markers(mod(0:n_types-1, numel(markers)) + 1);

% distinct colors for the 4 experiment repetitions (used in combined plots)
exp_colors = lines(cfg.num_experiments);

num_n       = length(cfg.n_range);
total_tasks = num_n * n_types * cfg.num_runs * cfg.num_experiments;
done_tasks  = 0;

%% ============================================================
% FOLDERS (created immediately, before any computation)
%% ============================================================

timestamp     = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
base_folder   = fullfile('results', timestamp);
latest_folder = fullfile('results', 'latest');
plots_folder  = fullfile(base_folder, 'plots');

if cfg.enable_saving
    if ~exist(base_folder, 'dir');   mkdir(base_folder);   end
    if ~exist(latest_folder, 'dir'); mkdir(latest_folder); end
    if ~exist(plots_folder, 'dir');  mkdir(plots_folder);  end
end

%% ============================================================
% STORAGE (across all experiments, types, n)
%   Indexing: all_*{e}{t, ni}
%% ============================================================

all_errors      = cell(cfg.num_experiments, 1);
all_rel_errors  = cell(cfg.num_experiments, 1);
all_unit_errors = cell(cfg.num_experiments, 1);
all_failed      = cell(cfg.num_experiments, 1);
all_u           = cell(cfg.num_experiments, 1);
all_v           = cell(cfg.num_experiments, 1);
all_S           = cell(cfg.num_experiments, 1);
all_Stilde      = cell(cfg.num_experiments, 1);

for e = 1:cfg.num_experiments
    all_errors{e}      = cell(n_types, num_n);
    all_rel_errors{e}  = cell(n_types, num_n);
    all_unit_errors{e} = cell(n_types, num_n);
    all_failed{e}       = cell(n_types, num_n);
    all_u{e}            = cell(n_types, num_n);
    all_v{e}            = cell(n_types, num_n);
    all_S{e}            = cell(n_types, num_n);
    all_Stilde{e}       = cell(n_types, num_n);
end

%% ============================================================
% LIVE PLOTS -- created once, updated + saved after EVERY run
%% ============================================================

fig_err = figure('Name', 'Live: Semiseparable Error', 'NumberTitle', 'off');
ax_err = axes(fig_err); hold(ax_err, 'on');

h_err = gobjects(n_types, 1);
for t = 1:n_types
    h_err(t) = scatter(ax_err, NaN, NaN, 80, ...
        'Marker', markers{t}, 'MarkerEdgeColor', colors(t, :), ...
        'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
        'DisplayName', cfg.type_names{t});
end

xlabel(ax_err, 'n'); ylabel(ax_err, '||S - \tilde{S}||_2');
title(ax_err, 'Semiseparable Error -- All Generators, All Experiments (live)');
grid(ax_err, 'on'); legend(ax_err, h_err, 'Location', 'eastoutside');

if cfg.enable_unitarity_check
    fig_unit = figure('Name', 'Live: Unitarity Error', 'NumberTitle', 'off');
    ax_unit = axes(fig_unit); hold(ax_unit, 'on');

    h_unit = gobjects(n_types, 1);
    for t = 1:n_types
        h_unit(t) = scatter(ax_unit, NaN, NaN, 80, ...
            'Marker', markers{t}, 'MarkerEdgeColor', colors(t, :), ...
            'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
            'DisplayName', cfg.type_names{t});
    end

    xlabel(ax_unit, 'n'); ylabel(ax_unit, '||U^\dagger U - I||_2');
    title(ax_unit, 'Unitarity Error -- All Generators, All Experiments (live)');
    grid(ax_unit, 'on'); legend(ax_unit, h_unit, 'Location', 'eastoutside');
end

drawnow;

%% ============================================================
% START TIMER
%% ============================================================

tic;


for ni = 1:num_n


for e = 1:cfg.num_experiments

    exp_folder = fullfile(base_folder, sprintf('exp_%d', e));
    runs_folder = fullfile(exp_folder, 'runs');
    if cfg.enable_saving
        if ~exist(exp_folder, 'dir'); mkdir(exp_folder); end
        if ~exist(runs_folder, 'dir'); mkdir(runs_folder); end
    end



        n = cfg.n_range(ni);
        N = 2^n;

        for t = 1:n_types

            errs       = zeros(cfg.num_runs, 1);
            rel_errs   = zeros(cfg.num_runs, 1);
            unit_errs  = zeros(cfg.num_runs, 1);
            failed     = false(cfg.num_runs, 1);
            u_cell     = cell(cfg.num_runs, 1);
            v_cell     = cell(cfg.num_runs, 1);
            S_cell     = cell(cfg.num_runs, 1);
            Stilde_cell = cell(cfg.num_runs, 1);

            for r = 1:cfg.num_runs

                rng(100000*e + 1000*n + 10*t + r);

                [u, v] = local_generate_uv(n, t);
                S      = tril(u*v') + triu(v*u', 1);
                U      = [];
                S_tilde = [];
                unit_err = NaN;

                try
                    circ = build_semiseparable_circuit(u, v);
                    U       = circ.matrix;
                    S_tilde = 2*sqrt(N) * U(1:N, 1:N);

                    errs(r)     = norm(S_tilde - S, 2);
                    normS       = norm(S, 2);
                    rel_errs(r) = errs(r) / max(normS, eps);

                    if cfg.enable_unitarity_check
                        Iu = eye(size(U, 1));
                        unit_err     = norm(U' * U - Iu, 2);
                        unit_errs(r) = unit_err;
                    end

                catch ME
                    warning('Run failed (exp=%d, n=%d, type=%d, run=%d): %s', ...
                        e, n, t, r, ME.message);
                    errs(r)      = NaN;
                    rel_errs(r)  = NaN;
                    unit_errs(r) = NaN;
                    failed(r)    = true;
                end

                u_cell{r}      = u;
                v_cell{r}      = v;
                S_cell{r}      = S;
                Stilde_cell{r} = S_tilde;

                % ---------------------------------------------------
                % SAVE IMMEDIATELY: this single run, right now.
                % ---------------------------------------------------
                run_abs_err   = errs(r);     
                run_rel_err   = rel_errs(r); 
                run_unit_err  = unit_errs(r);
                run_failed    = failed(r);   
                if cfg.enable_saving
                    run_file = fullfile(runs_folder, ...
                        sprintf('n%d_type%d_run%d.mat', n, t, r));
                    save(run_file, 'u', 'v', 'S', 'S_tilde', ...
                        'run_abs_err', 'run_rel_err', 'run_unit_err', ...
                        'run_failed', 'n', 't', 'r', 'e');
                end

                % ---------------------------------------------------
                % LIVE PLOT UPDATE: add this run's point right now,
                % redraw on screen, and re-save the overall plots so
                % the files on disk always reflect the latest state.
                % ---------------------------------------------------
                if ~failed(r)
                    scatter(ax_err, n, errs(r), 80, ...
                        'Marker', markers{t}, 'MarkerEdgeColor', colors(t, :), ...
                        'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
                        'HandleVisibility', 'off');
                    if cfg.enable_unitarity_check
                        scatter(ax_unit, n, unit_errs(r), 80, ...
                            'Marker', markers{t}, 'MarkerEdgeColor', colors(t, :), ...
                            'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
                            'HandleVisibility', 'off');
                    end
                    drawnow;
                end

                if cfg.enable_saving
                    save_figure_pair(fig_err, plots_folder, latest_folder, ...
                        'scatter_overlay_all', cfg.save_formats);
                    if cfg.enable_unitarity_check
                        save_figure_pair(fig_unit, plots_folder, latest_folder, ...
                            'scatter_unitarity_overall', cfg.save_formats);
                    end
                end

                done_tasks = done_tasks + 1;
                print_progress(done_tasks, total_tasks, e, n, t, r, cfg.num_runs);
            end

            all_errors{e}{t, ni}      = errs;
            all_rel_errors{e}{t, ni}  = rel_errs;
            all_unit_errors{e}{t, ni} = unit_errs;
            all_failed{e}{t, ni}      = failed;
            all_u{e}{t, ni}           = u_cell;
            all_v{e}{t, ni}           = v_cell;
            all_S{e}{t, ni}           = S_cell;
            all_Stilde{e}{t, ni}      = Stilde_cell;

            % ---------------------------------------------------
            % SAVE IMMEDIATELY: running snapshot of everything so
            % far in this experiment, overwritten after every
            % (n, type) block completes.
            % ---------------------------------------------------
            if cfg.enable_saving
                partial_errors      = all_errors{e};      
                partial_rel_errors  = all_rel_errors{e};   
                partial_unit_errors = all_unit_errors{e};  
                partial_failed      = all_failed{e};       
                partial_u           = all_u{e};            
                partial_v           = all_v{e};             
                partial_S           = all_S{e};             
                partial_Stilde      = all_Stilde{e};         
                save(fullfile(exp_folder, 'partial_results.mat'), ...
                    'partial_errors', 'partial_rel_errors', 'partial_unit_errors', ...
                    'partial_failed', 'partial_u', 'partial_v', 'partial_S', ...
                    'partial_Stilde', 'cfg', 'timestamp', 'e');
            end
        end
    end

    % -------------------------------------------------------------
    % SAVE IMMEDIATELY: full results for this completed experiment.
    % -------------------------------------------------------------
    if cfg.enable_saving
        exp_errors      = all_errors{e};      
        exp_rel_errors  = all_rel_errors{e};   
        exp_unit_errors = all_unit_errors{e};  
        exp_failed      = all_failed{e};       
        exp_u           = all_u{e};            
        exp_v           = all_v{e};             
        exp_S           = all_S{e};             
        exp_Stilde      = all_Stilde{e};         
        save(fullfile(exp_folder, 'results.mat'), ...
            'exp_errors', 'exp_rel_errors', 'exp_unit_errors', 'exp_failed', ...
            'exp_u', 'exp_v', 'exp_S', 'exp_Stilde', 'cfg', 'timestamp', 'e');

        % per-experiment summary CSV (no mean columns)
        exp_summary_rows = {};
        for t = 1:n_types
            for ni = 1:num_n
                n  = cfg.n_range(ni);
                ee = all_errors{e}{t, ni};
                re = all_rel_errors{e}{t, ni};
                ue = all_unit_errors{e}{t, ni};
                ok = ~all_failed{e}{t, ni};

                exp_summary_rows(end+1, :) = { ...
                    e, cfg.type_names{t}, n, 2^n, ...
                    sum(ok), sum(~ok), ...
                    std(ee(ok)),  median(ee(ok)), ...
                    std(re(ok)),  median(re(ok)), ...
                    std(ue(ok)),  median(ue(ok)) }; 
            end
        end
        exp_summary_table = cell2table(exp_summary_rows, 'VariableNames', ...
            {'Experiment','Type','n','N','NumOK','NumFailed', ...
             'StdAbsErr','MedianAbsErr', ...
             'StdRelErr','MedianRelErr', ...
             'StdUnitErr','MedianUnitErr'});
        writetable(exp_summary_table, fullfile(exp_folder, 'summary.csv'));
    end
end

fprintf('\n\nDone.\n');

%% ============================================================
% COMBINED SUMMARY (all experiments, NO mean column) -> CSV
%% ============================================================

combined_rows = {};
for e = 1:cfg.num_experiments
    for t = 1:n_types
        for ni = 1:num_n
            n  = cfg.n_range(ni);
            ee = all_errors{e}{t, ni};
            re = all_rel_errors{e}{t, ni};
            ue = all_unit_errors{e}{t, ni};
            ok = ~all_failed{e}{t, ni};

            combined_rows(end+1, :) = { ...
                e, cfg.type_names{t}, n, 2^n, ...
                sum(ok), sum(~ok), ...
                std(ee(ok)),  median(ee(ok)), ...
                std(re(ok)),  median(re(ok)), ...
                std(ue(ok)),  median(ue(ok)) }; 
        end
    end
end

combined_summary_table = cell2table(combined_rows, 'VariableNames', ...
    {'Experiment','Type','n','N','NumOK','NumFailed', ...
     'StdAbsErr','MedianAbsErr', ...
     'StdRelErr','MedianRelErr', ...
     'StdUnitErr','MedianUnitErr'});

if cfg.enable_saving
    writetable(combined_summary_table, fullfile(base_folder, 'combined_summary.csv'));
    writetable(combined_summary_table, fullfile(latest_folder, 'combined_summary.csv'));
end

disp(combined_summary_table);

%% ============================================================
% PLOTTING already done live during the run -- fig_err (and fig_unit,
% if cfg.enable_unitarity_check was true) hold the final state and
% were re-saved after every single point (if cfg.enable_saving was
% true). Left open on screen; close manually or programmatically as
% desired.
%% ============================================================

%% ============================================================
% FINAL SAVE (all experiments together, mirrored to latest)
%% ============================================================

if cfg.enable_saving
    save(fullfile(base_folder, 'all_results.mat'), ...
        'all_errors', 'all_rel_errors', 'all_unit_errors', 'all_failed', ...
        'all_u', 'all_v', 'all_S', 'all_Stilde', 'cfg', 'timestamp');
    save(fullfile(latest_folder, 'all_results.mat'), ...
        'all_errors', 'all_rel_errors', 'all_unit_errors', 'all_failed', ...
        'all_u', 'all_v', 'all_S', 'all_Stilde', 'cfg', 'timestamp');
end

%% ============================================================
% END TIMER
%% ============================================================

elapsed = toc;
fprintf('\nExperiment completed in %.2f seconds (%.2f min)\n', elapsed, elapsed/60);
fprintf('Timestamp: %s\n', timestamp);
if cfg.enable_saving
    fprintf('Results folder: %s\n', base_folder);
else
    fprintf('Saving was disabled (cfg.enable_saving = false) -- nothing written to disk.\n');
end

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function print_progress(done, total, e, n, t, r, num_runs)
    percent  = 100 * done / total;
    bar_len  = 30;
    filled   = round(bar_len * done / total);
    bar_str  = [repmat('#', 1, filled), repmat('-', 1, bar_len - filled)];

    fprintf('\r[%s] %6.2f%% | exp=%d | n=%d | type=%d | run=%d/%d', ...
        bar_str, percent, e, n, t, r, num_runs);

    if done == total
        fprintf('\n');
    end
end

function save_figure_pair(fig, base_folder, latest_folder, name, formats)
    for i = 1:numel(formats)
        fmt = formats{i};
        switch fmt
            case 'fig'
                savefig(fig, fullfile(base_folder, [name '.fig']));
                savefig(fig, fullfile(latest_folder, [name '.fig']));
            otherwise
                saveas(fig, fullfile(base_folder, [name '.' fmt]));
                saveas(fig, fullfile(latest_folder, [name '.' fmt]));
        end
    end
end