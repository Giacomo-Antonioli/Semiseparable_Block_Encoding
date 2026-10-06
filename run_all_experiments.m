%RUN_ALL_SEMISEPARABLE_BENCHMARKS_IMPROVED
% Benchmark script with:
%   - progress bar
%   - full experiment
%   - mirrored plots (timestamp + latest)
%   - per-run diagnostics (absolute + relative error)
%   - summary statistics (mean/std/median) exported to CSV
%   - robust to per-run circuit-build failures (recorded, not fatal)
%
% Computes semiseparable block-encoding errors.

clear;
clc;

%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.n_range      = 2:3;
cfg.num_runs     = 1;
cfg.type_names   = {'Random','Tridiagonal','Exponential'};
cfg.save_formats = {'png','eps','fig'};   % remove entries to save disk/time
cfg.use_parfor   = false;                  % set true if Parallel Computing Toolbox is available

n_types = numel(cfg.type_names);
colors  = lines(n_types);
markers = {'o','s','d','^','v','p','h'};   % extend automatically if n_types > 3
markers = markers(mod(0:n_types-1, numel(markers)) + 1);

%% ============================================================
% TOTAL WORK
%% ============================================================

num_n       = length(cfg.n_range);
total_tasks = num_n * n_types * cfg.num_runs;
done_tasks  = 0;

%% ============================================================
% STORAGE
%% ============================================================

all_errors     = cell(n_types, num_n);   % ||S - S_tilde||_2
all_rel_errors = cell(n_types, num_n);   % ||S - S_tilde||_2 / ||S||_2
all_failed     = cell(n_types, num_n);   % logical mask of failed runs

%% ============================================================
% START TIMER
%% ============================================================

tic;

%% ============================================================
% MAIN LOOP
%% ============================================================

for ni = 1:num_n

    n = cfg.n_range(ni);
    N = 2^n;

    for t = 1:n_types

        errs      = zeros(cfg.num_runs, 1);
        rel_errs  = zeros(cfg.num_runs, 1);
        failed    = false(cfg.num_runs, 1);

        for r = 1:cfg.num_runs

            rng(1000*n + 10*t + r);

            [u, v] = local_generate_uv(n, t);

            try
                circ = build_semiseparable_circuit(u, v);
                U    = circ.matrix;
                norm(u)
                norm(v)
                S       = tril(u*v') + triu(v*u', 1);
                S_tilde = 2*sqrt(N) * U(1:N, 1:N);
                disp(full(S))
                disp(full(S_tilde))
                errs(r)     = norm(S_tilde - S, 2);
                disp(norm(S_tilde - S, 2))
                normS       = norm(S, 2);
                rel_errs(r) = errs(r) / max(normS, eps);

            catch ME
                warning('Run failed (n=%d, type=%d, run=%d): %s', ...
                    n, t, r, ME.message);
                errs(r)     = NaN;
                rel_errs(r) = NaN;
                failed(r)   = true;
            end

            done_tasks = done_tasks + 1;
            print_progress(done_tasks, total_tasks, n, t, r, cfg.num_runs);
        end

        all_errors{t, ni}     = errs;
        all_rel_errors{t, ni} = rel_errs;
        all_failed{t, ni}     = failed;
    end
end

fprintf('\n\nDone.\n');

%% ============================================================
% FOLDERS
%% ============================================================

timestamp     = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
base_folder   = fullfile('results', timestamp);
latest_folder = fullfile('results', 'latest');

if ~exist(base_folder, 'dir');   mkdir(base_folder);   end
if ~exist(latest_folder, 'dir'); mkdir(latest_folder); end

%% ============================================================
% SAVE RAW DATA
%% ============================================================

save(fullfile(base_folder, 'results.mat'), ...
    'all_errors', 'all_rel_errors', 'all_failed', 'cfg', 'timestamp');
save(fullfile(latest_folder, 'results.mat'), ...
    'all_errors', 'all_rel_errors', 'all_failed', 'cfg', 'timestamp');

%% ============================================================
% SUMMARY STATISTICS -> CSV
%% ============================================================

summary_rows = {};
for t = 1:n_types
    for ni = 1:num_n
        n = cfg.n_range(ni);
        e = all_errors{t, ni};
        re = all_rel_errors{t, ni};
        ok = ~all_failed{t, ni};

        summary_rows(end+1, :) = { ...
            cfg.type_names{t}, n, 2^n, ...
            sum(ok), sum(~ok), ...
            mean(e(ok)), std(e(ok)), median(e(ok)), ...
            mean(re(ok)), std(re(ok)), median(re(ok)) }; 
    end
end

summary_table = cell2table(summary_rows, 'VariableNames', ...
    {'Type','n','N','NumOK','NumFailed', ...
     'MeanAbsErr','StdAbsErr','MedianAbsErr', ...
     'MeanRelErr','StdRelErr','MedianRelErr'});

writetable(summary_table, fullfile(base_folder, 'summary.csv'));
writetable(summary_table, fullfile(latest_folder, 'summary.csv'));

disp(summary_table);

%% ============================================================
% PLOTTING
%% ============================================================

% 1. Individual type plots
for t = 1:n_types
    fig = make_scatter_fig(cfg.n_range, all_errors(t, :), ...
        colors(t, :), markers{t}, cfg.type_names{t}, ...
        cfg.type_names{t}, false);
    save_figure_pair(fig, base_folder, latest_folder, ...
        sprintf('scatter_type_%d', t), cfg.save_formats);
    close(fig);
end

% 2. Overlay plot (all types)
%
% Legend handles are built from dedicated proxy markers (plotted at NaN,
% so they're invisible in the axes) rather than inferred from the first
% real data point of each series. Tying legend style to the first data
% point breaks if that point is missing (e.g. a failed run got dropped)
% or if handle/plot order shifts, which is what caused the color/marker
% mismatches in the legend.
fig = figure; hold on;

h = gobjects(n_types, 1);
for t = 1:n_types
    h(t) = scatter(NaN, NaN, 80, ...
        'Marker', markers{t}, ...
        'MarkerEdgeColor', colors(t, :), ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5, ...
        'DisplayName', cfg.type_names{t});
end

for t = 1:n_types
    for ni = 1:num_n
        n = cfg.n_range(ni);
        y = all_errors{t, ni};
        y = y(~all_failed{t, ni});
        x = n * ones(size(y));

        scatter(x, y, 80, ...
            'Marker', markers{t}, ...
            'MarkerEdgeColor', colors(t, :), ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
end

xlabel('n');
ylabel('||S - \tilde{S}||_2');
title('Semiseparable Error (All Generators)');
grid on;
legend(h, 'Location', 'best');

save_figure_pair(fig, base_folder, latest_folder, 'scatter_overlay', cfg.save_formats);
close(fig);

%% ============================================================
% END TIMER
%% ============================================================

elapsed = toc;
fprintf('\nExperiment completed in %.2f seconds (%.2f min)\n', elapsed, elapsed/60);
fprintf('Timestamp: %s\n', timestamp);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function print_progress(done, total, n, t, r, num_runs)
    percent  = 100 * done / total;
    bar_len  = 30;
    filled   = round(bar_len * done / total);
    bar_str  = [repmat('#', 1, filled), repmat('-', 1, bar_len - filled)];

    fprintf('\r[%s] %6.2f%% | n=%d | type=%d | run=%d/%d', ...
        bar_str, percent, n, t, r, num_runs);

    if done == total
        fprintf('\n');
    end
end

function fig = make_scatter_fig(n_range, err_row, color, marker, ttl, ~, ~)
    fig = figure; hold on;
    for ni = 1:numel(n_range)
        n = n_range(ni);
        y = err_row{ni};
        scatter(n * ones(size(y)), y, 80, ...
            'Marker', marker, ...
            'MarkerEdgeColor', color, ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.5);
    end
    title(ttl);
    xlabel('n');
    ylabel('||S - \tilde{S}||_2');
    grid on;
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

function [u, v] = local_generate_uv(n, type)

    N = 2^n;

    switch type

        case 1
            u = randn(N, 1); u = u / norm(u);
            v = randn(N, 1); v = v / norm(v);

        case 2
            T = diag(randn(N, 1)) + ...
                diag(randn(N-1, 1), 1) + ...
                diag(randn(N-1, 1), -1);
            
            S=inv(T);
            [u, v] = uv(S);
            v = v';
            u = u / norm(u);
            v = v / norm(v);
            u=u';
            v=v';

        case 3
            [u,v]=gen_ornstein_uhlenbeck(n, 3.0);

        otherwise
            error('local_generate_uv:badType', 'Unknown type %d', type);
    end
end

