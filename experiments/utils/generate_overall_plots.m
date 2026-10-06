function generate_overall_plots(base_folder, plots_folder, latest_folder, cfg, markers, colors)
    fprintf('Compiling data points to construct structural master figures...\n');
    n_types = numel(cfg.type_names);
    
    fig_final_err = figure('Name', 'Overall: Semiseparable Error Portfolio', 'NumberTitle', 'off', 'Visible', 'off');
    ax_final_err = axes(fig_final_err); hold(ax_final_err, 'on');
    
    h_leg = gobjects(n_types, 1);
    for t = 1:n_types
        h_leg(t) = scatter(ax_final_err, NaN, NaN, 80, 'Marker', markers{t}, ...
            'MarkerEdgeColor', colors(t, :), 'LineWidth', 1.5, 'DisplayName', cfg.type_names{t});
    end
    
    for e = 1:cfg.num_experiments
        exp_folder = fullfile(base_folder, sprintf('exp_%d', e));
        runs_folder = fullfile(exp_folder, 'runs');
        
        if ~exist(runs_folder, 'dir'); continue; end
        files = dir(fullfile(runs_folder, 'n*.mat'));
        
        for f = 1:numel(files)
            data = load(fullfile(runs_folder, files(f).name));
            if ~data.failed
                scatter(ax_final_err, 2^data.n, data.errs, 80, ...
                    'Marker', markers{data.t_save}, 'MarkerEdgeColor', colors(data.t_save, :), ...
                    'MarkerFaceColor', 'none', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
    end
    
  xlabel(ax_final_err, 'N');
ylabel(ax_final_err, '$\|S - \tilde{S}\|_2$', 'Interpreter', 'latex');
title(ax_final_err, 'Error');

% x-axis: logarithmic scale with powers of two labels
set(ax_final_err, 'XScale', 'log');
xticks(ax_final_err, 2.^cfg.n_range);
xticklabels(ax_final_err, arrayfun(@(x) sprintf('%d',2^x), ...
    cfg.n_range, 'UniformOutput', false));
set(ax_final_err, 'TickLabelInterpreter', 'latex');

grid(ax_final_err, 'on');
legend(ax_final_err, h_leg, 'Location', 'eastoutside');
    
    save_figure_pair(fig_final_err, plots_folder, latest_folder, 'consolidated_error_master', cfg.save_formats);
    close(fig_final_err);
    fprintf('Overall master figure rendering processed successfully.\n');
end
