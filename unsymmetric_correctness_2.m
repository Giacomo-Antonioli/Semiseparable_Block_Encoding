% BENCHMARK SCRIPT: Unsymmetric Random Plot Generator (10-Point Star)
% Now runs 3 trials per N, saves all data per trial, plots the MEAN
% error per N, and prints min/max of each generator per trial.
clear;
clc;
rng(46);

%% ============================================================
% CONFIGURATION
%% ============================================================
cfg.n_range        = 2:6;                               % N = 2^n -> 4, 8, 16, 32, 64, 128
cfg.n_trials       = 1;                                  % Number of trials per N
cfg.enable_saving  = true;                               % Saves .mat workspace data (generators/matrices)
cfg.save_dir       = fullfile('results_unsymmetric_random', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));

num_n         = length(cfg.n_range);
errors_trials = NaN(num_n, cfg.n_trials);   % raw error per N per trial
errors_mean   = NaN(num_n, 1);              % mean error per N (used for plotting)

% TikZ formatting style using the 10-pointed star mark
latex_color   = 'magenta';
latex_marker  = 'star, mark options={star points=10}';
latex_label   = 'Random Unsymmetric (mean of 3 trials)';

if cfg.enable_saving && ~exist(cfg.save_dir, 'dir')
    mkdir(cfg.save_dir);
end

%% ============================================================
% LIVE PLOTS INITIALIZATION (Will NOT be saved to disk)
%% ============================================================
fig_err = figure('Name', 'Unsymmetric Random Correctness', 'NumberTitle', 'off');
ax_err = axes(fig_err); hold(ax_err, 'on');

h_err = scatter(ax_err, NaN, NaN, 80, ...
    'Marker', 'p', 'MarkerEdgeColor', [0.8 0 0.5], ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
    'DisplayName', latex_label);

set(ax_err, 'XScale', 'log', 'XLim', [2, 256]);
xticks(ax_err, 2.^cfg.n_range);
xticklabels(ax_err, string(2.^cfg.n_range));
xlabel(ax_err, 'N'); ylabel(ax_err, '||S - S_{tilde}||_2 (mean of trials)');
title(ax_err, 'Unsymmetric Random Block Encoding Absolute Error');
grid(ax_err, 'on'); legend(ax_err, h_err, 'Location', 'northwest');
drawnow;

%% ============================================================
% CORE BENCHMARK LOOP
%% ============================================================
fprintf('Starting Unsymmetric Random Verification via genera_semisep...\n\n');

for ni = 1:num_n
    n = cfg.n_range(ni);
    N = 2^n;
    fprintf('==================== Processing N = %d ====================\n', N);

    for t = 1:cfg.n_trials
        fprintf('--- N = %d | Trial %d/%d ---\n', N, t, cfg.n_trials);

        % 1. Fetch random generators and matrix directly from your function
        [S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);

        % 1b. Print min/max of each generator for this trial
        fprintf('  Generator ranges (Trial %d, N=%d):\n', t, N);
        fprintf('    u: min = %.6e, max = %.6e, max/min = %.6e\n', min(abs(u(:))), max(abs(u(:))),max(abs(u(:)))/min(abs(u(:))));
        fprintf('    v: min = %.6e, max = %.6e, max/min = %.6e\n', min(abs(v(:))), max(abs(v(:))),max(abs(v(:)))/min(abs(v(:))));
        fprintf('    x: min = %.6e, max = %.6e, max/min = %.6e\n', min(abs(x(:))), max(abs(x(:))),max(abs(x(:)))/min(abs(x(:))));
        fprintf('    y: min = %.6e, max = %.6e, max/min = %.6e\n', min(abs(y(:))), max(abs(y(:))),max(abs(y(:)))/min(abs(y(:))));
        fprintf('    S: min = %.6e, max = %.6e, max/min = %.6e\n', min(abs(S(:))), max(abs(S(:))),max(abs(S(:)))/min(abs(S(:))));
        fprintf('   cond(S) = %.6e', cond(S));
        

        % 2. Build the quantum circuit
        try
            circ = build_unsymmetric_semiseparable_circuit(u, v, x, y,nu,nv,nx,ny);
            dim = 2^circ.nbQubits;
            psi = zeros(dim, 1);
            S_tilde = zeros(N);

            % 3. Extract matrix column by column via simulation
            for j = 1:N
                j
                psi(:) = 0;
                psi(j) = 1;
                sim = circ.simulate(psi);
                out = sim.states;
                S_tilde(:, j) = 2 * sqrt(N) * out(1:N);
            end

            % 4. Compute absolute spectral error
            err = norm(S_tilde - S, 2);
        catch ME
            warning('Execution failed for N=%d, Trial=%d: %s', N, t, ME.message);
            err = NaN;
            S_tilde = [];
        end

        errors_trials(ni, t) = err;
        fprintf('  Trial %d error: %.6e\n', t, err);

        % 5. Save variables to disk for this trial if enabled
        if cfg.enable_saving && ~isnan(err)
            save_file = fullfile(cfg.save_dir, ...
                sprintf('N%d_trial%d_unsymmetric_random.mat', N, t));
            save(save_file, 'u', 'v', 'x', 'y', 'S', 'S_tilde', 'err', 'N', 't', '-v7.3');
        end
    end

    % 6. Compute mean error across trials for this N
    mean_err = mean(errors_trials(ni, :), 'omitnan');
    errors_mean(ni) = mean_err;
    fprintf('N = %d | Mean error over %d trials: %.6e\n\n', N, cfg.n_trials, mean_err);

    % 7. Update Live Plot with the MEAN error only (No image saving executed)
    if ~isnan(mean_err)
        scatter(ax_err, N, mean_err, 80, ...
            'Marker', 'p', 'MarkerEdgeColor', [0.8 0 0.5], ...
            'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        drawnow;
    end
end

%% ============================================================
% SAVE SUMMARY OF ALL TRIALS AND MEANS
%% ============================================================
if cfg.enable_saving
    summary_file = fullfile(cfg.save_dir, 'summary_all_N.mat');
    save(summary_file, 'cfg', 'errors_trials', 'errors_mean', '-v7.3');
end

%% ============================================================
% LATEX FORMATTED OUTPUT GENERATOR (uses MEAN error per N)
%% ============================================================
fprintf('\n============================================================\n');
fprintf('LATEX COORD BLOCK: Copy and paste this into your tikzaxis environment\n');
fprintf('============================================================\n\n');

fprintf('  %% --- %s Plot ---\n', latex_label);
fprintf('  \\addplot[%s, thick, mark=%s, only marks] coordinates {\n  ', ...
        latex_color, latex_marker);

for ni = 1:num_n
    N_val = 2^cfg.n_range(ni);
    fprintf('(%d, %.6e) ', N_val, errors_mean(ni));
end

fprintf('\n  };\n  \\addlegendentry{%s}\n\n', latex_label);
fprintf('============================================================\n');