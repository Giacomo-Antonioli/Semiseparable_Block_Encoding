% BENCHMARK SCRIPT: Single Unsymmetric Random Plot Generator (10-Point Star)
clear;
clc;
rng(46);

%% ============================================================
% CONFIGURATION
%% ============================================================
cfg.n_range        = 2:4;                               % N = 2^n -> 4, 8, 16, 32, 64, 128
cfg.enable_saving  = true;                              % Saves .mat workspace data (generators/matrices)
cfg.save_dir       = fullfile('results_unsymmetric_random', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));

num_n   = length(cfg.n_range);
errors  = zeros(num_n, 1);

% TikZ formatting style using the 10-pointed star mark
latex_color   = 'magenta';
latex_marker  = 'star, mark options={star points=10}';
latex_label   = 'Random Unsymmetric';

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
xlabel(ax_err, 'N'); ylabel(ax_err, '||S - S_{tilde}||_2');
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
    fprintf('Processing N = %d...\n', N);
    
    % 1. Fetch random generators and matrix directly from your function
    [S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);

    u=u/nu;
    v=v/nv;
    x=x/nx;
    y=y/ny;

    
    % 2. Build the quantum circuit
    try
        [circ,t1,t2] = build_unsymmetric_semiseparable_circuit(u, v, x, y,nu,nv,nx,ny);
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
            S_tilde(:, j) = out(1:N);
        end
        
        % 4. Compute absolute spectral error
        err = norm(S/(sqrt(N)*(nu*nv+nx*ny))-S_tilde, 2);
    catch ME
        warning('Execution failed for N=%d: %s', N, ME.message);
        err = NaN;
        S_tilde = [];
    end
    
    errors(ni) = err;
    
    % 5. Save variables to disk if enabled
    if cfg.enable_saving && ~isnan(err)
        save_file = fullfile(cfg.save_dir, sprintf('N%d_unsymmetric_random.mat', N));
        save(save_file, 'u', 'v', 'x', 'y', 'S', 'S_tilde', 'err', 'N', '-v7.3');
    end
    
    % 6. Update Live Plot (No image saving executed)
    if ~isnan(err)
        scatter(ax_err, N, err, 80, ...
            'Marker', 'p', 'MarkerEdgeColor', [0.8 0 0.5], ...
            'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        drawnow;
    end
end

%% ============================================================
% LATEX FORMATTED OUTPUT GENERATOR
%% ============================================================
fprintf('\n============================================================\n');
fprintf('LATEX COORD BLOCK: Copy and paste this into your tikzaxis environment\n');
fprintf('============================================================\n\n');

fprintf('  %% --- %s Plot ---\n', latex_label);
fprintf('  \\addplot[%s, thick, mark=%s, only marks] coordinates {\n  ', ...
        latex_color, latex_marker);

for ni = 1:num_n
    N_val = 2^cfg.n_range(ni);
    fprintf('(%d, %.6e) ', N_val, errors(ni));
end

fprintf('\n  };\n  \\addlegendentry{%s}\n\n', latex_label);
fprintf('============================================================\n');