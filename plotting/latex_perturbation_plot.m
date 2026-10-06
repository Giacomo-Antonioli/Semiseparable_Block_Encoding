%% MATLAB Script: Read type*_eps_sweep.mat / type*_N_sweep.mat and Generate LaTeX Plots
%
% Reads directly from the "runs" folder produced by experiments/symmetric_perturbation.m
% (files named type1_eps_sweep.mat, type1_N_sweep.mat, type2_..., etc.)
% and never re-runs any trials. Each file contains, per type t:
%
%   eps_list (or N_list), t, type_err (M x num_trials), type_bound_ok, ...
%
% We take err(i) = max(type_err(i,:)) across trials for each sweep point,
% matching how the benchmark script computed err_eps(ei,t)/err_N(ni,t).
%
clear; clc;
root = setup_paths();

%% --- Values not stored in the per-type files (only in cfg) ---
% Adjust these if your run used different settings than the defaults
% in experiments/symmetric_perturbation.m.
BOUND_CONST = 2*sqrt(2);   % cfg.bound_const
EPS_FIXED   = 1e-6;        % cfg.eps_fixed (used as the reference line in the N-plot)
N_FIXED     = 2^3;         % 2^cfg.n_fixed (used only in the eps-plot title)

type_names = {'Random', 'Tridiagonal', 'Exponential'};
colors = {'blue', 'red', 'green!60!black', 'orange', 'purple'};
marks  = {'*', 'square*', 'triangle*', 'diamond*', 'pentagon*'};

%% 1. Open the GUI Directory Chooser (point this at the "runs" folder)
selected_dir = uigetdir(fullfile(root, 'results'), ...
    'Select the "runs" folder containing type*_eps_sweep.mat / type*_N_sweep.mat');
if isequal(selected_dir, 0)
    disp('User selected Cancel.');
    return;
end
fprintf('Scanning directory: %s\n\n', selected_dir);

%% 2. Read the epsilon-sweep files
eps_files = dir(fullfile(selected_dir, 'type*_eps_sweep.mat'));
if isempty(eps_files)
    error('No files matching "type*_eps_sweep.mat" found in the selected folder.');
end

eps_by_type = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:length(eps_files)
    filepath = fullfile(eps_files(k).folder, eps_files(k).name);
    d = load(filepath, 't', 'eps_list', 'type_err');
    err_vec = max(d.type_err, [], 2)';   % 1 x M, max over trials
    eps_by_type(d.t) = struct('x', d.eps_list, 'err', err_vec);
end

%% 3. Read the N-sweep files
N_files = dir(fullfile(selected_dir, 'type*_N_sweep.mat'));
if isempty(N_files)
    error('No files matching "type*_N_sweep.mat" found in the selected folder.');
end

N_by_type = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:length(N_files)
    filepath = fullfile(N_files(k).folder, N_files(k).name);
    d = load(filepath, 't', 'N_list', 'type_err');
    err_vec = max(d.type_err, [], 2)';   % 1 x M, max over trials
    N_by_type(d.t) = struct('x', d.N_list, 'err', err_vec);
end

%% 4. Build the "Error vs epsilon" figure
tex_eps = build_scatter_figure( ...
    eps_by_type, type_names, colors, marks, ...
    struct( ...
        'xlabel',      '$\varepsilon$', ...
        'ylabel',      '$\|S-\tilde S\|_2$', ...
        'log_basis_x', 10, ...
        'caption',     sprintf('Perturbation error versus $\\varepsilon$ ($N=%d$).', N_FIXED), ...
        'label',       'fig:err_vs_eps', ...
        'ref_line',    struct('type', 'diag', 'const', BOUND_CONST, ...
                               'legend', '$2\sqrt{2}\varepsilon$') ...
    ));

%% 5. Build the "Error vs N" figure
tex_N = build_scatter_figure( ...
    N_by_type, type_names, colors, marks, ...
    struct( ...
        'xlabel',      '$N$', ...
        'ylabel',      '$\|S-\tilde S\|_2$', ...
        'log_basis_x', 2, ...
        'caption',     sprintf('Perturbation error versus $N$ ($\\varepsilon=%.1e$).', EPS_FIXED), ...
        'label',       'fig:err_vs_N', ...
        'ref_line',    struct('type', 'const', 'const', BOUND_CONST*EPS_FIXED, ...
                               'legend', '$2\sqrt{2}\varepsilon$') ...
    ));

%% 6. Output the results
full_latex = [tex_eps, sprintf('\n\n'), tex_N];

fprintf('--- GENERATED LATEX CODE ---\n\n');
disp(full_latex);
fprintf('----------------------------\n');

tex_filename = fullfile(selected_dir, 'perturbation_plots.tex');
fid = fopen(tex_filename, 'w');
if fid ~= -1
    fprintf(fid, '%s', full_latex);
    fclose(fid);
    fprintf('Successfully saved LaTeX block to:\n%s\n', tex_filename);
else
    warning('Could not write out the .tex file automatically.');
end


%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function tex = build_scatter_figure(data_map, type_names, colors, marks, opts)
% data_map : containers.Map keyed by type index t -> struct('x', x_list, 'err', err_vec)

keys_list = sort(cell2mat(data_map.keys));

tex = sprintf([...
    '\\begin{figure}[htbp]\n', ...
    '  \\centering\n', ...
    '  \\begin{tikzpicture}\n', ...
    '  \\begin{axis}[\n', ...
    '      width=0.7\\linewidth,\n', ...
    '      height=7.0cm,\n', ...
    '      xmode=log,\n', ...
    '      ymode=log,\n', ...
    '      log basis x=%d,\n', ...
    '      xlabel={%s},\n', ...
    '      ylabel={%s},\n', ...
    '      grid=both,\n', ...
    '      grid style={dotted},\n', ...
    '      minor tick num=0,\n', ...
    '      enlarge x limits=false,\n', ...
    '      legend style={\n', ...
    '          font=\\footnotesize,\n', ...
    '          draw=none,\n', ...
    '          fill=none,\n', ...
    '          at={(0.02,0.98)},\n', ...
    '          anchor=north west\n', ...
    '      }\n', ...
    '  ]\n\n'], opts.log_basis_x, opts.xlabel, opts.ylabel);

all_x_first = [];
all_x_last  = [];

for t = keys_list
    entry = data_map(t);

    % Sort by x (crucial for clean PGFPlots lines)
    [x_sorted, idx] = sort(entry.x);
    y_sorted = entry.err(idx);

    if t <= numel(type_names)
        t_label = type_names{t};
    else
        t_label = sprintf('Type %d', t);
    end

    c_idx = mod(t-1, numel(colors)) + 1;
    m_idx = mod(t-1, numel(marks)) + 1;

    tex = [tex, sprintf('  %% --- %s ---\n', t_label)];
    tex = [tex, sprintf('  \\addplot[%s, thick, mark=%s] coordinates {\n  ', ...
                          colors{c_idx}, marks{m_idx})];
    for idx2 = 1:numel(x_sorted)
        tex = [tex, sprintf('(%.6e, %e) ', x_sorted(idx2), y_sorted(idx2))];
    end
    tex = [tex, sprintf('\n  };\n  \\addlegendentry{%s}\n\n', t_label)];

    all_x_first = [all_x_first, x_sorted(1)];
    all_x_last  = [all_x_last, x_sorted(end)];
end

% Reference line (theoretical bound), spanning the overall x range
if isfield(opts, 'ref_line') && ~isempty(opts.ref_line)
    x0 = min(all_x_first);
    x1 = max(all_x_last);
    switch opts.ref_line.type
        case 'diag'
            y0 = opts.ref_line.const * x0;
            y1 = opts.ref_line.const * x1;
        case 'const'
            y0 = opts.ref_line.const;
            y1 = opts.ref_line.const;
    end
    tex = [tex, sprintf('  %% --- Theoretical bound ---\n')];
    tex = [tex, sprintf('  \\addplot[black, dashed, thick, mark=none] coordinates {\n  ')];
    tex = [tex, sprintf('(%.6e, %e) (%.6e, %e) ', x0, y0, x1, y1)];
    tex = [tex, sprintf('\n  };\n  \\addlegendentry{%s}\n\n', opts.ref_line.legend)];
end

tex = [tex, sprintf([...
    '  \\end{axis}\n', ...
    '  \\end{tikzpicture}\n', ...
    '  \\caption{%s}\n', ...
    '  \\label{%s}\n', ...
    '\\end{figure}\n'], opts.caption, opts.label)];

end