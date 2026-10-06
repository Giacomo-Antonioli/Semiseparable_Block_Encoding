%% MATLAB Script: Load Directory of .mat runs and Generate LaTeX Scatter Plot
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
root = setup_paths();

% 1. Open the GUI Directory Chooser
selected_dir = uigetdir(fullfile(root, 'results'), 'Select the folder containing your .mat files (e.g., the "runs" folder)');
if isequal(selected_dir, 0)
    disp('User selected Cancel.');
    return;
end
fprintf('Scanning directory: %s\n\n', selected_dir);

% 2. Find all relevant .mat files in the selected directory
mat_files = dir(fullfile(selected_dir, 'n*_type*_run*.mat'));
if isempty(mat_files)
    error('No files matching "n*_type*_run*.mat" found in the selected folder.');
end

% Configuration for the plot (Matches your benchmark cfg.type_names)
type_names = {'Random', 'Tridiagonal', 'Exponential'}; 
colors = {'blue', 'red', 'green!60!black', 'orange', 'purple'};
marks  = {'*', 'square*', 'triangle*', 'diamond*', 'pentagon*'};

% Prepare a structure to hold data grouped by type
max_type = 0;
plot_data = struct('N', {}, 'errs', {}, 'rel_errs', {});

% 3. Extract data and group by type (t_save)
for k = 1:length(mat_files)
    filepath = fullfile(mat_files(k).folder, mat_files(k).name);
    
    % Load only the necessary variables to save memory/time
    data = load(filepath, 'n', 't_save', 'errs', 'rel_errs');
    
    t = data.t_save;
    N = 2^(data.n); % Convert n to matrix size N
    
    if t > max_type
        max_type = t;
    end
    
    % Initialize arrays for this type if they don't exist yet
    if length(plot_data) < t || isempty(plot_data(t).N)
        plot_data(t).N = [];
        plot_data(t).errs = [];
        plot_data(t).rel_errs = [];
    end
    
    % Append the data
    plot_data(t).N(end+1) = N;
    plot_data(t).errs(end+1) = data.errs;
    if isfield(data, 'rel_errs')
        plot_data(t).rel_errs(end+1) = data.rel_errs;
    end
end

% 4. Sort the data for each type by N (Crucial for clean PGFPlots lines)
for t = 1:max_type
    if length(plot_data) >= t && ~isempty(plot_data(t).N)
        [sorted_N, idx] = sort(plot_data(t).N);
        plot_data(t).N = sorted_N;
        plot_data(t).errs = plot_data(t).errs(idx);
        if ~isempty(plot_data(t).rel_errs)
            plot_data(t).rel_errs = plot_data(t).rel_errs(idx);
        end
    end
end

% 5. Generate the LaTeX PGFPlots code (Now wrapped in a figure block)
latex_code = sprintf([...
    '\\begin{figure}[htbp]\n', ...
    '  \\centering\n', ...
    '  \\begin{tikzpicture}\n', ...
    '  \\begin{axis}[\n', ...
    '      width=0.7\\linewidth,\n', ...
    '      height=7.0cm,\n', ...
    '      xmode=log,\n', ...
    '      ymode=log,\n', ...
    '      log basis x=2,\n', ...
    '      xlabel={$N$},\n', ...
    '      ylabel={Absolute Error},\n', ...
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
    '  ]\n\n']);

for t = 1:max_type
    if length(plot_data) >= t && ~isempty(plot_data(t).N)
        
        % Assign styling 
        c_idx = mod(t-1, length(colors)) + 1;
        m_idx = mod(t-1, length(marks)) + 1;
        
        % Get the label name
        if t <= length(type_names)
            t_label = type_names{t};
        else
            t_label = sprintf('Type %d', t);
        end
        
        % Write the \addplot header
        latex_code = [latex_code, sprintf('  %% --- %s Plot ---\n', t_label)];
        latex_code = [latex_code, sprintf('  \\addplot[%s, thick, mark=%s] coordinates {\n  ', colors{c_idx}, marks{m_idx})];
        
        % Write the coordinates
        for idx = 1:length(plot_data(t).N)
            latex_code = [latex_code, sprintf('(%d, %e) ', plot_data(t).N(idx), plot_data(t).errs(idx))];
        end
        
        % Close the \addplot and add the legend entry
        latex_code = [latex_code, sprintf('\n  };\n  \\addlegendentry{%s}\n\n', t_label)];
    end
end

% Close the axis, tikzpicture, and figure wrapper safely
latex_code = [latex_code, sprintf([...
    '  \\end{axis}\n', ...
    '  \\end{tikzpicture}\n', ...
    '  \\caption{Semiseparable reconstruction error scaling over matrix dimension $N$.}\n', ...
    '  \\label{fig:recon_error}\n', ...
    '\\end{figure}\n'])];

% 6. Output the results
fprintf('--- GENERATED LATEX CODE ---\n\n');
disp(latex_code);
fprintf('----------------------------\n');

% Save to a .tex file in the selected directory
tex_filename = fullfile(selected_dir, 'overall_scatter_plot.tex');
fid = fopen(tex_filename, 'w');
if fid ~= -1
    fprintf(fid, '%s', latex_code);
    fclose(fid);
    fprintf('Successfully saved LaTeX block to:\n%s\n', tex_filename);
else
    warning('Could not write out the .tex file automatically.');
end