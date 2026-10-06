% SUBNORMALIZATION_VS_FABLE
%
% Calculates the subnormalization factors (r_my and r_fable) for custom 
% semiseparable block-encoding circuit constructions compared to FABLE.
% Bypasses all circuit generation and simulation.
%
% Subnormalization factor - r(N) = alpha / ||S||_2, comparing the
% block-encoding "looseness" of each construction against FABLE, 
% as a function of matrix size N.
%
% Requires on the MATLAB path:
%   - local_generate_uv.m
%   - build_semiseparable_circuit_exp.m (only to fetch updated u, v factors)
%
% OUTPUT: in addition to the interactive MATLAB figure, this script now
% also writes a grouped PGFPlots/LaTeX figure directly from the results
% table T (no external files needed -- everything is already in memory),
% laid out the same way as the MATLAB subplot(2,2,...) grid: Random
% (top-left), Tridiagonal (bottom-left), Exponential (spanning the full
% right column).

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
root = setup_paths();

%% ---- Configuration --------------------------------------------------
dims = 2:10;     % log2(N): matrix sizes N = 4, 8, ...
types = 1:3;    % semiseparable-generator types (see local_generate_uv)

model_names = {'Unitary', ...
               'Mottonnen', ...
               'Exp_tailored'};
type_names = {'Random', ...
               'Tridiagonal', ...
               'Exponential'};

% Simplified results struct just for subnormalization
results = struct('dim', {}, 'N', {}, 'type', {}, ...
    'circ_name', {}, 'alpha_my', {}, 'alpha_fable', {}, ...
    'r_my', {}, 'r_fable', {});

%% ---- Main sweep -------------------------------------------------------
for dim = dims
    N = 2^dim;
    for type = types
        [u, v] = local_generate_uv(dim, type);
        
        % Decide which circuit constructions (k) to run based on type
        % 0: Unitary
        % 1: Mottonnen
        % 2: Exp_tailored
        if type == 1 || type == 2
            models_to_run = [0, 1];
        elseif type == 3
            models_to_run = 2;
        else
            error('Unknown type');
        end
        
        for k = models_to_run
            switch k
                case 0
                    u_k = u; v_k = v;
                case 1
                    u_k = u; v_k = v;
                case 2
                    % Discard the circuit (~), just grab the updated factors
                    [~, v_k, u_k] = build_semiseparable_circuit_exp(u, v);
                otherwise
                    error('test:invalidModel', ...
                        ['Invalid value of k = %d. ' ...
                         'Expected k in {0,1,2}.'], k);
            end
            
            % Ground-truth semiseparable matrix for THIS variant's (u, v)
            S = tril(u_k * v_k') + triu(v_k * u_k', 1);
            
            % --- Subnormalization-factor calculation ---------------------
            alpha_my = 2 * sqrt(N);   % subnormalization used by custom circ
            alpha_fable = N;          % subnormalization used by FABLE
            
            normS_op = norm(full(S));       % operator (spectral) norm
            r_my = alpha_my / normS_op;
            r_fable = alpha_fable / normS_op;
            
            % --- Record --------------------------------------------------
            results(end+1) = struct( ...                        
                'dim', dim, 'N', N, 'type', type_names{type},  ...
                'circ_name', model_names{k+1}, ...
                'alpha_my', alpha_my, 'alpha_fable', alpha_fable, ...
                'r_my', r_my, 'r_fable', r_fable);
        end
    end
end

%% ---- Report -------------------------------------------------------------
T = struct2table(results);
disp(T);
fprintf('\nExtracted subnormalization factors for %d configurations.\n', height(T));
save(fullfile(root, 'results', 'subnormalization_vs_fable', 'semiseparable_vs_fable_subnorm_results.mat'), 'results', 'T');

%% ---- Subnormalization-factor comparison plot (interactive, MATLAB) -----
figure('Name', 'Subnormalization ratio r(N) vs N');
for type = types
    
    % Subplot setup: Type 1 top-left, Type 2 bottom-left, Type 3 spans right column
    if type == 1
        subplot(2, 2, 1);
    elseif type == 2
        subplot(2, 2, 3);
    elseif type == 3
        subplot(2, 2, [2, 4]);
    end
    
    idx = strcmp(T.type, type_names{type});
    
    % Scatter plot grouped by the specific model (circ_name)
    gscatter(T.N(idx), T.r_my(idx), T.circ_name(idx));
    hold on;
    
    % Overlay the FABLE reference line
    Ns = unique(T.N(idx));
    rf = arrayfun(@(n) mean(T.r_fable(idx & T.N == n)), Ns);
    plot(Ns, rf, 'k--o', 'DisplayName', 'FABLE');
    
    xlabel('N'); ylabel('r(N) = \alpha / ||S||_2');
    title(sprintf('Type: %s', type_names{type}));
    legend('Location', 'best');
end


%% ============================================================
% LaTeX (PGFPlots) grouped figure -- same 2x2-with-merged-column
% layout as the MATLAB subplot above, built directly from T.
%% ============================================================

style_map = containers.Map();
style_map('Unitary')      = struct('color', 'blue',           'mark', '*');
style_map('Mottonnen')    = struct('color', 'red',            'mark', 'square*');
style_map('Exp_tailored') = struct('color', 'green!60!black', 'mark', 'triangle*');
style_map('FABLE')        = struct('color', 'black',          'mark', 'o');

N_all_vals = unique(T.N);
xtick_str  = strjoin(arrayfun(@(n) sprintf('%d', n), N_all_vals, 'UniformOutput', false), ',');

panel_tex = cell(1, numel(types));

for type = types
    idx = strcmp(T.type, type_names{type});
    Tt  = T(idx, :);

    % Small panels (Random, Tridiagonal) get half height; the panel that
    % spans both rows (Exponential) gets roughly double + the gap between them
    if type == 3
        height_str = '9.3cm';
    else
        height_str = '4.5cm';
    end

    tex = sprintf([...
        '  \\begin{tikzpicture}\n', ...
        '  \\begin{axis}[\n', ...
        '      width=\\linewidth,\n', ...
        '      height=%s,\n', ...
        '      xmode=log,\n', ...
        '      log basis x=2,\n', ...
        '      xmin=%d, xmax=%d,\n', ...
        '      xtick={%s},\n', ...
        '      xticklabels={%s},\n', ...
        '      xlabel={$N$},\n', ...
        '      ylabel={$r(N) = \\alpha / \\|S\\|_2$},\n', ...
        '      xlabel near ticks,\n', ...
        '      ylabel near ticks,\n', ...
        '      xlabel style={font=\\tiny},\n', ...
        '      ylabel style={font=\\tiny},\n', ...
        '      title style={font=\\tiny, yshift=-2pt},\n', ...
        '      tick label style={font=\\tiny},\n', ...
        '      grid=both,\n', ...
        '      grid style={dotted},\n', ...
        '      minor tick num=0,\n', ...
        '      enlarge x limits=false,\n', ...
        '      title={Type: %s},\n', ...
        '      legend style={\n', ...
        '          font=\\tiny,\n', ...
        '          draw=none,\n', ...
        '          fill=none,\n', ...
        '          at={(0.02,0.02)},\n', ...
        '          anchor=south west,\n', ...
        '      },\n', ...
        '  ]\n\n'], ...
        height_str, min(N_all_vals), max(N_all_vals), ...
        xtick_str, xtick_str, type_names{type});

    circ_list = unique(Tt.circ_name, 'stable');
    for c = 1:numel(circ_list)
        cname = circ_list{c};
        sub = Tt(strcmp(Tt.circ_name, cname), :);
        [Ns_sorted, sidx] = sort(sub.N);
        rs_sorted = sub.r_my(sidx);

        st = style_map(cname);
        legend_label = strrep(cname, '_', '\_');

        % Scatter only (no connecting line) for the circuit variants
        tex = [tex, sprintf('  %% --- %s ---\n', cname)];
        tex = [tex, sprintf('  \\addplot[%s, only marks, mark=%s, mark size=2pt] coordinates {\n  ', ...
                              st.color, st.mark)];
        for kk = 1:numel(Ns_sorted)
            tex = [tex, sprintf('(%d, %e) ', Ns_sorted(kk), rs_sorted(kk))];
        end
        tex = [tex, sprintf('\n  };\n  \\addlegendentry{%s}\n\n', legend_label)];
    end

    % FABLE reference: mean r_fable at each N present for this type
    Ns_fable = unique(Tt.N);
    rf_vals  = arrayfun(@(n) mean(Tt.r_fable(Tt.N == n)), Ns_fable);
    [Ns_fable, sidx] = sort(Ns_fable);
    rf_vals = rf_vals(sidx);

    st = style_map('FABLE');
    tex = [tex, sprintf('  %% --- FABLE reference ---\n')];
    tex = [tex, sprintf('  \\addplot[%s, thick, mark=%s] coordinates {\n  ', ...
                          st.color, st.mark)];
    for kk = 1:numel(Ns_fable)
        tex = [tex, sprintf('(%d, %e) ', Ns_fable(kk), rf_vals(kk))];
    end
    tex = [tex, sprintf('\n  };\n  \\addlegendentry{FABLE}\n\n')];

    tex = [tex, sprintf('  \\end{axis}\n  \\end{tikzpicture}\n')];

    panel_tex{type} = tex;
end

full_latex = sprintf([...
    '%% Requires \\usepackage{pgfplots} (\\pgfplotsset{compat=1.18}) in the preamble\n', ...
    '\\begin{figure}[htbp]\n', ...
    '\\centering\n', ...
    '\\begin{minipage}[t]{0.48\\textwidth}\n', ...
    '  \\centering\n', ...
    '%s', ...
    '  \\vspace{0.3cm}\n\n', ...
    '%s', ...
    '\\end{minipage}\n', ...
    '\\hfill\n', ...
    '\\begin{minipage}[t]{0.48\\textwidth}\n', ...
    '  \\centering\n', ...
    '%s', ...
    '\\end{minipage}\n', ...
    '\\caption{Subnormalization ratio $r(N)=\\alpha/\\|S\\|_2$ versus $N$, compared against FABLE.}\n', ...
    '\\label{fig:subnorm}\n', ...
    '\\end{figure}\n'], ...
    panel_tex{1}, panel_tex{2}, panel_tex{3});

fprintf('\n--- GENERATED LATEX CODE ---\n\n');
disp(full_latex);
fprintf('----------------------------\n');

tex_filename = fullfile(root, 'results', 'subnormalization_vs_fable', 'semiseparable_vs_fable_subnorm_plot.tex');
fid = fopen(tex_filename, 'w');
if fid ~= -1
    fprintf(fid, '%s', full_latex);
    fclose(fid);
    fprintf('Successfully saved LaTeX block to:\n%s\n', tex_filename);
else
    warning('Could not write out the .tex file automatically.');
end