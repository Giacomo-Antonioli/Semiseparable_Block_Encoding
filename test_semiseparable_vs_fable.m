% test_semiseparable_vs_fable.m
%
% Compares custom semiseparable block-encoding circuit constructions
% (build_semiseparable_circuit, build_semiseparable_circuit2,
% build_semiseparable_circuit3) against the reference FABLE block
% encoding, for:
%
%   1. Correctness  - relative Frobenius error of the reconstructed
%                      matrix against the ground-truth semiseparable S,
%                      for both circuits, and against each other.
%   2. Resources     - gate count, circuit depth, qubit count, and
%                      single-/multi-controlled gate breakdown.
%   3. Subnormalization factor - r(N) = alpha / ||S||_2, comparing the
%                      block-encoding "looseness" of each construction
%                      against FABLE, as a function of matrix size N.
%
% Requires on the MATLAB path:
%   - fable-qclab (the `fable` reference implementation)
%   - local_generate_uv.m
%   - build_semiseparable_circuit.m / _circuit2.m / _circuit3.m
%   - A circuit class exposing getGateCount, getDepth, getQubitCount,
%     getControlledStats (as methods, e.g. on qclab.QCircuit), and
%     assumed to be implemented identically for both `circ` and the
%     circuit object returned by `fable(...)`.
%
% Outputs:
%   - Console table of per-configuration results.
%   - semiseparable_vs_fable_results.mat containing `results` (struct
%     array) and `T` (table) for later analysis.
%   - Two summary figures (gate count vs. N, subnormalization vs. N).

clear; clc;
addpath('fable/fable-qclab/');

%% ---- Configuration --------------------------------------------------
dims = 2:4;     % log2(N): matrix sizes N = 4, 8, ...
types = 1:3;    % semiseparable-generator types (see local_generate_uv)
tol = 1e-10;    % correctness tolerance (relative Frobenius error)

model_names = {'Unitary', ...
               'Mottonnen', ...
               'Exp_tailored'};
type_names = {'Random', ...
               'Tridiagonal', ...
               'Exponential'};

results = struct('dim', {}, 'N', {}, 'type', {}, ...
    'circ_name', {}, 'alpha_my', {}, 'alpha_fable', {}, ...
    'r_my', {}, 'r_fable', {}, ...
    'gates_my', {}, 'gates_fable', {}, ...
    'depth_my', {}, 'depth_fable', {}, 'qubits_my', {}, ...
    'qubits_fable', {}, 'single_ctrl_my', {}, 'multi_ctrl_my', {},'single_ctrl_fable',{},'multi_ctrl_fable',{});

%% ---- Main sweep -------------------------------------------------------
for dim = dims
    N = 2^dim
    for type = types
        [u, v] = local_generate_uv(dim, type);
        
        % Decide which circuit constructions (k) to run based on type
        % 0: build_semiseparable_circuit
        % 1: build_semiseparable_circuit2
        % 2: build_semiseparable_circuit3
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
                    circ = build_semiseparable_circuit(u, v);
                    u_k = u; v_k = v;
                case 1
                    circ = build_semiseparable_circuit2(u, v);
                    u_k = u; v_k = v;
                case 2
                    [circ, v_k, u_k] = build_semiseparable_circuit3(u, v);
                otherwise
                    error('test:invalidModel', ...
                        ['Invalid value of k = %d. ' ...
                         'Expected k in {0,1,2}.'], k);
            end
            
            % Ground-truth semiseparable matrix for THIS variant's
            % (u, v), since build_semiseparable_circuit3 may return
            % updated factors.
            S = tril(u_k * v_k') + triu(v_k * u_k', 1);
            
            % Reference FABLE block encoding of the same matrix.
            fableCirc = fable(S, 'cutoff', 0, false);
            
            % --- Reconstruct S from both circuits, column by column ---
            alpha_my = 2 * sqrt(N);   % subnormalization used by circ
            alpha_fable = N;          % subnormalization used by FABLE
            my_S = zeros(N);
            Fable_S = zeros(N);
            psiMy = zeros(2^circ.nbQubits, 1);
            psiFable = zeros(2^fableCirc.nbQubits, 1);
            
            for col = 1:N
                psiMy(:) = 0; psiMy(col) = 1;
                psiFable(:) = 0; psiFable(col) = 1;
                
                simMy = circ.simulate(psiMy);
                outMy = simMy.states;
                my_S(:, col) = alpha_my * outMy(1:N);
                
                simFable = fableCirc.simulate(psiFable);
                outFable = simFable.states;
                Fable_S(:, col) = alpha_fable * outFable(1:N);
            end
            
            % --- Correctness metrics ------------------------------------
            normS_fro = norm(S, 'fro');
            rel_err_my = norm(my_S - S, 'fro') / normS_fro
            rel_err_fable = norm(Fable_S - S, 'fro') / normS_fro
            rel_err_my_v_fable = norm(my_S - Fable_S, 'fro') / normS_fro
            
            if rel_err_my > tol
                warning('test:mismatchMy', ...
                    ['dim=%d type=%d k=%d: circuit reconstruction of ' ...
                     'S off by relative error %.3e (tol %.1e)'], ...
                    dim, type, k, rel_err_my, tol);
            end
            if rel_err_fable > tol
                warning('test:mismatchFable', ...
                    ['dim=%d type=%d k=%d: FABLE reconstruction of ' ...
                     'S off by relative error %.3e (tol %.1e)'], ...
                    dim, type, k, rel_err_fable, tol);
            end
            
            % --- Subnormalization-factor comparison ----------------------
            % r(N) = alpha / ||S||_2 : how loose the block encoding is
            % relative to the matrix's own operator norm (smaller is
            % better / tighter).
            normS_op = norm(S);   % operator (spectral) norm
            r_my = alpha_my / normS_op;
            r_fable = alpha_fable / normS_op;
            
            % --- Resource metrics ------------------------------------------
            gates_my = circ.getGateCount();
            gates_fable = fableCirc.getGateCount();
            depth_my = circ.getDepth();
            depth_fable = fableCirc.getDepth();
            qubits_my = circ.getQubitCount();
            qubits_fable = fableCirc.getQubitCount();
            ctrl_my = circ.getControlledStats();
            ctrl_fable = fableCirc.getControlledStats(); % <- Add this line
            
            % --- Record ------------------------------------------------------
            results(end+1) = struct( ...                        
                'dim', dim, 'N', N, 'type', type_names{type},  ...
                'circ_name', model_names{k+1}, ...
                'alpha_my', alpha_my, 'alpha_fable', alpha_fable, ...
                'r_my', r_my, 'r_fable', r_fable, ...
                'gates_my', gates_my, 'gates_fable', gates_fable, ...
                'depth_my', depth_my, 'depth_fable', depth_fable, ...
                'qubits_my', qubits_my, 'qubits_fable', qubits_fable, ...
                'single_ctrl_my', ctrl_my.single, ...
                'multi_ctrl_my', ctrl_my.multi,...
                'single_ctrl_fable', ctrl_fable.single, ... % <- Add this line
                'multi_ctrl_fable', ctrl_fable.multi);       
        end
    end
end

%% ---- Report -------------------------------------------------------------
T = struct2table(results);
disp(T);
fprintf('\nAll %d configurations checked (tol = %.1e).\n', ...
    height(T), tol);
save('semiseparable_vs_fable_results.mat', 'results', 'T');

%% ---- Gate-count comparison plot ------------------------------------------
%% ---- Resource Plots (Including FABLE) -----------------------------------
type_names = {'Random', 'Tridiagonal', 'Exponential'};
types = 1:3;

% Configuration matrix for the four resource plots:
% {Figure Name, Y-Axis Label, Custom Data Column, FABLE Data Column, Plot FABLE line?}
plot_configs = { ...
    {'Gate count vs N', 'gate count', 'gates_my', 'gates_fable', true}, ...
    {'Circuit Depth vs N', 'circuit depth', 'depth_my', 'depth_fable', true}, ...
    {'Single-Controlled Gates vs N', 'single-controlled gate count', 'single_ctrl_my', 'single_ctrl_fable', true}, ...
    {'Multi-Controlled Gates vs N', 'multi-controlled gate count', 'multi_ctrl_my', 'multi_ctrl_fable', true} ...
};

for m = 1:numel(plot_configs)
    config = plot_configs{m};
    figure('Name', config{1});
    
    for type = types
        % Custom Layout: Type 1 top-left, Type 2 bottom-left, Type 3 full right column
        if type == 1
            subplot(2, 2, 1);
        elseif type == 2
            subplot(2, 2, 3);
        elseif type == 3
            subplot(2, 2, [2, 4]);
        end
        
        idx = strcmp(T.type, type_names{type});
        
        % Scatter plot for custom constructions grouped by model name
        gscatter(T.N(idx), T.(config{3})(idx), T.circ_name(idx));
        hold on;
        
        % Overlay reference FABLE line
        if config{5}
            Ns = unique(T.N(idx));
            fable_mean = arrayfun(@(n) mean(T.(config{4})(idx & T.N == n)), Ns);
            plot(Ns, fable_mean, 'k--o', 'DisplayName', 'FABLE');
        end
        
        xlabel('N'); ylabel(config{2});
        title(sprintf('Type: %s', type_names{type}));
        legend('Location', 'best');
    end
end

%% ---- Subnormalization-factor comparison plot -----------------------------
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