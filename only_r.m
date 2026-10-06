% test_semiseparable_vs_fable.m
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
%   - build_semiseparable_circuit3.m (only to fetch updated u, v factors)

clear; clc;

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
                    [~, v_k, u_k] = build_semiseparable_circuit3(u, v);
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
save('semiseparable_vs_fable_subnorm_results.mat', 'results', 'T');

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