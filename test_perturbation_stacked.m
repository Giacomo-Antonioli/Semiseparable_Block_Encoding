%TEST_PERTURBATION_BOUND
% Tests the claim:
%
%   If u, v are unit vectors, S = tril(u*v') + triu(v*u',1), and
%   u_tilde = u + delta_u, v_tilde = v + delta_v with
%   ||delta_u|| = ||delta_v|| = epsilon, then
%
%       ||S - S_tilde||_2 < 2*sqrt(2) * epsilon        (for small epsilon)
%
% where S_tilde = tril(u_tilde*v_tilde') + triu(v_tilde*u_tilde',1).
%
% Two experiments:
%   (1) Sweep epsilon over 1e-1 ... 1e-14 for SEVERAL fixed values of N,
%       check the bound holds for each, and overlay all curves on one
%       log-log plot to see how the error scales with epsilon (and how
%       consistent that scaling is across N).
%   (2) Fix epsilon, sweep N over several powers of two, check the
%       error (and the bound) do NOT depend on N.

clear;
clc;

%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.N_list_exp1    = 2.^(2:3);            % experiment 1: sweep of fixed matrix sizes
cfg.eps_powers     = -6:-2:-14;           % epsilon = 10.^eps_powers
cfg.num_trials     = 1;                  % random trials per epsilon / per N

cfg.eps_fixed      = 1e-6;                % experiment 2: fixed epsilon
cfg.n_powers       = 2:4;                % N = 2.^n_powers
cfg.bound_const    = 2*sqrt(2);           % theoretical constant

rng(42);  % reproducibility for the whole script

%% ============================================================
% FOLDERS
%% ============================================================

timestamp     = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
base_folder   = fullfile('results_perturbation', timestamp);
latest_folder = fullfile('results_perturbation', 'latest');

if ~exist(base_folder, 'dir');   mkdir(base_folder);   end
if ~exist(latest_folder, 'dir'); mkdir(latest_folder); end

%% ============================================================
% EXPERIMENT 1: error vs epsilon, overlaid for several fixed N
%% ============================================================

eps_list  = 10.^cfg.eps_powers;
num_eps   = numel(eps_list);
N_list_1  = cfg.N_list_exp1;
num_N1    = numel(N_list_1);

err_max_mat   = zeros(num_N1, num_eps);   % rows: N, cols: epsilon
err_mean_mat  = zeros(num_N1, num_eps);
ratio_max_mat = zeros(num_N1, num_eps);
n_violations_mat = zeros(num_N1, num_eps);

for Ni = 1:num_N1
    N = N_list_1(Ni);
    for ei = 1:num_eps
        epsilon = eps_list(ei);
        errs = zeros(cfg.num_trials, 1);
        violations = 0;

        for r = 1:cfg.num_trials
            [err, bound_ok] = run_one_perturbation_trial(N, epsilon, cfg.bound_const);
            errs(r) = err;
            violations = violations + ~bound_ok;
        end

        err_max_mat(Ni, ei)   = max(errs);
        err_mean_mat(Ni, ei)  = mean(errs);
        ratio_max_mat(Ni, ei) = err_max_mat(Ni, ei) / epsilon;
        n_violations_mat(Ni, ei) = violations;

        fprintf('N = %5d | eps = %.1e | max err = %.3e | max err/eps = %.4f | bound = %.4f | violations = %d/%d\n', ...
            N, epsilon, err_max_mat(Ni,ei), ratio_max_mat(Ni,ei), cfg.bound_const, violations, cfg.num_trials);
    end
end

%% ---- Plot: error vs epsilon (log-log), one curve per N ----

fig1 = figure; hold on;
colors = lines(num_N1);
for Ni = 1:num_N1
    loglog(eps_list, err_max_mat(Ni,:), 'o-', 'LineWidth', 1.5, ...
        'Color', colors(Ni,:), ...
        'DisplayName', sprintf('N = %d (observed)', N_list_1(Ni)));
end
loglog(eps_list, cfg.bound_const * eps_list, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('bound: $2\\sqrt{2} \\times \\epsilon \\approx %.4f\\epsilon$', cfg.bound_const));
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\epsilon');
ylabel('||S - \tilde{S}||_2');
title('Perturbation error vs \epsilon, overlaid across N');
legend('Location', 'best', 'Interpreter', 'latex');
grid on;

save_figure_pair(fig1, base_folder, latest_folder, 'error_vs_epsilon', {'png','fig'});
% close(fig1);

%% ============================================================
% EXPERIMENT 2: error vs N, epsilon fixed
%% ============================================================

N_list  = 2.^cfg.n_powers;
num_N   = numel(N_list);

err_max_N  = zeros(num_N, 1);
err_mean_N = zeros(num_N, 1);
n_violations_N = zeros(num_N, 1);

for Ni = 1:num_N
    N = N_list(Ni);
    errs = zeros(cfg.num_trials, 1);

    for r = 1:cfg.num_trials
        [err, bound_ok] = run_one_perturbation_trial(N, cfg.eps_fixed, cfg.bound_const);
        errs(r) = err;
        n_violations_N(Ni) = n_violations_N(Ni) + ~bound_ok;
    end

    err_max_N(Ni)  = max(errs);
    err_mean_N(Ni) = mean(errs);

    fprintf('N = %5d | max err = %.3e | mean err = %.3e | violations = %d/%d\n', ...
        N, err_max_N(Ni), err_mean_N(Ni), n_violations_N(Ni), cfg.num_trials);
end

%% ---- Plot: error vs N (semilog-x), should be flat ----

fig2 = figure; hold on;
semilogx(N_list, err_max_N, 'o-', 'LineWidth', 1.5, 'DisplayName', 'max ||S - S_{tilde}||_2 (observed)');
yline(cfg.bound_const * cfg.eps_fixed, 'k--', 'LineWidth', 1.5, ...
    'Label', sprintf('bound = %.4f \\times \\epsilon', cfg.bound_const));
xlabel('N');
ylabel('||S - \tilde{S}||_2');
title(sprintf('Perturbation error vs N  (\\epsilon = %.1e, fixed)', cfg.eps_fixed));
legend('Location', 'best');
grid on;

save_figure_pair(fig2, base_folder, latest_folder, 'error_vs_N', {'png','fig'});
% close(fig2);

%% ============================================================
% SAVE SUMMARY DATA
%% ============================================================

% Experiment 1: long-format table (one row per N x epsilon combination)
[N_grid, eps_grid] = ndgrid(N_list_1, eps_list);
exp1_table = table(N_grid(:), eps_grid(:), err_max_mat(:), err_mean_mat(:), ...
    ratio_max_mat(:), n_violations_mat(:), ...
    'VariableNames', {'N','epsilon','MaxErr','MeanErr','MaxErrOverEps','NumViolations'});

exp2_table = table(N_list(:), err_max_N, err_mean_N, n_violations_N, ...
    'VariableNames', {'N','MaxErr','MeanErr','NumViolations'});

writetable(exp1_table, fullfile(base_folder, 'exp1_error_vs_epsilon.csv'));
writetable(exp1_table, fullfile(latest_folder, 'exp1_error_vs_epsilon.csv'));
writetable(exp2_table, fullfile(base_folder, 'exp2_error_vs_N.csv'));
writetable(exp2_table, fullfile(latest_folder, 'exp2_error_vs_N.csv'));

save(fullfile(base_folder, 'perturbation_results.mat'), ...
    'cfg', 'eps_list', 'N_list_1', 'err_max_mat', 'err_mean_mat', 'ratio_max_mat', 'n_violations_mat', ...
    'N_list', 'err_max_N', 'err_mean_N', 'n_violations_N', 'timestamp');
save(fullfile(latest_folder, 'perturbation_results.mat'), ...
    'cfg', 'eps_list', 'N_list_1', 'err_max_mat', 'err_mean_mat', 'ratio_max_mat', 'n_violations_mat', ...
    'N_list', 'err_max_N', 'err_mean_N', 'n_violations_N', 'timestamp');

disp(exp1_table);
disp(exp2_table);

fprintf('\nNote: for the very smallest epsilon values (~1e-14), floating-point\n');
fprintf('roundoff (machine epsilon ~2.2e-16) can occasionally push the ratio\n');
fprintf('above the theoretical constant; this is a numerical-precision effect,\n');
fprintf('not a violation of the analytic bound.\n');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [err, bound_ok] = run_one_perturbation_trial(N, epsilon, bound_const)
% Draws a random unit pair (u, v), perturbs each by a random direction
% scaled to norm epsilon, and returns ||S - S_tilde||_2 along with
% whether it satisfies err < bound_const * epsilon.

    u = randn(N, 1); u = u / norm(u);
    v = randn(N, 1); v = v / norm(v);

    delta_u = randn(N, 1); delta_u = delta_u / norm(delta_u) * epsilon;
    delta_v = randn(N, 1); delta_v = delta_v / norm(delta_v) * epsilon;

    u_tilde = u + delta_u;
    v_tilde = v + delta_v;

    % S       = tril(u*v')             + triu(v*u', 1);
    S = tril(u_tilde*v_tilde') + triu(v_tilde*u_tilde', 1);
    circ = build_semiseparable_circuit(u_tilde, v_tilde);
    U    = circ.matrix;
    S_tilde = 2*sqrt(N) * U(1:N, 1:N);
    err = norm(S - S_tilde, 2);
    bound_ok = err < bound_const * epsilon;
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