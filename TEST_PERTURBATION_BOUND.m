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
%   (1) Fix N = 32, sweep epsilon over 1e-1 ... 1e-14, check the bound
%       holds and see how the error scales with epsilon.
%   (2) Fix epsilon, sweep N over several powers of two, check the
%       error (and the bound) do NOT depend on N.

clear;
clc;

%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.N_fixed        = 4;                  % experiment 1: fixed matrix size
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
% EXPERIMENT 1: error vs epsilon, N fixed
%% ============================================================

eps_list  = 10.^cfg.eps_powers;
num_eps   = numel(eps_list);

err_max  = zeros(num_eps, 1);
err_mean = zeros(num_eps, 1);
ratio_max = zeros(num_eps, 1);   % error / epsilon, should stay < bound_const
n_violations_eps = zeros(num_eps, 1);

for ei = 1:num_eps
    epsilon = eps_list(ei);
    errs = zeros(cfg.num_trials, 1);

    for r = 1:cfg.num_trials
        [err, bound_ok] = run_one_perturbation_trial(cfg.N_fixed, epsilon, cfg.bound_const);
        errs(r) = err;
        n_violations_eps(ei) = n_violations_eps(ei) + ~bound_ok;
    end

    err_max(ei)   = max(errs);
    err_mean(ei)  = mean(errs);
    ratio_max(ei) = err_max(ei) / epsilon;

    fprintf('eps = %.1e | max err = %.3e | max err/eps = %.4f | bound = %.4f | violations = %d/%d\n', ...
        epsilon, err_max(ei), ratio_max(ei), cfg.bound_const, n_violations_eps(ei), cfg.num_trials);
end

%% ---- Plot: error vs epsilon (log-log) ----

fig1 = figure; hold on;
loglog(eps_list, err_max, 'o-', 'LineWidth', 1.5, 'DisplayName', 'max ||S - S_{tilde}||_2 (observed)');
loglog(eps_list, cfg.bound_const * eps_list, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('bound: %.4f \\times \\epsilon', cfg.bound_const));
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\epsilon');
ylabel('||S - \tilde{S}||_2');
title(sprintf('Perturbation error vs \\epsilon  (N = %d)', cfg.N_fixed));
legend('Location', 'best');
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

exp1_table = table(eps_list(:), err_max, err_mean, ratio_max, n_violations_eps, ...
    'VariableNames', {'epsilon','MaxErr','MeanErr','MaxErrOverEps','NumViolations'});
exp2_table = table(N_list(:), err_max_N, err_mean_N, n_violations_N, ...
    'VariableNames', {'N','MaxErr','MeanErr','NumViolations'});

writetable(exp1_table, fullfile(base_folder, 'exp1_error_vs_epsilon.csv'));
writetable(exp1_table, fullfile(latest_folder, 'exp1_error_vs_epsilon.csv'));
writetable(exp2_table, fullfile(base_folder, 'exp2_error_vs_N.csv'));
writetable(exp2_table, fullfile(latest_folder, 'exp2_error_vs_N.csv'));

save(fullfile(base_folder, 'perturbation_results.mat'), ...
    'cfg', 'eps_list', 'err_max', 'err_mean', 'ratio_max', 'n_violations_eps', ...
    'N_list', 'err_max_N', 'err_mean_N', 'n_violations_N', 'timestamp');
save(fullfile(latest_folder, 'perturbation_results.mat'), ...
    'cfg', 'eps_list', 'err_max', 'err_mean', 'ratio_max', 'n_violations_eps', ...
    'N_list', 'err_max_N', 'err_mean_N', 'n_violations_N', 'timestamp');

disp(exp1_table);
disp(exp2_table);



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