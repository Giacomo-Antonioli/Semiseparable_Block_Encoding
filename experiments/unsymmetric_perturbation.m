% UNSYMMETRIC_PERTURBATION
%
% Tightness test for the UNSYMMETRIC semiseparable block encoding,
% restricted to the "Random" generator case only.
%
% Given (from generateUnsymm):
%
%   [S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);
%   u = u/nu;  v = v/nv;  x = x/nx;  y = y/ny;
%   S_exact = S / (sqrt(N)*(nu*nv + nx*ny));
%
% perturb the (normalized) generators:
%
%   u_tilde = u + delta_u,   ||delta_u|| = epsilon
%   v_tilde = v + delta_v,   ||delta_v|| = epsilon
%   x_tilde = x + delta_x,   ||delta_x|| = epsilon
%   y_tilde = y + delta_y,   ||delta_y|| = epsilon
%
% build the circuit with the PERTURBED generators (norms nu,nv,nx,ny are
% kept as originally generated, only the direction vectors are perturbed):
%
%   [circ,t1,t2] = build_unsymmetric_semiseparable_circuit( ...
%                       u_tilde,v_tilde,x_tilde,y_tilde,nu,nv,nx,ny);
%
% extract S_tilde column by column via simulation (as in
% demo_unsymmetric_semiseparable), and compute:
%
%   err = || S_exact - S_tilde ||_2
%
% Reference bound (Theorem thm:unsym, unsymmetric one-pair block encoding):
%
%   w_lo = nu*nv,  w_up = nx*ny            (norms of the ORIGINAL generators)
%   alpha = sqrt(N)*(w_lo + w_up)
%
%   || S - alpha*S_tilde_raw ||_2 <= 2*sqrt(w_lo^2 + w_up^2) * epsilon
%
% Since this script compares the ALREADY alpha-normalized quantities
% (S_exact = S/alpha  vs.  S_tilde = raw circuit output = S_tilde_raw),
% the bound is applied in the same normalized scale:
%
%   || S_exact - S_tilde ||_2 <= (2*sqrt(w_lo^2 + w_up^2)/alpha) * epsilon
%
% w_lo, w_up (hence alpha and the bound) depend on the actual random
% generators drawn in each trial, so the bound is computed PER TRIAL
% rather than being a single fixed constant.
%
% -------------------------------------------------------------------------
% SIMPLIFICATIONS REQUESTED:
%   - only the "Random" generator case (no Tridiagonal/Exponential loop)
%   - a single trial per (epsilon,N) point (cfg.num_trials = 1)
%   - the figures are live-updated with drawnow as points are computed,
%     but are NOT saved to disk
%   - for every point, the same disp(...) diagnostics as
%     demo_unsymmetric_semiseparable are printed (S_exact, S_tilde, error, and
%     the elementwise ratio S_exact./S_tilde)
% -------------------------------------------------------------------------

clear;
clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
root = setup_paths();

rng(42);


%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.n_fixed     = 3;              % N = 2^3 = 8 for the epsilon sweep
cfg.n_range     = 2:4;            % N sweep (2,4,8) for the N experiment

cfg.eps_powers  = -6:-2:-14;      % epsilon sweep exponents
cfg.eps_fixed   = 1e-6;           % epsilon held fixed for the N sweep

cfg.num_trials  = 1;              % single repetition, as requested

% NOTE: unlike the symmetric case, the bound is NOT a fixed constant.
% It depends on w_lo = nu*nv, w_up = nx*ny of the actual generators drawn
% each trial, so it is computed inside run_one_perturbation_trial_unsymm
% and returned together with the error.

cfg.verbose_matrices = true;      % print S_exact / S_tilde / err / ratio


%% ============================================================
% EXPERIMENT 1: ERROR VS EPSILON  (N fixed)
%% ============================================================

eps_list = 10.^cfg.eps_powers;
N_fixed  = 2^cfg.n_fixed;

err_eps   = zeros(length(eps_list),1);
bound_eps = zeros(length(eps_list),1);
viol_eps  = zeros(length(eps_list),1);

% ---- live figure 1 --------------------------------------------------
fig1 = figure('Name','Live: Error vs epsilon (Random, unsymmetric)', ...
              'NumberTitle','off');
ax1 = axes(fig1); hold(ax1,'on');

h1 = plot(ax1, NaN, NaN, '-o', 'LineWidth',1.5, 'MarkerSize',8, ...
          'DisplayName','Random (error)');
h1b = plot(ax1, NaN, NaN, 'k--s', 'LineWidth',1.5, 'MarkerSize',6, ...
          'DisplayName','$2\sqrt{w_{lo}^2+w_{up}^2}\,\varepsilon$');

set(ax1,'XScale','log','YScale','log');
xlabel(ax1,'$\varepsilon$','Interpreter','latex');
ylabel(ax1,'$\|S_{\mathrm{exact}}-\tilde S\|_2$','Interpreter','latex');
title(ax1, sprintf('Unsymmetric perturbation error (N=%d, Random)', N_fixed));
grid(ax1,'on');
legend(ax1,'Location','eastoutside','Interpreter','latex');
drawnow;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 1 (Random, unsymmetric): error vs epsilon\n');
fprintf('====================================================\n');

for ei = 1:length(eps_list)

    epsilon = eps_list(ei);

    [err,bound_val,bound_ok,u,v,x,y,u_tilde,v_tilde,x_tilde,y_tilde, ...
        S_exact,S_tilde,w_lo,w_up,alpha] = ...
        run_one_perturbation_trial_unsymm(N_fixed, epsilon);

    err_eps(ei)   = err;
    bound_eps(ei) = bound_val;
    viol_eps(ei)  = ~bound_ok;

    fprintf(['eps=%8.1e | error=%12.4e | bound=%12.4e | bound_ok=%d ', ...
        '| w_lo=%.4e w_up=%.4e alpha=%.4e\n'], ...
        epsilon, err, bound_val, bound_ok, w_lo, w_up, alpha);

    if cfg.verbose_matrices
        disp('S_exact:');
        disp(S_exact);
        disp('S_tilde:');
        disp(S_tilde);
        fprintf('Error ||S_exact - S_tilde||_2 = %.6e\n', err);
        fprintf('Theorem bound              = %.6e\n', bound_val);
        disp('S_exact ./ S_tilde:');
        disp(S_exact ./ S_tilde);
    end

    % live update
    set(h1,  'XData', eps_list(1:ei), 'YData', err_eps(1:ei)');
    set(h1b, 'XData', eps_list(1:ei), 'YData', bound_eps(1:ei)');
    drawnow;

end


%% ============================================================
% EXPERIMENT 2: ERROR VS N  (epsilon fixed)
%% ============================================================

N_list  = 2.^cfg.n_range;

err_N   = zeros(length(N_list),1);
bound_N = zeros(length(N_list),1);
viol_N  = zeros(length(N_list),1);

% ---- live figure 2 --------------------------------------------------
fig2 = figure('Name','Live: Error vs N (Random, unsymmetric)', ...
              'NumberTitle','off');
ax2 = axes(fig2); hold(ax2,'on');

h2 = plot(ax2, NaN, NaN, '-o', 'LineWidth',1.5, 'MarkerSize',8, ...
          'DisplayName','Random (error)');
h2b = plot(ax2, NaN, NaN, 'k--s', 'LineWidth',1.5, 'MarkerSize',6, ...
          'DisplayName','$2\sqrt{w_{lo}^2+w_{up}^2}\,\varepsilon$');

set(ax2,'XScale','log');
xticks(ax2, N_list);
xticklabels(ax2, arrayfun(@(x)sprintf('$%d$',2^x), ...
    cfg.n_range, 'UniformOutput',false));
set(ax2,'TickLabelInterpreter','latex');

xlabel(ax2,'$N$','Interpreter','latex');
ylabel(ax2,'$\|S_{\mathrm{exact}}-\tilde S\|_2$','Interpreter','latex');
title(ax2, sprintf('Unsymmetric perturbation error ($\\epsilon=%.1e$, Random)', ...
      cfg.eps_fixed), 'Interpreter','latex');
grid(ax2,'on');
legend(ax2,'Location','eastoutside','Interpreter','latex');
drawnow;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 2 (Random, unsymmetric): error vs N\n');
fprintf('====================================================\n');

for ni = 1:length(N_list)

    N = N_list(ni);

    [err,bound_val,bound_ok,u,v,x,y,u_tilde,v_tilde,x_tilde,y_tilde, ...
        S_exact,S_tilde,w_lo,w_up,alpha] = ...
        run_one_perturbation_trial_unsymm(N, cfg.eps_fixed);

    err_N(ni)   = err;
    bound_N(ni) = bound_val;
    viol_N(ni)  = ~bound_ok;

    fprintf(['N=%6d | error=%12.4e | bound=%12.4e | bound_ok=%d ', ...
        '| w_lo=%.4e w_up=%.4e alpha=%.4e\n'], ...
        N, err, bound_val, bound_ok, w_lo, w_up, alpha);

    if cfg.verbose_matrices
        disp('S_exact:');
        disp(S_exact);
        disp('S_tilde:');
        disp(S_tilde);
        fprintf('Error ||S_exact - S_tilde||_2 = %.6e\n', err);
        fprintf('Theorem bound              = %.6e\n', bound_val);
        disp('S_exact ./ S_tilde:');
        disp(S_exact ./ S_tilde);
    end

    % live update
    set(h2,  'XData', N_list(1:ni), 'YData', err_N(1:ni)');
    set(h2b, 'XData', N_list(1:ni), 'YData', bound_N(1:ni)');
    drawnow;

end


fprintf('\n');
fprintf('====================================================\n');
fprintf(' Unsymmetric (Random) perturbation tightness test done\n');
fprintf('====================================================\n');


%% ============================================================
% LATEX FIGURE GENERATION (pgfplots, same style as the reference doc)
%% ============================================================
%
% Produces a two-subfigure figure (error vs epsilon, error vs N), each
% with the single "Random" curve and the per-trial theorem bound curve
% 2*sqrt(w_lo^2+w_up^2)*epsilon (plotted through the actual
% computed bound values, since -- unlike the symmetric 2*sqrt(2)*epsilon
% case -- it is not a single fixed constant here).

tex_str = generate_latex_figure( ...
    eps_list, err_eps, bound_eps, N_fixed, ...
    N_list,   err_N,   bound_N,   cfg.eps_fixed);

tex_file = fullfile(root,'results','unsymmetric_perturbation','perturbation_unsymm_random.tex');
fid = fopen(tex_file,'w');
fprintf(fid,'%s', tex_str);
fclose(fid);

fprintf('\nLaTeX figure written to: %s\n\n', tex_file);
fprintf('%s\n', tex_str);





