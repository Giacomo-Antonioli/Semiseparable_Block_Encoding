% SYMMETRIC_PERTURBATION
%
% Tests the perturbation of the generators through the QCLAB block encoding.
%
% Given:
%
%   S(u,v) = tril(u*v') + triu(v*u',1)
%
% perturb:
%
%   u_tilde = u + delta_u
%   v_tilde = v + delta_v
%
% with:
%
%   ||delta_u|| = ||delta_v|| = epsilon
%
% and compute:
%
%   S_tilde = 2*sqrt(N)*U(1:N,1:N)
%
% where U is obtained from the QCLAB semiseparable circuit built with
% (u_tilde,v_tilde).
%
% Error:
%
%   ||S(u,v)-S_tilde||_2
%
% Reference:
%
%   2*sqrt(2)*epsilon
%
% -------------------------------------------------------------------------
%
% SAVING: instead of writing one .mat file per individual (eps/N, trial)
% run, all trial data for a given generator type is accumulated in memory
% and written ONCE, at the end of that type's iteration -- i.e. the save
% happens in the outermost loop (the type loop) rather than the innermost
% one. This applies to both Experiment 1 (error vs epsilon) and
% Experiment 2 (error vs N).
%
% LIVE PLOTTING: both figures are created before their respective
% experiment loops start and are updated in real time (drawnow) as each
% new epsilon/N data point is computed, so you can watch the curves build
% up on screen while the script runs. Figures are saved to disk once,
% after each experiment's loop completes.
%
% -------------------------------------------------------------------------

clear;
clc;
root = setup_paths();

rng(42);


%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.n_fixed          = 3;              % N=32 for epsilon sweep
cfg.n_range          = 1:3;           % N sweep

cfg.type_names       = {'Random','Tridiagonal','Exponential'};
cfg.types            = 1:numel(cfg.type_names);

cfg.eps_powers       = -6:-2:-14;
cfg.eps_fixed        = 1e-6;

cfg.num_trials       = 1;

cfg.bound_const      = 2*sqrt(2);

cfg.save_formats     = {'png','fig'};


n_types = numel(cfg.types);

colors = lines(n_types);
markers = {'o','s','d','^','v','p','h'};



%% ============================================================
% FOLDERS
%% ============================================================

timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');

base_folder   = fullfile(root,'results','symmetric_perturbation',timestamp);
latest_folder = fullfile(root,'results','symmetric_perturbation','latest');

plots_folder  = fullfile(base_folder,'plots');
latest_plots  = fullfile(latest_folder,'plots');


if ~exist(base_folder,'dir'); mkdir(base_folder); end
if ~exist(latest_folder,'dir'); mkdir(latest_folder); end

if ~exist(plots_folder,'dir'); mkdir(plots_folder); end
if ~exist(latest_plots,'dir'); mkdir(latest_plots); end


runs_folder = fullfile(base_folder,'runs');

if ~exist(runs_folder,'dir')
    mkdir(runs_folder);
end



%% ============================================================
% EXPERIMENT 1
% ERROR VS EPSILON
%% ============================================================


eps_list = 10.^cfg.eps_powers;


err_eps = zeros(length(eps_list),n_types);
viol_eps = zeros(length(eps_list),n_types);


%% ----------------------------------------------------------------
% LIVE FIGURE 1 -- created now, updated in real time during the loop
%% ----------------------------------------------------------------

fig1 = figure('Name','Live: Error vs epsilon','NumberTitle','off');
ax1 = axes(fig1); hold(ax1,'on');

h1 = gobjects(n_types,1);
for t = 1:n_types
    h1(t) = plot(ax1, NaN, NaN, '-o', ...
        'LineWidth',1.5, 'MarkerSize',8, ...
        'Color',colors(t,:), 'DisplayName',cfg.type_names{t});
end

plot(ax1, eps_list, cfg.bound_const*eps_list, 'k--', ...
    'LineWidth',1.5, 'DisplayName','$2\sqrt{2}\varepsilon$');

set(ax1,'XScale','log','YScale','log');
xlabel(ax1,'$\varepsilon$','Interpreter','latex');
ylabel(ax1,'$\|S-\tilde S\|_2$','Interpreter','latex');
title(ax1, sprintf('Perturbation error (N=%d)', 2^cfg.n_fixed));
grid(ax1,'on');
legend(ax1,'Location','eastoutside','Interpreter','latex');

drawnow;


fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 1: error versus epsilon\n');
fprintf('====================================================\n');


for t = 1:n_types


    fprintf('\nGenerator: %s\n',cfg.type_names{t});

    % ------------------------------------------------------------
    % Per-type accumulation storage (filled across all eps/trials
    % for this type, saved ONCE at the end of this outer iteration)
    % ------------------------------------------------------------
    type_u        = cell(length(eps_list), cfg.num_trials);
    type_v        = cell(length(eps_list), cfg.num_trials);
    type_u_tilde  = cell(length(eps_list), cfg.num_trials);
    type_v_tilde  = cell(length(eps_list), cfg.num_trials);
    type_S        = cell(length(eps_list), cfg.num_trials);
    type_S_tilde  = cell(length(eps_list), cfg.num_trials);
    type_err      = zeros(length(eps_list), cfg.num_trials);
    type_bound_ok = false(length(eps_list), cfg.num_trials);


    for ei = 1:length(eps_list)


        epsilon = eps_list(ei);


        errors = zeros(cfg.num_trials,1);



        for r = 1:cfg.num_trials


            [err,bound_ok,u,v,u_tilde,v_tilde,S,S_tilde] = ...
                run_one_perturbation_trial(...
                cfg.n_fixed,...
                t,...
                epsilon,...
                cfg.bound_const);



            errors(r)=err;

            type_u{ei,r}        = u;
            type_v{ei,r}        = v;
            type_u_tilde{ei,r}  = u_tilde;
            type_v_tilde{ei,r}  = v_tilde;
            type_S{ei,r}        = S;
            type_S_tilde{ei,r}  = S_tilde;
            type_err(ei,r)      = err;
            type_bound_ok(ei,r) = bound_ok;


            if ~bound_ok
                viol_eps(ei,t)=viol_eps(ei,t)+1;
            end


        end


        err_eps(ei,t)=max(errors);



        fprintf(...
        'eps=%8.1e | error=%12.4e | violations=%d/%d\n',...
        epsilon,...
        err_eps(ei,t),...
        viol_eps(ei,t),...
        cfg.num_trials);


        % --------------------------------------------------------
        % LIVE PLOT UPDATE: extend this type's curve with the point
        % just computed and redraw right now.
        % --------------------------------------------------------
        set(h1(t), 'XData', eps_list(1:ei), 'YData', err_eps(1:ei,t)');
        drawnow;


    end

    % ------------------------------------------------------------
    % SAVE: single consolidated file for this type, written once at
    % the end of the outermost (type) loop.
    % ------------------------------------------------------------
    save(fullfile(runs_folder, sprintf('type%d_eps_sweep.mat', t)), ...
        'type_u','type_v','type_u_tilde','type_v_tilde', ...
        'type_S','type_S_tilde','type_err','type_bound_ok', ...
        'eps_list','t','-v7.3');

end


save_figure_pair(fig1,...
    plots_folder,...
    latest_plots,...
    'error_vs_epsilon',...
    cfg.save_formats);



%% ============================================================
% EXPERIMENT 2
% ERROR VS N
%% ============================================================


N_list = 2.^cfg.n_range;


err_N = zeros(length(N_list),n_types);
viol_N = zeros(length(N_list),n_types);


%% ----------------------------------------------------------------
% LIVE FIGURE 2 -- created now, updated in real time during the loop
%% ----------------------------------------------------------------

fig2 = figure('Name','Live: Error vs N','NumberTitle','off');
ax2 = axes(fig2); hold(ax2,'on');

h2 = gobjects(n_types,1);
for t = 1:n_types
    h2(t) = plot(ax2, NaN, NaN, '-o', ...
        'LineWidth',1.5, 'MarkerSize',8, ...
        'Color',colors(t,:), 'DisplayName',cfg.type_names{t});
end

yline(ax2, cfg.bound_const*cfg.eps_fixed, 'k--', ...
    'LineWidth',1.5, 'DisplayName','$2\sqrt{2}\varepsilon$');

set(ax2,'XScale','log');
xticks(ax2, N_list);
xticklabels(ax2, arrayfun(@(x)sprintf('$%d$',2^x), ...
    cfg.n_range, 'UniformOutput',false));
set(ax2,'TickLabelInterpreter','latex');

xlabel(ax2,'$N$','Interpreter','latex');
ylabel(ax2,'$\|S-\tilde S\|_2$','Interpreter','latex');
title(ax2, sprintf('Perturbation error ($\\epsilon=%.1e$)', ...
      cfg.eps_fixed), 'Interpreter','latex');
grid(ax2,'on');
legend(ax2,'Location','eastoutside','Interpreter','latex');

drawnow;  


fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 2: error versus N\n');
fprintf('====================================================\n');



for t=1:n_types


    fprintf('\nGenerator: %s\n',cfg.type_names{t});

    % ------------------------------------------------------------
    % Per-type accumulation storage (filled across all N/trials for
    % this type, saved ONCE at the end of this outer iteration)
    % ------------------------------------------------------------
    type_u        = cell(length(N_list), cfg.num_trials);
    type_v        = cell(length(N_list), cfg.num_trials);
    type_u_tilde  = cell(length(N_list), cfg.num_trials);
    type_v_tilde  = cell(length(N_list), cfg.num_trials);
    type_S        = cell(length(N_list), cfg.num_trials);
    type_S_tilde  = cell(length(N_list), cfg.num_trials);
    type_err      = zeros(length(N_list), cfg.num_trials);
    type_bound_ok = false(length(N_list), cfg.num_trials);


    for ni=1:length(N_list)


        n = cfg.n_range(ni);


        errors=zeros(cfg.num_trials,1);



        for r=1:cfg.num_trials


            [err,bound_ok,u,v,u_tilde,v_tilde,S,S_tilde] = ...
                run_one_perturbation_trial(...
                n,...
                t,...
                cfg.eps_fixed,...
                cfg.bound_const);



            errors(r)=err;

            type_u{ni,r}        = u;
            type_v{ni,r}        = v;
            type_u_tilde{ni,r}  = u_tilde;
            type_v_tilde{ni,r}  = v_tilde;
            type_S{ni,r}        = S;
            type_S_tilde{ni,r}  = S_tilde;
            type_err(ni,r)      = err;
            type_bound_ok(ni,r) = bound_ok;



            if ~bound_ok
                viol_N(ni,t)=viol_N(ni,t)+1;
            end


        end



        err_N(ni,t)=max(errors);



        fprintf(...
        'N=%6d | error=%12.4e | violations=%d/%d\n',...
        N_list(ni),...
        err_N(ni,t),...
        viol_N(ni,t),...
        cfg.num_trials);


        % --------------------------------------------------------
        % LIVE PLOT UPDATE: extend this type's curve with the point
        % just computed and redraw right now.
        % --------------------------------------------------------
        set(h2(t), 'XData', N_list(1:ni), 'YData', err_N(1:ni,t)');
        drawnow;


    end

    % ------------------------------------------------------------
    % SAVE: single consolidated file for this type, written once at
    % the end of the outermost (type) loop.
    % ------------------------------------------------------------
    save(fullfile(runs_folder, sprintf('type%d_N_sweep.mat', t)), ...
        'type_u','type_v','type_u_tilde','type_v_tilde', ...
        'type_S','type_S_tilde','type_err','type_bound_ok', ...
        'N_list','t','-v7.3');

end


save_figure_pair(fig2,...
    plots_folder,...
    latest_plots,...
    'error_vs_N',...
    cfg.save_formats);




%% ============================================================
% CSV SUMMARY
%% ============================================================


summary_eps = [];


for t = 1:n_types

    for ei = 1:length(eps_list)


        summary_eps = [summary_eps;
            t,...
            eps_list(ei),...
            err_eps(ei,t),...
            viol_eps(ei,t)];


    end

end



T_eps = array2table(summary_eps,...
    'VariableNames',...
    {'Type','epsilon','MaxError','Violations'});


writetable(T_eps,...
    fullfile(base_folder,...
    'epsilon_summary.csv'));

writetable(T_eps,...
    fullfile(latest_folder,...
    'epsilon_summary.csv'));




summary_N = [];


for t = 1:n_types

    for ni = 1:length(N_list)


        summary_N = [summary_N;
            t,...
            N_list(ni),...
            err_N(ni,t),...
            viol_N(ni,t)];


    end

end



T_N = array2table(summary_N,...
    'VariableNames',...
    {'Type','N','MaxError','Violations'});



writetable(T_N,...
    fullfile(base_folder,...
    'N_summary.csv'));

writetable(T_N,...
    fullfile(latest_folder,...
    'N_summary.csv'));




%% ============================================================
% SAVE COMPLETE DATA
%% ============================================================


save(fullfile(base_folder,...
    'perturbation_qclab_results.mat'),...
    'cfg',...
    'eps_list',...
    'err_eps',...
    'viol_eps',...
    'N_list',...
    'err_N',...
    'viol_N',...
    'timestamp');



save(fullfile(latest_folder,...
    'perturbation_qclab_results.mat'),...
    'cfg',...
    'eps_list',...
    'err_eps',...
    'viol_eps',...
    'N_list',...
    'err_N',...
    'viol_N',...
    'timestamp');



disp(T_eps);
disp(T_N);



fprintf('\n');
fprintf('====================================================\n');
fprintf(' Perturbation benchmark completed\n');
fprintf(' Results saved in:\n %s\n',base_folder);
fprintf('====================================================\n');




%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================


