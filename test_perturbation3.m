% TEST_PERTURBATION_QCLAB
%
% Perturbation benchmark for QCLAB semiseparable block encoding.
%
% Computation unchanged:
%
%   S(u,v) = tril(u*v') + triu(v*u',1)
%
%   u_tilde = u + delta_u
%   v_tilde = v + delta_v
%
%   S_tilde = 2*sqrt(N)*U(1:N,1:N)
%
% Error:
%
%   ||S(u,v)-S_tilde||_2
%
% Bound:
%
%   2*sqrt(2)*epsilon
%
% Saving/plotting:
%   - incremental run files
%   - disk based reconstruction of final plots
%   - live curves during execution
%
% -------------------------------------------------------------------------

clear;
clc;
rng(42);

%% ============================================================
% CONFIGURATION
%% ============================================================

cfg.n_fixed          = 3;
cfg.n_range          = 2:3;
cfg.type_names       = {'Random','Tridiagonal','Exponential'};
cfg.types            = 1:numel(cfg.type_names);

cfg.eps_powers       = -6:-2:-14;
cfg.eps_fixed        = 1e-6;

cfg.num_trials       = 1;
cfg.bound_const      = 2*sqrt(2);

cfg.save_formats     = {'png','fig'};

cfg.enable_saving       = true;
cfg.generate_plots_only = false;
cfg.load_timestamp      = '';

n_types = numel(cfg.type_names);

colors = lines(n_types);

markers = {'o','s','d','^','v','p','h'};
markers = markers(mod(0:n_types-1,numel(markers))+1);


%% ============================================================
% EXECUTION MODE
%% ============================================================

if ~cfg.generate_plots_only

%% ============================================================
% FOLDERS
%% ============================================================

timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');

base_folder   = fullfile('results_perturbation',timestamp);
latest_folder = fullfile('results_perturbation','latest');
plots_folder  = fullfile(base_folder,'plots');

if cfg.enable_saving
    if ~exist(base_folder,'dir'); mkdir(base_folder); end
    if ~exist(latest_folder,'dir'); mkdir(latest_folder); end
    if ~exist(plots_folder,'dir'); mkdir(plots_folder); end
end


%% ============================================================
% LIVE FIGURE 1
% ERROR VS EPSILON
%% ============================================================

eps_list = 10.^cfg.eps_powers;

fig1 = figure('Name','Live: Error vs epsilon','NumberTitle','off');
ax1 = axes(fig1);
hold(ax1,'on');

h_eps = gobjects(n_types,1);

for t = 1:n_types
    h_eps(t)=plot(ax1,NaN,NaN,'-o',...
        'LineWidth',1.5,...
        'MarkerSize',8,...
        'Color',colors(t,:),...
        'DisplayName',cfg.type_names{t});
end

plot(ax1,...
    eps_list,...
    cfg.bound_const*eps_list,...
    'k--',...
    'LineWidth',1.5,...
    'DisplayName','$2\sqrt{2}\varepsilon$');


set(ax1,'XScale','log','YScale','log');

eps_ticks = sort(eps_list);

xticks(ax1,eps_ticks);

xticklabels(ax1,...
    arrayfun(@(x)sprintf('$10^{%d}$',x),...
    flip(cfg.eps_powers),...
    'UniformOutput',false));
set(ax1,'TickLabelInterpreter','latex');

xlabel(ax1,'$\varepsilon$','Interpreter','latex');
ylabel(ax1,'$\|S-\tilde S_{\mathrm{QCLAB}}\|_2$','Interpreter','latex');

title(ax1,...
    sprintf('QCLAB perturbation error ($N=%d$)',2^cfg.n_fixed),...
    'Interpreter','latex');

grid(ax1,'on');

legend(ax1,...
    'Location','eastoutside',...
    'Interpreter','latex');

drawnow;


%% ============================================================
% LIVE FIGURE 2
% ERROR VS N
%% ============================================================

N_list = 2.^cfg.n_range;

fig2 = figure('Name','Live: Error vs N','NumberTitle','off');
ax2 = axes(fig2);
hold(ax2,'on');

h_N = gobjects(n_types,1);

for t = 1:n_types
    h_N(t)=plot(ax2,NaN,NaN,'-o',...
        'LineWidth',1.5,...
        'MarkerSize',8,...
        'Color',colors(t,:),...
        'DisplayName',cfg.type_names{t});
end

yline(ax2,...
    cfg.bound_const*cfg.eps_fixed,...
    'k--',...
    'LineWidth',1.5,...
    'DisplayName','$2\sqrt{2}\varepsilon$');


set(ax2,'XScale','log');

xticks(ax2,N_list);

xticklabels(ax2,...
    arrayfun(@(x)sprintf('$2^{%d}$',x),...
    cfg.n_range,'UniformOutput',false));

set(ax2,'TickLabelInterpreter','latex');

xlabel(ax2,'$N$','Interpreter','latex');

ylabel(ax2,...
    '$\|S-\tilde S_{\mathrm{QCLAB}}\|_2$',...
    'Interpreter','latex');

title(ax2,...
    sprintf('QCLAB perturbation error ($\\epsilon=%.1e$)',...
    cfg.eps_fixed),...
    'Interpreter','latex');

grid(ax2,'on');

legend(ax2,...
    'Location','eastoutside',...
    'Interpreter','latex');

drawnow;
%% ============================================================
% EXPERIMENT 1
% ERROR VS EPSILON
%% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 1: error versus epsilon\n');
fprintf('====================================================\n');

num_eps = length(eps_list);
total_tasks = num_eps*n_types*cfg.num_trials;
done_tasks = 0;

for ei = 1:num_eps

    epsilon = eps_list(ei);
    eps_block_data = cell(n_types,1);

    for t = 1:n_types

        errors = zeros(cfg.num_trials,1);
        bound_ok_all = false(cfg.num_trials,1);

        u_cell = cell(cfg.num_trials,1);
        v_cell = cell(cfg.num_trials,1);
        u_tilde_cell = cell(cfg.num_trials,1);
        v_tilde_cell = cell(cfg.num_trials,1);
        S_cell = cell(cfg.num_trials,1);
        Stilde_cell = cell(cfg.num_trials,1);

        for r = 1:cfg.num_trials

            [err,bound_ok,...
                u,v,u_tilde,v_tilde,S,S_tilde] = ...
                run_one_perturbation_trial(...
                cfg.n_fixed,...
                t,...
                epsilon,...
                cfg.bound_const);

            errors(r) = err;
            bound_ok_all(r) = bound_ok;

            u_cell{r} = u;
            v_cell{r} = v;
            u_tilde_cell{r} = u_tilde;
            v_tilde_cell{r} = v_tilde;
            S_cell{r} = S;
            Stilde_cell{r} = S_tilde;

            done_tasks = done_tasks + 1;

            print_progress(...
                done_tasks,...
                total_tasks,...
                ei,...
                epsilon,...
                t,...
                r,...
                cfg.num_trials);
        end

        eps_block_data{t}.errors = errors;
        eps_block_data{t}.bound_ok = bound_ok_all;

        eps_block_data{t}.u = u_cell;
        eps_block_data{t}.v = v_cell;
        eps_block_data{t}.u_tilde = u_tilde_cell;
        eps_block_data{t}.v_tilde = v_tilde_cell;
        eps_block_data{t}.S = S_cell;
        eps_block_data{t}.S_tilde = Stilde_cell;

        set(h_eps(t),...
            'XData',eps_list(1:ei),...
            'YData',[get(h_eps(t),'YData'),max(errors)]);

        drawnow;
    end


    %% --------------------------------------------------------
    % INCREMENTAL SAVE
    %% --------------------------------------------------------

    if cfg.enable_saving

        runs_folder = fullfile(base_folder,'eps_runs');

        if ~exist(runs_folder,'dir')
            mkdir(runs_folder);
        end

        for t = 1:n_types

            data = eps_block_data{t};

            for r = 1:cfg.num_trials

                u = data.u{r};
                v = data.v{r};

                u_tilde = data.u_tilde{r};
                v_tilde = data.v_tilde{r};

                S = data.S{r};
                S_tilde = data.S_tilde{r};

                err = data.errors(r);
                bound_ok = data.bound_ok(r);

                run_file = fullfile(...
                    runs_folder,...
                    sprintf('eps_%d_type_%d_run_%d.mat',...
                    ei,t,r));

                save(run_file,...
                    'u','v',...
                    'u_tilde','v_tilde',...
                    'S','S_tilde',...
                    'err','bound_ok',...
                    'epsilon','t','r',...
                    '-v7.3');
            end
        end
    end
end


save_figure_pair(...
    fig1,...
    plots_folder,...
    latest_folder,...
    'error_vs_epsilon',...
    cfg.save_formats);
%% ============================================================
% EXPERIMENT 2
% ERROR VS N
%% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' Experiment 2: error versus N\n');
fprintf('====================================================\n');

num_N = length(N_list);

total_tasks = num_N*n_types*cfg.num_trials;
done_tasks = 0;

for ni = 1:num_N

    n = cfg.n_range(ni);

    N = N_list(ni);

    N_block_data = cell(n_types,1);

    for t = 1:n_types

        errors = zeros(cfg.num_trials,1);
        bound_ok_all = false(cfg.num_trials,1);

        u_cell = cell(cfg.num_trials,1);
        v_cell = cell(cfg.num_trials,1);
        u_tilde_cell = cell(cfg.num_trials,1);
        v_tilde_cell = cell(cfg.num_trials,1);
        S_cell = cell(cfg.num_trials,1);
        Stilde_cell = cell(cfg.num_trials,1);


        for r = 1:cfg.num_trials

            [err,bound_ok,...
                u,v,u_tilde,v_tilde,S,S_tilde] = ...
                run_one_perturbation_trial(...
                n,...
                t,...
                cfg.eps_fixed,...
                cfg.bound_const);


            errors(r) = err;
            bound_ok_all(r) = bound_ok;

            u_cell{r} = u;
            v_cell{r} = v;

            u_tilde_cell{r} = u_tilde;
            v_tilde_cell{r} = v_tilde;

            S_cell{r} = S;
            Stilde_cell{r} = S_tilde;


            done_tasks = done_tasks + 1;

            print_progress(...
                done_tasks,...
                total_tasks,...
                n,...
                N,...
                t,...
                r,...
                cfg.num_trials);

        end


        N_block_data{t}.errors = errors;
        N_block_data{t}.bound_ok = bound_ok_all;

        N_block_data{t}.u = u_cell;
        N_block_data{t}.v = v_cell;

        N_block_data{t}.u_tilde = u_tilde_cell;
        N_block_data{t}.v_tilde = v_tilde_cell;

        N_block_data{t}.S = S_cell;
        N_block_data{t}.S_tilde = Stilde_cell;


        set(h_N(t),...
            'XData',N_list(1:ni),...
            'YData',[get(h_N(t),'YData'),max(errors)]);


        drawnow;

    end


    %% --------------------------------------------------------
    % INCREMENTAL SAVE
    %% --------------------------------------------------------

    if cfg.enable_saving

        runs_folder = fullfile(base_folder,'N_runs');

        if ~exist(runs_folder,'dir')
            mkdir(runs_folder);
        end


        for t = 1:n_types

            data = N_block_data{t};


            for r = 1:cfg.num_trials

                u = data.u{r};
                v = data.v{r};

                u_tilde = data.u_tilde{r};
                v_tilde = data.v_tilde{r};

                S = data.S{r};
                S_tilde = data.S_tilde{r};

                err = data.errors(r);
                bound_ok = data.bound_ok(r);


                run_file = fullfile(...
                    runs_folder,...
                    sprintf('N_%d_type_%d_run_%d.mat',...
                    N,t,r));


                save(run_file,...
                    'u','v',...
                    'u_tilde','v_tilde',...
                    'S','S_tilde',...
                    'err','bound_ok',...
                    'n','N','t','r',...
                    '-v7.3');

            end
        end
    end
end


save_figure_pair(...
    fig2,...
    plots_folder,...
    latest_folder,...
    'error_vs_N',...
    cfg.save_formats);
%% ============================================================
% FINAL PLOTS FROM SAVED DATA
%% ============================================================

if cfg.enable_saving || cfg.generate_plots_only

    if cfg.generate_plots_only

        if isempty(cfg.load_timestamp)
            base_folder = fullfile('results_perturbation','latest');
        else
            base_folder = fullfile('results_perturbation',cfg.load_timestamp);
        end

        plots_folder = fullfile(base_folder,'plots');
        latest_folder = fullfile('results_perturbation','latest');

    end

    generate_overall_plots(...
        base_folder,...
        plots_folder,...
        latest_folder,...
        cfg,...
        colors,...
        markers);

end


%% ============================================================
% CSV SUMMARIES
%% ============================================================

if cfg.enable_saving && ~cfg.generate_plots_only

    generate_csv_summaries(...
        base_folder,...
        latest_folder,...
        cfg);

end

end

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================


function generate_overall_plots(base_folder,plots_folder,latest_folder,cfg,colors,markers)

fprintf('Generating final plots from saved data...\n');

n_types = numel(cfg.type_names);


%% ------------------------------------------------------------
% EPSILON PLOT
%% ------------------------------------------------------------

fig = figure(...
    'Name','Final: Error vs epsilon',...
    'NumberTitle','off',...
    'Visible','off');

ax = axes(fig);
hold(ax,'on');


for t = 1:n_types

    eps_vals = [];
    err_vals = [];

    files = dir(fullfile(...
        base_folder,...
        'eps_runs',...
        sprintf('eps_*_type_%d_run_*.mat',t)));


    for k = 1:length(files)

        data = load(fullfile(...
            files(k).folder,...
            files(k).name));

        eps_vals(end+1) = data.epsilon;
        err_vals(end+1) = data.err;

    end


    [eps_vals,idx] = sort(eps_vals);
    err_vals = err_vals(idx);


    plot(ax,...
        eps_vals,...
        err_vals,...
        '-o',...
        'LineWidth',1.5,...
        'MarkerSize',8,...
        'Color',colors(t,:),...
        'DisplayName',cfg.type_names{t});

end


eps_list = 10.^cfg.eps_powers;


plot(ax,...
    eps_list,...
    cfg.bound_const*eps_list,...
    'k--',...
    'LineWidth',1.5,...
    'DisplayName','$2\sqrt{2}\varepsilon$');


set(ax,'XScale','log','YScale','log');

xticks(ax,eps_list);

xticklabels(ax,...
    arrayfun(@(x)sprintf('$10^{%d}$',x),...
    cfg.eps_powers,...
    'UniformOutput',false));

set(ax,'TickLabelInterpreter','latex');

xlabel(ax,'$\varepsilon$','Interpreter','latex');

ylabel(ax,...
    '$\|S-\tilde S_{\mathrm{QCLAB}}\|_2$',...
    'Interpreter','latex');


grid(ax,'on');

legend(ax,...
    'Location','eastoutside',...
    'Interpreter','latex');


save_figure_pair(...
    fig,...
    plots_folder,...
    latest_folder,...
    'final_error_vs_epsilon',...
    cfg.save_formats);


close(fig);



%% ------------------------------------------------------------
% N PLOT
%% ------------------------------------------------------------

fig = figure(...
    'Name','Final: Error vs N',...
    'NumberTitle','off',...
    'Visible','off');


ax = axes(fig);
hold(ax,'on');


for t = 1:n_types

    N_vals = [];
    err_vals = [];


    files = dir(fullfile(...
        base_folder,...
        'N_runs',...
        sprintf('N_*_type_%d_run_*.mat',t)));


    for k = 1:length(files)

        data = load(fullfile(...
            files(k).folder,...
            files(k).name));

        N_vals(end+1) = data.N;
        err_vals(end+1) = data.err;

    end


    [N_vals,idx] = sort(N_vals);
    err_vals = err_vals(idx);


    plot(ax,...
        N_vals,...
        err_vals,...
        '-o',...
        'LineWidth',1.5,...
        'MarkerSize',8,...
        'Color',colors(t,:),...
        'DisplayName',cfg.type_names{t});

end


N_list = 2.^cfg.n_range;


yline(ax,...
    cfg.bound_const*cfg.eps_fixed,...
    'k--',...
    'LineWidth',1.5,...
    'DisplayName','$2\sqrt{2}\varepsilon$');


set(ax,'XScale','log');


xticks(ax,N_list);

xticklabels(ax,...
    arrayfun(@(x)sprintf('$2^{%d}$',x),...
    cfg.n_range,...
    'UniformOutput',false));


set(ax,'TickLabelInterpreter','latex');


xlabel(ax,'$N$','Interpreter','latex');

ylabel(ax,...
    '$\|S-\tilde S_{\mathrm{QCLAB}}\|_2$',...
    'Interpreter','latex');


grid(ax,'on');

legend(ax,...
    'Location','eastoutside',...
    'Interpreter','latex');


save_figure_pair(...
    fig,...
    plots_folder,...
    latest_folder,...
    'final_error_vs_N',...
    cfg.save_formats);


close(fig);

fprintf('Final plots completed.\n');

end



function generate_csv_summaries(base_folder,latest_folder,cfg)

rows = {};

n_types = numel(cfg.type_names);


%% epsilon summary

for t = 1:n_types

    files = dir(fullfile(...
        base_folder,...
        'eps_runs',...
        sprintf('eps_*_type_%d_run_*.mat',t)));

    for k = 1:length(files)

        d = load(fullfile(...
            files(k).folder,...
            files(k).name));

        rows(end+1,:) = {...
            'epsilon',...
            cfg.type_names{t},...
            d.epsilon,...
            d.err,...
            d.bound_ok};

    end
end


T = cell2table(rows,...
    'VariableNames',...
    {'Experiment','Type','Value','Error','BoundOK'});


writetable(T,...
    fullfile(base_folder,'summary.csv'));

writetable(T,...
    fullfile(latest_folder,'summary.csv'));

end



function print_progress(done,total,a,b,t,r,num_runs)

percent = 100*done/total;

bar_len = 30;

filled = round(bar_len*done/total);

bar_str = [...
    repmat('#',1,filled),...
    repmat('-',1,bar_len-filled)];


fprintf('\r[%s] %6.2f%% | step=%d | value=%g | type=%d | run=%d/%d',...
    bar_str,...
    percent,...
    a,...
    b,...
    t,...
    r,...
    num_runs);


if done==total
    fprintf('\n');
end

end



function save_figure_pair(fig,base_folder,latest_folder,name,formats)


for i = 1:numel(formats)

    fmt = formats{i};

    switch fmt

        case 'fig'

            savefig(fig,...
                fullfile(base_folder,[name '.fig']));

            savefig(fig,...
                fullfile(latest_folder,[name '.fig']));

        otherwise

            saveas(fig,...
                fullfile(base_folder,[name '.' fmt]));

            saveas(fig,...
                fullfile(latest_folder,[name '.' fmt]));

    end

end

end

function [err,bound_ok,u,v,u_tilde,v_tilde,S,S_tilde] = ...
    run_one_perturbation_trial(n,type,epsilon,bound_const)

N = 2^n;

% Generate original semiseparable generators
[u,v] = local_generate_uv(n,type);

u = u(:);
v = v(:);

% Normalize generators
u = u/norm(u);
v = v/norm(v);

% Random perturbations
delta_u = randn(N,1);
delta_u = epsilon*delta_u/norm(delta_u);

delta_v = randn(N,1);
delta_v = epsilon*delta_v/norm(delta_v);

u_tilde = u + delta_u;
v_tilde = v + delta_v;


% Original semiseparable matrix
S = tril(u*v') + triu(v*u',1);


% QCLAB circuit generated from perturbed generators
circ = build_semiseparable_circuit(u_tilde,v_tilde);

U = circ.matrix;


% Extract encoded matrix
S_tilde = 2*sqrt(N)*U(1:N,1:N);


% Perturbation error with respect to the circuit
err = norm(S-S_tilde,2);


% Bound verification
bound_ok = err < bound_const*epsilon;

end
