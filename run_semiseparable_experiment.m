function results_out = run_semiseparable_experiment(n_range, num_runs)
%RUN_SEMISEPARABLE_EXPERIMENT
% Full benchmark for semiseparable block-encoding circuits.
%
%
% Inputs:
%   n_range   : vector of n values (e.g. 2:6)
%   num_runs  : repetitions per configuration (default = 5)
%
% Output:
%   results_out : struct with all errors + metadata

    if nargin < 2
        num_runs = 5;
    end

    %% ============================================================
    % SETUP
    %% ============================================================

    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

    base_folder = fullfile('results', timestamp);
    latest_folder = fullfile('results', 'latest');

    mkdir(base_folder);
    mkdir(latest_folder);

    type_names = {'Random','Tridiagonal','Exponential'};
    colors = lines(3);
    markers = {'o','s','d'};

    all_errors = cell(3, length(n_range));
    %% ============================================================


    for ni = 1:length(n_range)

        n = n_range(ni);
        N = 2^n;

        for t = 1:3   % ALL GENERATORS

            all_errors{t,ni} = zeros(num_runs,1);

            for r = 1:num_runs

                rng(1000*n + 10*t + r);

                [u,v] = local_generate_uv(n, t);

                circ = build_semiseparable_circuit(u, v);
                U = circ.matrix;

                S = tril(u*v') + triu(v*u',1);

                S_tilde = 2*sqrt(N) * U(1:N,1:N);
                S
                S_tilde
                all_errors{t,ni}(r) = norm(S_tilde - S, 2);

            end
        end
    end

    %% ============================================================
    % SAVE DATA
    %% ============================================================

    save(fullfile(base_folder, 'results.mat'), ...
        'all_errors','n_range','timestamp');

    save(fullfile(latest_folder, 'results.mat'), ...
        'all_errors','n_range','timestamp');

    %% ============================================================
    % PLOTTING SETUP
    %% ============================================================

    %% ---------------- INDIVIDUAL PLOTS ----------------

    for t = 1:3

        fig = figure; hold on;

        for ni = 1:length(n_range)

            n = n_range(ni);
            y = all_errors{t,ni};

            scatter(n*ones(size(y)), y, 80, ...
                'filled', ...
                'MarkerFaceColor', colors(t,:), ...
                'Marker', markers{t});
        end

        title(type_names{t});
        xlabel('n');
        ylabel('||S - \tilde{S}||_2');
        grid on;

        name = ['scatter_type_' num2str(t)];

        saveas(fig, fullfile(base_folder, [name '.png']));
        saveas(fig, fullfile(base_folder, [name '.eps']));
        savefig(fig, fullfile(base_folder, [name '.fig']));

        saveas(fig, fullfile(latest_folder, [name '.png']));
        saveas(fig, fullfile(latest_folder, [name '.eps']));
        savefig(fig, fullfile(latest_folder, [name '.fig']));

        close(fig);
    end

    %% ---------------- OVERLAY PLOT ----------------

    fig = figure; hold on;

    for t = 1:3
        for ni = 1:length(n_range)

            n = n_range(ni);
            y = all_errors{t,ni};

            scatter(n*ones(size(y)), y, 80, ...
                'filled', ...
                'MarkerFaceColor', colors(t,:), ...
                'Marker', markers{t});
        end
    end

    xlabel('n');
    ylabel('||S - \tilde{S}||_2');
    title('Semiseparable Block-Encoding Error (All Generators)');
    grid on;
    legend(type_names);

    saveas(fig, fullfile(base_folder, 'scatter_overlay.png'));
    saveas(fig, fullfile(base_folder, 'scatter_overlay.eps'));
    savefig(fig, fullfile(base_folder, 'scatter_overlay.fig'));

    saveas(fig, fullfile(latest_folder, 'scatter_overlay.png'));
    saveas(fig, fullfile(latest_folder, 'scatter_overlay.eps'));
    savefig(fig, fullfile(latest_folder, 'scatter_overlay.fig'));

    close(fig);

    %% ============================================================
    % OUTPUT
    %% ============================================================

    results_out.all_errors = all_errors;
    results_out.n_range = n_range;
    results_out.timestamp = timestamp;
end

%% ============================================================
% LOCAL GENERATOR FUNCTION
%% ============================================================

function [u,v] = local_generate_uv(n, type)

    N = 2^n;

    switch type

        case 1
            u = randn(N,1); u = u/norm(u);
            v = randn(N,1); v = v/norm(v);

        case 2
            T = diag(randn(N,1)) + ...
                diag(randn(N-1,1),1) + ...
                diag(randn(N-1,1),-1);
            
            S=inv(T);
           
            [u,v] = uv(S);
            v=v';

        case 3
            gamma = cumsum(randn(N,1));

            u = exp(gamma);
            v = exp(-gamma);

            u = u/norm(u);
            v = v/norm(v);

        otherwise
            error('Invalid generator type');
    end
end
