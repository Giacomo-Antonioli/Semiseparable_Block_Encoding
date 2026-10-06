% plot_saved_semiseparable_results.m
%
% This script loads saved semiseparable generators and error results.
%
% It can work in two modes:
%
%   recompute_results = true
%
%       Load the generators using SemiseparableGeneratorLoader,
%       reconstruct S, rebuild the circuit, extract S_tilde again,
%       and recompute:
%
%           || S - alpha S_tilde ||_2
%
%       and
%
%           || S^{-1} - (alpha S_tilde)^{-1} ||_2.
%
%
%   recompute_results = false
%
%       Load the already computed errors from:
%
%           semiseparable_single_spectral_errors.txt
%
%       and plot them directly.
%
%
% Required files in the selected folder:
%
%   semiseparable_single_generators.txt
%   semiseparable_single_spectral_errors.txt
%
% Required class:
%
%   SemiseparableGeneratorLoader.m
%
% Saved plots:
%
%   loaded_forward_spectral_error.png
%   loaded_inverse_spectral_error.png

clear;
clc;
close all;

% -------------------------------------------------------------------------
% User flag
% -------------------------------------------------------------------------

recompute_results = true;

% Set to true if you want to recompute the errors from the generators.
% Set to false if you only want to reload the saved errors and plot them.

% -------------------------------------------------------------------------
% Select results folder using GUI
% -------------------------------------------------------------------------

results_folder = uigetdir(pwd, 'Select the results folder');

if isequal(results_folder, 0)
    error('No folder selected.');
end

generators_filename = fullfile(results_folder, ...
    'semiseparable_single_generators.txt');

results_filename = fullfile(results_folder, ...
    'semiseparable_single_spectral_errors.txt');

if ~isfile(generators_filename)
    error('Generator file not found:\n%s', generators_filename);
end

if ~isfile(results_filename)
    error('Results file not found:\n%s', results_filename);
end

fprintf('Selected results folder:\n%s\n\n', results_folder);

% -------------------------------------------------------------------------
% Load generators
% -------------------------------------------------------------------------

loader = SemiseparableGeneratorLoader(generators_filename);

fprintf('Loaded %d generator configurations.\n\n', ...
    loader.numGenerators());

% -------------------------------------------------------------------------
% Either recompute or load existing results
% -------------------------------------------------------------------------

if recompute_results

    fprintf('Recomputing errors from saved generators.\n\n');

    all_N = [];
    all_test_index = [];
    all_alpha = [];
    all_error_norm2 = [];
    all_inverse_error_norm2 = [];
    all_rcond_S = [];
    all_rcond_S_extracted = [];

    for q = 1:loader.numGenerators()

        [x, y, N, test_index] = loader.get(q);

        n = log2(N);
        alpha = 2^(n/2 + 1);

        fprintf('Recomputing index %d / %d: N = %d, test = %d\n', ...
            q, loader.numGenerators(), N, test_index);

        % -----------------------------------------------------------------
        % Original semiseparable matrix
        % -----------------------------------------------------------------

        S = tril(x * y') + triu(y * x', 1);

        % -----------------------------------------------------------------
        % Build circuit
        % -----------------------------------------------------------------

        circ = build_semiseparable_circuit(x, y);

        % Depending on your QCLAB version, use one of:
        %
        %   nbQubits = circ.nbQubits;
        %   nbQubits = circ.nbQubits();
        %

        nbQubits = circ.nbQubits;

        dim = 2^nbQubits;

        % -----------------------------------------------------------------
        % Reconstruct first N columns by simulation
        % -----------------------------------------------------------------

        first_N_cols = zeros(dim, N);

        progress_msg = '';

        for j = 1:N

            fprintf(repmat('\b', 1, length(progress_msg)));

            progress_msg = sprintf('    Simulating column %d / %d', j, N);
            fprintf('%s', progress_msg);

            basis_index = j - 1;

            input_state = dec2bin(basis_index, nbQubits);

            sim = circ.simulate(input_state);

            state = sim.states();

            first_N_cols(:, j) = state(:);

        end

        fprintf('\n');

        % -----------------------------------------------------------------
        % Extract block and scale it
        % -----------------------------------------------------------------

        S_tilde = full(first_N_cols(1:N, 1:N));

        S_extracted = alpha * S_tilde;

        % -----------------------------------------------------------------
        % Forward error
        % -----------------------------------------------------------------

        error_norm2 = norm(S - S_extracted, 2);

        % -----------------------------------------------------------------
        % Inverse error
        % -----------------------------------------------------------------

        rcond_S = rcond(S);
        rcond_S_extracted = rcond(S_extracted);

        if rcond_S == 0 || rcond_S_extracted == 0

            inverse_error_norm2 = NaN;

            warning(['Singular matrix detected for N = %d, test = %d. ', ...
                     'Inverse error set to NaN.'], ...
                     N, test_index);

        else

            S_inv = S \ eye(N);
            S_extracted_inv = S_extracted \ eye(N);

            inverse_error_norm2 = norm(S_inv - S_extracted_inv, 2);

        end

        % -----------------------------------------------------------------
        % Store recomputed data
        % -----------------------------------------------------------------

        all_N(end + 1, 1) = N;
        all_test_index(end + 1, 1) = test_index;
        all_alpha(end + 1, 1) = alpha;
        all_error_norm2(end + 1, 1) = error_norm2;
        all_inverse_error_norm2(end + 1, 1) = inverse_error_norm2;
        all_rcond_S(end + 1, 1) = rcond_S;
        all_rcond_S_extracted(end + 1, 1) = rcond_S_extracted;

    end

    % ---------------------------------------------------------------------
    % Save recomputed results to a separate file
    % ---------------------------------------------------------------------

    recomputed_results_filename = fullfile(results_folder, ...
        'semiseparable_single_spectral_errors_recomputed.txt');

    fid = fopen(recomputed_results_filename, 'w');

    fprintf(fid, ['N test_index alpha error_norm2 ', ...
                  'inverse_error_norm2 rcond_S rcond_S_extracted\n']);

    for q = 1:length(all_N)

        fprintf(fid, '%d %d %.16e %.16e %.16e %.16e %.16e\n', ...
            all_N(q), ...
            all_test_index(q), ...
            all_alpha(q), ...
            all_error_norm2(q), ...
            all_inverse_error_norm2(q), ...
            all_rcond_S(q), ...
            all_rcond_S_extracted(q));

    end

    fclose(fid);

    fprintf('\nSaved recomputed results to:\n%s\n\n', ...
        recomputed_results_filename);

else

    fprintf('Loading errors from existing results file.\n\n');

    T = readtable(results_filename, ...
        'FileType', 'text', ...
        'Delimiter', ' ', ...
        'MultipleDelimsAsOne', true);

    required_columns = {'N', 'test_index', 'alpha', ...
                        'error_norm2', 'inverse_error_norm2'};

    for c = 1:numel(required_columns)

        if ~ismember(required_columns{c}, T.Properties.VariableNames)
            error('Missing required column in results file: %s', ...
                required_columns{c});
        end

    end

    all_N = T.N;
    all_test_index = T.test_index;
    all_alpha = T.alpha;
    all_error_norm2 = T.error_norm2;
    all_inverse_error_norm2 = T.inverse_error_norm2;

    if ismember('rcond_S', T.Properties.VariableNames)
        all_rcond_S = T.rcond_S;
    else
        all_rcond_S = NaN(size(all_N));
    end

    if ismember('rcond_S_extracted', T.Properties.VariableNames)
        all_rcond_S_extracted = T.rcond_S_extracted;
    else
        all_rcond_S_extracted = NaN(size(all_N));
    end

end

% -------------------------------------------------------------------------
% Sort data by N and test_index
% -------------------------------------------------------------------------

data_table = table( ...
    all_N, ...
    all_test_index, ...
    all_alpha, ...
    all_error_norm2, ...
    all_inverse_error_norm2, ...
    all_rcond_S, ...
    all_rcond_S_extracted, ...
    'VariableNames', { ...
        'N', ...
        'test_index', ...
        'alpha', ...
        'error_norm2', ...
        'inverse_error_norm2', ...
        'rcond_S', ...
        'rcond_S_extracted'});

data_table = sortrows(data_table, {'N', 'test_index'});

all_N = data_table.N;
all_test_index = data_table.test_index;
all_alpha = data_table.alpha;
all_error_norm2 = data_table.error_norm2;
all_inverse_error_norm2 = data_table.inverse_error_norm2;
all_rcond_S = data_table.rcond_S;
all_rcond_S_extracted = data_table.rcond_S_extracted;

completed_Ns = unique(all_N, 'stable');

% -------------------------------------------------------------------------
% Plot forward spectral error
% -------------------------------------------------------------------------

fig_forward = figure;
hold on;

scatter(all_N, all_error_norm2, ...
    90, ...
    'o', ...
    'LineWidth', 1.2);

box on;
grid on;

ax = gca;

set(ax, 'XScale', 'log');
set(ax, 'YScale', 'linear');

xticks(completed_Ns);
xticklabels(string(completed_Ns));

xlim([completed_Ns(1) / sqrt(2), completed_Ns(end) * sqrt(2)]);

y_data = all_error_norm2(isfinite(all_error_norm2));

if isempty(y_data)

    ylim([0, 1]);

else

    y_min = min(y_data);
    y_max = max(y_data);

    if y_min == y_max
        ylim([0, y_max + 1]);
    else
        margin = 0.12 * (y_max - y_min);
        ylim([max(0, y_min - margin), y_max + margin]);
    end

end

ax.GridAlpha = 0.20;
ax.MinorGridAlpha = 0.05;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

xlabel('$N$', 'Interpreter', 'latex');

ylabel('$\|S - \alpha \widetilde{S}\|_2$', ...
    'Interpreter', 'latex');
% 
% legend({'single configurations'}, ...
%     'Interpreter', 'latex', ...
%     'Location', 'best');

drawnow;

forward_plot_filename = fullfile(results_folder, ...
    'loaded_forward_spectral_error.png');

saveas(fig_forward, forward_plot_filename);

fprintf('Saved forward-error plot to:\n%s\n\n', forward_plot_filename);

% -------------------------------------------------------------------------
% Plot inverse spectral error
% -------------------------------------------------------------------------

fig_inverse = figure;
hold on;

scatter(all_N, all_inverse_error_norm2, ...
    90, ...
    'o', ...
    'LineWidth', 1.2);

box on;
grid on;

ax = gca;

set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');

xticks(completed_Ns);
xticklabels(string(completed_Ns));

xlim([completed_Ns(1) / sqrt(2), completed_Ns(end) * sqrt(2)]);

y_data = all_inverse_error_norm2( ...
    isfinite(all_inverse_error_norm2) & all_inverse_error_norm2 > 0);

if isempty(y_data)

    ylim([1e-16, 1]);

else

    y_min = min(y_data);
    y_max = max(y_data);

    if y_min == y_max
        ylim([y_min / 10, y_max * 10]);
    else
        ylim([y_min / 2, y_max * 2]);
    end

end

ax.GridAlpha = 0.20;
ax.MinorGridAlpha = 0.05;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

xlabel('$N$', 'Interpreter', 'latex');

ylabel('$\|S^{-1} - (\alpha \widetilde{S})^{-1}\|_2$', ...
    'Interpreter', 'latex');
% 
% legend({'single configurations'}, ...
%     'Interpreter', 'latex', ...
%     'Location', 'best');

drawnow;

inverse_plot_filename = fullfile(results_folder, ...
    'loaded_inverse_spectral_error.png');

saveas(fig_inverse, inverse_plot_filename);

fprintf('Saved inverse-error plot to:\n%s\n\n', inverse_plot_filename);

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------

fprintf('\nSummary\n\n');

fprintf([' N        test      alpha          ', ...
         '||S-alpha*S_tilde||_2       ', ...
         '||inv(S)-inv(alpha*S_tilde)||_2\n']);

fprintf('------------------------------------------------------------------------------------------\n');

for q = 1:length(all_N)

    fprintf('%-8d %-8d %.6e   %.6e                 %.6e\n', ...
        all_N(q), ...
        all_test_index(q), ...
        all_alpha(q), ...
        all_error_norm2(q), ...
        all_inverse_error_norm2(q));

end

fprintf('\nDone.\n');