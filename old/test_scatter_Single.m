% test_semiseparable_single_error_plot_simulate.m
%
% Test all configurations for N = 2,4,8,16,32,64,128.
% For each N, run 4 random configurations.
%
% For each configuration, reconstruct the raw encoded block S_tilde and
% compute
%
%   || S - alpha S_tilde ||_2,
%
% where alpha depends on N:
%
%   alpha = 2^(log2(N)/2 + 1).
%
% The script updates and saves the scatter plot after each completed N.
% It also saves the error data and the generators x,y to text files.
%
% Results are saved in:
%
%   Results_singles_N2_to_N128_<timestamp>

clear;
clc;
close all;

rng(42);

% -------------------------------------------------------------------------
% Test configuration
% -------------------------------------------------------------------------

Ns = 2.^(1:3);      % N = 2, 4, 8, 16, 32, 64, 128
num_tests = 4;

% -------------------------------------------------------------------------
% Output folder and files
% -------------------------------------------------------------------------

N_range_string = sprintf('N%d_to_N%d', Ns(1), Ns(end));
timestamp_string = datestr(now, 'yyyymmdd_HHMMSS');

results_folder = sprintf('Results_singles_%s_%s', ...
    N_range_string, timestamp_string);

if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

output_filename = fullfile(results_folder, ...
    'semiseparable_single_spectral_errors.txt');

generators_filename = fullfile(results_folder, ...
    'semiseparable_single_generators.txt');

fid = fopen(output_filename, 'w');
fprintf(fid, 'N test_index alpha error_norm2\n');
fclose(fid);

fid = fopen(generators_filename, 'w');
fprintf(fid, 'N test_index entry_index x y\n');
fclose(fid);

fprintf("Saving results in folder:\n%s\n\n", results_folder);

% -------------------------------------------------------------------------
% Storage for plot
% -------------------------------------------------------------------------

all_N = [];
all_test_index = [];
all_alpha = [];
all_error_norm2 = [];

% -------------------------------------------------------------------------
% Create figure once
% -------------------------------------------------------------------------

fig = figure;

% -------------------------------------------------------------------------
% Main test loop
% -------------------------------------------------------------------------

for iN = 1:length(Ns)

    N = Ns(iN);
    n = log2(N);

    % N-dependent scaling factor
    alpha = 2^(n/2 + 1);

    fprintf("Testing N = %d, alpha = %.6e\n", N, alpha);

    for k = 1:num_tests

        fprintf("  Configuration %d / %d\n", k, num_tests);

        % -----------------------------------------------------------------
        % Random normalized vectors
        % -----------------------------------------------------------------

        x = randn(N, 1);
        x = x / norm(x);

        y = randn(N, 1);
        y = y / norm(y);

        % -----------------------------------------------------------------
        % Save generators immediately
        % -----------------------------------------------------------------

        fid_gen = fopen(generators_filename, 'a');

        for ell = 1:N
            fprintf(fid_gen, '%d %d %d %.16e %.16e\n', ...
                N, ...
                k, ...
                ell, ...
                x(ell), ...
                y(ell));
        end

        fclose(fid_gen);

        % -----------------------------------------------------------------
        % Target semiseparable matrix
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

            % Clear previous progress message from the command window
            fprintf(repmat('\b', 1, length(progress_msg)));

            progress_msg = sprintf('    Simulating column %d / %d', j, N);
            fprintf('%s', progress_msg);

            basis_index = j - 1;

            input_state = dec2bin(basis_index, nbQubits);

            sim = circ.simulate(input_state);

            state = sim.states();

            first_N_cols(:, j) = state(:);

        end

        % Move to next line after finishing the column loop
        fprintf('\n');

        % -----------------------------------------------------------------
        % Extract raw encoded top-left N x N block
        % -----------------------------------------------------------------

        S_tilde = full(first_N_cols(1:N, 1:N));

        % -----------------------------------------------------------------
        % Compute spectral-norm error
        % -----------------------------------------------------------------

        error_norm2 = norm(S - alpha * S_tilde, 2);

        % -----------------------------------------------------------------
        % Store data
        % -----------------------------------------------------------------

        all_N(end + 1, 1) = N;
        all_test_index(end + 1, 1) = k;
        all_alpha(end + 1, 1) = alpha;
        all_error_norm2(end + 1, 1) = error_norm2;

        % -----------------------------------------------------------------
        % Save error data immediately after each configuration
        % -----------------------------------------------------------------

        fid = fopen(output_filename, 'a');

        fprintf(fid, '%d %d %.16e %.16e\n', ...
            N, ...
            k, ...
            alpha, ...
            error_norm2);

        fclose(fid);

    end

    % ---------------------------------------------------------------------
    % Update and save plot after each completed N
    % ---------------------------------------------------------------------

    figure(fig);
    clf(fig);
    hold on;


% ---------------------------------------------------------------------
% Scatter data only, no jitter, no mean
% ---------------------------------------------------------------------

scatter(all_N, all_error_norm2, ...
    90, ...
    'r*', ...
    'LineWidth', 1.2);

    % ---------------------------------------------------------------------
    % Axes and labels
    % ---------------------------------------------------------------------

    box on;
    grid on;

    ax = gca;

    set(ax, 'XScale', 'log');
    set(ax, 'YScale', 'linear');

    completed_Ns = Ns(1:iN);

    xticks(completed_Ns);
    xticklabels(string(completed_Ns));

    xlim([Ns(1) / sqrt(2), Ns(iN) * sqrt(2)]);

    y_min = min(all_error_norm2);
    y_max = max(all_error_norm2);

    if y_min == y_max
        ylim([0, y_max + 1]);
    else
        margin = 0.12 * (y_max - y_min);
        ylim([max(0, y_min - margin), y_max + margin]);
    end

    % Lighter and less cluttered grid
    ax.GridAlpha = 0.20;
    ax.MinorGridAlpha = 0.05;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

    xlabel('$N$', 'Interpreter', 'latex');

    ylabel('$\|S - \alpha \widetilde{S}\|_2$', ...
        'Interpreter', 'latex');

    legend({'single configurations'}, ...
        'Interpreter', 'latex', ...
        'Location', 'best');

    % No graph title.

    drawnow;

    current_plot_filename = fullfile(results_folder, ...
        sprintf('single_spectral_error_up_to_N%d.png', N));

    saveas(fig, current_plot_filename);

    final_plot_filename = fullfile(results_folder, ...
        'single_spectral_error_current.png');

    saveas(fig, final_plot_filename);

    fprintf("Saved updated plot up to N = %d:\n%s\n\n", ...
        N, current_plot_filename);

end

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------

fprintf("\nSummary\n\n");
fprintf(" N        test      alpha          ||S - alpha S_tilde||_2\n");
fprintf("----------------------------------------------------------------\n");

for idx = 1:length(all_N)

    fprintf("%-8d %-8d %.6e   %.6e\n", ...
        all_N(idx), ...
        all_test_index(idx), ...
        all_alpha(idx), ...
        all_error_norm2(idx));

end

fprintf("\nSaved all error data to:\n%s\n", output_filename);
fprintf("Saved all generators to:\n%s\n", generators_filename);
fprintf("Saved final plot in folder:\n%s\n", results_folder);