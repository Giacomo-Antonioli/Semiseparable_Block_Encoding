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
% and also
%
%   || S^{-1} - (alpha S_tilde)^{-1} ||_2,
%
% and also the off-tridiagonal inverse comparison
%
%   || off(S^{-1}) - off((alpha S_tilde)^{-1}) ||_F,
%
% where
%
%   off(A) = A - tridiag(A),
%
% and alpha depends on N:
%
%   alpha = 2^(log2(N)/2 + 1).
%
% The script updates and saves three scatter plots after each completed N:
%
%   1. forward spectral error;
%   2. inverse spectral error;
%   3. off-tridiagonal inverse comparison.
%
% It also saves the error data to a text file and saves the real and
% extracted matrices to .mat files.
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

Ns = 2.^(2:6);      % N = 2, 4, 8, 16
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

matrices_folder = fullfile(results_folder, 'saved_matrices');

if ~exist(matrices_folder, 'dir')
    mkdir(matrices_folder);
end

output_filename = fullfile(results_folder, ...
    'semiseparable_single_spectral_errors.txt');

fid = fopen(output_filename, 'w');

fprintf(fid, ['N test_index alpha error_norm2 ', ...
              'inverse_error_norm2 ', ...
              'off_tridiagonal_inverse_error_fro ', ...
              'rcond_S rcond_S_extracted\n']);

fclose(fid);

fprintf("Saving results in folder:\n%s\n\n", results_folder);

% -------------------------------------------------------------------------
% Storage for plots
% -------------------------------------------------------------------------

all_N = [];
all_test_index = [];
all_alpha = [];
all_error_norm2 = [];
all_inverse_error_norm2 = [];
all_off_tridiagonal_inverse_error_fro = [];
all_rcond_S = [];
all_rcond_S_extracted = [];

% -------------------------------------------------------------------------
% Create figures once
% -------------------------------------------------------------------------

fig_error = figure;
fig_inverse_error = figure;
fig_off_tridiagonal_error = figure;

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

        % Scaled extracted matrix
        S_extracted = alpha * S_tilde;

        % -----------------------------------------------------------------
        % Save real and extracted matrices
        % -----------------------------------------------------------------

        matrix_filename = fullfile(matrices_folder, ...
            sprintf('matrices_N%d_test%d.mat', N, k));

        save(matrix_filename, ...
            'S', ...
            'S_tilde', ...
            'S_extracted', ...
            'alpha', ...
            'N', ...
            'k');

        % -----------------------------------------------------------------
        % Compute spectral-norm error
        % -----------------------------------------------------------------

        error_norm2 = norm(S - S_extracted, 2);

        % -----------------------------------------------------------------
        % Compute inverse spectral-norm error and off-tridiagonal comparison
        % -----------------------------------------------------------------
        %
        % Instead of explicitly using inv(S), we compute inverse matrices
        % through linear solves:
        %
        %   S_inv = S \ eye(N)
        %
        % This is numerically preferable to inv(S).

        rcond_S = rcond(S);
        rcond_S_extracted = rcond(S_extracted);

        if rcond_S == 0 || rcond_S_extracted == 0

            inverse_error_norm2 = NaN;
            off_tridiagonal_inverse_error_fro = NaN;

            warning(['Singular matrix detected for N = %d, test = %d. ', ...
                     'Inverse errors set to NaN.'], ...
                     N, k);

        else

            S_inv = S \ eye(N);
            S_extracted_inv = S_extracted \ eye(N);

            inverse_error_norm2 = norm(S_inv - S_extracted_inv, 2);

            off_S_inv = S_inv - tridiag_part(S_inv);
            off_S_extracted_inv = S_extracted_inv - ...
                tridiag_part(S_extracted_inv);

            off_tridiagonal_inverse_error_fro = norm( ...
                off_S_inv - off_S_extracted_inv, 'fro');

        end

        % -----------------------------------------------------------------
        % Store data
        % -----------------------------------------------------------------

        all_N(end + 1, 1) = N;
        all_test_index(end + 1, 1) = k;
        all_alpha(end + 1, 1) = alpha;
        all_error_norm2(end + 1, 1) = error_norm2;
        all_inverse_error_norm2(end + 1, 1) = inverse_error_norm2;
        all_off_tridiagonal_inverse_error_fro(end + 1, 1) = ...
            off_tridiagonal_inverse_error_fro;
        all_rcond_S(end + 1, 1) = rcond_S;
        all_rcond_S_extracted(end + 1, 1) = rcond_S_extracted;

        % -----------------------------------------------------------------
        % Save error data immediately after each configuration
        % -----------------------------------------------------------------

        fid = fopen(output_filename, 'a');

        fprintf(fid, '%d %d %.16e %.16e %.16e %.16e %.16e %.16e\n', ...
            N, ...
            k, ...
            alpha, ...
            error_norm2, ...
            inverse_error_norm2, ...
            off_tridiagonal_inverse_error_fro, ...
            rcond_S, ...
            rcond_S_extracted);

        fclose(fid);

    end

    % ---------------------------------------------------------------------
    % Update and save forward-error plot after each completed N
    % ---------------------------------------------------------------------

    figure(fig_error);
    clf(fig_error);
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

    completed_Ns = Ns(1:iN);

    xticks(completed_Ns);
    xticklabels(string(completed_Ns));

    xlim([Ns(1) / sqrt(2), Ns(iN) * sqrt(2)]);

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

    drawnow;

    current_plot_filename = fullfile(results_folder, ...
        sprintf('single_spectral_error_up_to_N%d.png', N));

    saveas(fig_error, current_plot_filename);

    final_plot_filename = fullfile(results_folder, ...
        'single_spectral_error_current.png');

    saveas(fig_error, final_plot_filename);

    fprintf("Saved updated forward-error plot up to N = %d:\n%s\n\n", ...
        N, current_plot_filename);

    % ---------------------------------------------------------------------
    % Update and save inverse-error plot after each completed N
    % ---------------------------------------------------------------------

    figure(fig_inverse_error);
    clf(fig_inverse_error);
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

    xlim([Ns(1) / sqrt(2), Ns(iN) * sqrt(2)]);

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

    drawnow;

    current_inverse_plot_filename = fullfile(results_folder, ...
        sprintf('single_inverse_spectral_error_up_to_N%d.png', N));

    saveas(fig_inverse_error, current_inverse_plot_filename);

    final_inverse_plot_filename = fullfile(results_folder, ...
        'single_inverse_spectral_error_current.png');

    saveas(fig_inverse_error, final_inverse_plot_filename);

    fprintf("Saved updated inverse-error plot up to N = %d:\n%s\n\n", ...
        N, current_inverse_plot_filename);

    % ---------------------------------------------------------------------
    % Update and save off-tridiagonal inverse comparison plot
    % ---------------------------------------------------------------------

    figure(fig_off_tridiagonal_error);
    clf(fig_off_tridiagonal_error);
    hold on;

    scatter(all_N, all_off_tridiagonal_inverse_error_fro, ...
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

    xlim([Ns(1) / sqrt(2), Ns(iN) * sqrt(2)]);

    y_data = all_off_tridiagonal_inverse_error_fro( ...
        isfinite(all_off_tridiagonal_inverse_error_fro) & ...
        all_off_tridiagonal_inverse_error_fro > 0);

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

    ylabel(['$\|\mathrm{off}(S^{-1}) - ', ...
            '\mathrm{off}((\alpha\widetilde{S})^{-1})\|_F$'], ...
        'Interpreter', 'latex');

    drawnow;

    current_off_plot_filename = fullfile(results_folder, ...
        sprintf('single_off_tridiagonal_inverse_error_up_to_N%d.png', N));

    saveas(fig_off_tridiagonal_error, current_off_plot_filename);

    final_off_plot_filename = fullfile(results_folder, ...
        'single_off_tridiagonal_inverse_error_current.png');

    saveas(fig_off_tridiagonal_error, final_off_plot_filename);

    fprintf("Saved updated off-tridiagonal inverse-error plot up to N = %d:\n%s\n\n", ...
        N, current_off_plot_filename);

end

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------

fprintf("\nSummary\n\n");

fprintf([' N        test      alpha          ', ...
         '||S - alpha S_tilde||_2       ', ...
         '||inv(S) - inv(alpha S_tilde)||_2       ', ...
         'off-tridiagonal inverse error\n']);

fprintf("---------------------------------------------------------------------------------------------------------------\n");

for idx = 1:length(all_N)

    fprintf("%-8d %-8d %.6e   %.6e                 %.6e                 %.6e\n", ...
        all_N(idx), ...
        all_test_index(idx), ...
        all_alpha(idx), ...
        all_error_norm2(idx), ...
        all_inverse_error_norm2(idx), ...
        all_off_tridiagonal_inverse_error_fro(idx));

end

fprintf("\nSaved all error data to:\n%s\n", output_filename);
fprintf("Saved all matrices to:\n%s\n", matrices_folder);
fprintf("Saved final forward-error plot in folder:\n%s\n", results_folder);
fprintf("Saved final inverse-error plot in folder:\n%s\n", results_folder);
fprintf("Saved final off-tridiagonal inverse-error plot in folder:\n%s\n", results_folder);

% -------------------------------------------------------------------------
% Helper function
% -------------------------------------------------------------------------

