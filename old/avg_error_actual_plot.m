% test_semiseparable_circuit_scaling_simulate.m
%
% Randomized scaling test for build_semiseparable_circuit using simulate().
%
% For each N = 2^p, this script runs num_tests random tests.
% It reconstructs the first N columns of the operator by simulating the
% circuit on the first N computational basis states.
%
% Results are saved in a timestamped folder:
%
%   Results_N<first>_to_N<last>_<timestamp>
%
% with one row per N in:
%
%   semiseparable_simulate_errors.txt

clear;
clc;
close all;

rng(42);

% -------------------------------------------------------------------------
% Test configuration
% -------------------------------------------------------------------------

% Powers of two to test
Ns = 2.^(1:3);       % N = 2, 4, 8, 16, 32

% Number of random configurations for each N
num_tests = 50;

% -------------------------------------------------------------------------
% Output folder and file
% -------------------------------------------------------------------------

N_range_string = sprintf('N%d_to_N%d', Ns(1), Ns(end));
timestamp_string = datestr(now, 'yyyymmdd_HHMMSS');

results_folder = sprintf('Results_%s_%s', N_range_string, timestamp_string);

if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

output_filename = fullfile(results_folder, ...
    'semiseparable_simulate_errors.txt');

% Create/reset output file and write header
fid = fopen(output_filename, 'w');

fprintf(fid, ['N num_tests ', ...
              'avg_abs_err_2 avg_rel_err_2 ', ...
              'avg_abs_err_fro avg_rel_err_fro ', ...
              'max_abs_err_2 max_abs_err_fro\n']);

fclose(fid);

fprintf("Saving results in folder:\n%s\n\n", results_folder);

% -------------------------------------------------------------------------
% Storage for plotted quantities
% -------------------------------------------------------------------------

avg_abs_err_2   = zeros(length(Ns), 1);
avg_rel_err_2   = zeros(length(Ns), 1);

avg_abs_err_fro = zeros(length(Ns), 1);
avg_rel_err_fro = zeros(length(Ns), 1);

max_abs_err_2   = zeros(length(Ns), 1);
max_abs_err_fro = zeros(length(Ns), 1);

% -------------------------------------------------------------------------
% Main test loop
% -------------------------------------------------------------------------

for iN = 1:length(Ns)

    N = Ns(iN);
    n = log2(N);

    abs_err_2   = zeros(num_tests, 1);
    rel_err_2   = zeros(num_tests, 1);

    abs_err_fro = zeros(num_tests, 1);
    rel_err_fro = zeros(num_tests, 1);

    fprintf("Testing N = %d\n", N);

    for k = 1:num_tests

        % Random normalized vectors
        x = randn(N, 1);
        x = x / norm(x);

        y = randn(N, 1);
        y = y / norm(y);

        % Target semiseparable matrix
        S = tril(x * y') + triu(y * x', 1);

        % Build circuit
        circ = build_semiseparable_circuit(x, y);

        % Number of qubits in the circuit.
        %
        % Depending on your QCLAB version, you may need one of:
        %
        %   nbQubits = circ.nbQubits;
        %   nbQubits = circ.nbQubits();
        %
        nbQubits = circ.nbQubits;

        dim = 2^nbQubits;

        % Reconstruct the first N columns of the full operator by simulation.
        %
        % Column j is obtained by simulating the circuit on computational
        % basis state |j-1>.
        first_N_cols = zeros(dim, N);

        for j = 1:N

            % MATLAB index j corresponds to computational basis index j-1
            basis_index = j - 1;

            % Binary string input of length nbQubits
            input_state = dec2bin(basis_index, nbQubits);

            % Simulate circuit
            sim = circ.simulate(input_state);

            % Extract final state from QSimulate object
            state = sim.states();

            % Store output state as column j
            first_N_cols(:, j) = state(:);

        end

        % Extract top-left N x N block and rescale block encoding
        myS = full(first_N_cols(1:N, 1:N) * 2^(n/2 + 1));

        % Spectral norm errors
        abs_err_2(k) = norm(myS - S, 2);
        rel_err_2(k) = abs_err_2(k) / norm(S, 2);

        % Frobenius norm errors
        abs_err_fro(k) = norm(myS - S, 'fro');
        rel_err_fro(k) = abs_err_fro(k) / norm(S, 'fro');

    end

    % Average errors for this N
    avg_abs_err_2(iN)   = mean(abs_err_2);
    avg_rel_err_2(iN)   = mean(rel_err_2);

    avg_abs_err_fro(iN) = mean(abs_err_fro);
    avg_rel_err_fro(iN) = mean(rel_err_fro);

    % Maximum errors for diagnostics
    max_abs_err_2(iN)   = max(abs_err_2);
    max_abs_err_fro(iN) = max(abs_err_fro);

    % Save one row immediately after completing this N
    fid = fopen(output_filename, 'a');

    fprintf(fid, '%d %d %.16e %.16e %.16e %.16e %.16e %.16e\n', ...
        N, ...
        num_tests, ...
        avg_abs_err_2(iN), ...
        avg_rel_err_2(iN), ...
        avg_abs_err_fro(iN), ...
        avg_rel_err_fro(iN), ...
        max_abs_err_2(iN), ...
        max_abs_err_fro(iN));

    fclose(fid);

    fprintf("Saved result for N = %d to:\n%s\n\n", ...
        N, output_filename);

end

% -------------------------------------------------------------------------
% Print summary table
% -------------------------------------------------------------------------

fprintf("\nSummary over %d random tests per N\n\n", num_tests);
fprintf(" N        avg abs 2      avg rel 2      avg abs fro    avg rel fro\n");
fprintf("---------------------------------------------------------------------\n");

for iN = 1:length(Ns)
    fprintf("%-8d %.6e   %.6e   %.6e   %.6e\n", ...
        Ns(iN), ...
        avg_abs_err_2(iN), ...
        avg_rel_err_2(iN), ...
        avg_abs_err_fro(iN), ...
        avg_rel_err_fro(iN));
end

fprintf("\nSaved all results to:\n%s\n", output_filename);
% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

% 1. Average absolute error - spectral norm
figure;
loglog(Ns, avg_abs_err_2, '-o', 'LineWidth', 1.5);
grid on;
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\frac{1}{M}\sum_{k=1}^{M} \|\widetilde{S}_k - S_k\|_2$', ...
    'Interpreter', 'latex');
title('Average absolute error in spectral norm', ...
    'Interpreter', 'latex');

saveas(gcf, fullfile(results_folder, ...
    'average_absolute_error_spectral.png'));

% 2. Average relative error - spectral norm
figure;
loglog(Ns, avg_rel_err_2, '-o', 'LineWidth', 1.5);
grid on;
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\frac{1}{M}\sum_{k=1}^{M} \frac{\|\widetilde{S}_k - S_k\|_2}{\|S_k\|_2}$', ...
    'Interpreter', 'latex');
title('Average relative error in spectral norm', ...
    'Interpreter', 'latex');

saveas(gcf, fullfile(results_folder, ...
    'average_relative_error_spectral.png'));

% 3. Average absolute error - Frobenius norm
figure;
loglog(Ns, avg_abs_err_fro, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\frac{1}{M}\sum_{k=1}^{M} \|\widetilde{S}_k - S_k\|_F$', ...
    'Interpreter', 'latex');
title('Average absolute error in Frobenius norm', ...
    'Interpreter', 'latex');

saveas(gcf, fullfile(results_folder, ...
    'average_absolute_error_frobenius.png'));

% 4. Average relative error - Frobenius norm
figure;
loglog(Ns, avg_rel_err_fro, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\frac{1}{M}\sum_{k=1}^{M} \frac{\|\widetilde{S}_k - S_k\|_F}{\|S_k\|_F}$', ...
    'Interpreter', 'latex');
title('Average relative error in Frobenius norm', ...
    'Interpreter', 'latex');

saveas(gcf, fullfile(results_folder, ...
    'average_relative_error_frobenius.png'));

% 5. Combined summary plot
figure;
loglog(Ns, avg_abs_err_2, '-o', 'LineWidth', 1.5);
hold on;
loglog(Ns, avg_rel_err_2, '--o', 'LineWidth', 1.5);
loglog(Ns, avg_abs_err_fro, '-s', 'LineWidth', 1.5);
loglog(Ns, avg_rel_err_fro, '--s', 'LineWidth', 1.5);
grid on;

xlabel('$N$', 'Interpreter', 'latex');
ylabel('Error', 'Interpreter', 'latex');
title(['Average reconstruction errors: ', ...
       '$\widetilde{S}_k = 2^{n/2+1} U_k(1:N,1:N)$'], ...
       'Interpreter', 'latex');

legend( ...
    '$\frac{1}{M}\sum_{k=1}^{M}\|\widetilde{S}_k-S_k\|_2$', ...
    '$\frac{1}{M}\sum_{k=1}^{M}\frac{\|\widetilde{S}_k-S_k\|_2}{\|S_k\|_2}$', ...
    '$\frac{1}{M}\sum_{k=1}^{M}\|\widetilde{S}_k-S_k\|_F$', ...
    '$\frac{1}{M}\sum_{k=1}^{M}\frac{\|\widetilde{S}_k-S_k\|_F}{\|S_k\|_F}$', ...
    'Interpreter', 'latex', ...
    'Location', 'best');

saveas(gcf, fullfile(results_folder, ...
    'average_errors_combined.png'));

fprintf('\nSaved plots in folder:\n%s\n', results_folder);