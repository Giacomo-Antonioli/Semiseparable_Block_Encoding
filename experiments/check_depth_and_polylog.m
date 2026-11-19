% CHECK_DEPTH_AND_POLYLOG  Independent check of two numbers in Sec. 8.3.
%
%   Part A. Raw (pre-transpilation) depth and size of the exponential-loader
%           block encoding, built in QCLAB, against the closed forms
%               depth n^2+2n+6,  size n^2+8n+3   (current Fix_gate)
%
%   Part B. Polylog fit  C*(log2 N)^p  of the TRANSPILED exponential-loader
%           depth and size. Transpilation is done by Qiskit, so the data are
%           read from results/fable_python/fable_comparison.csv (written by
%           python/semiseparable_vs_fable_qiskit.py).
%           The fit is a least-squares line of log(y) against log(log2 N)
%           over N = 4 ... 1024; R^2 is on the log scale.
%
%   Run from MATLAB:  >> check_depth_and_polylog

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
root = setup_paths();

%% ---- Part A: raw depth and size of the exponential loader ------------
fprintf('Part A: raw depth/size of build_semiseparable_circuit_exp (QCLAB)\n');
fprintf('%3s %6s | %8s %10s | %8s %10s\n', 'n', 'N', ...
        'depth', 'n^2+2n+6', 'size', 'n^2+8n+3');
nA = 2:8;
dep = zeros(size(nA)); siz = zeros(size(nA));
for k = 1:numel(nA)
    n = nA(k); N = 2^n;
    circ = build_semiseparable_circuit_exp(ones(N,1), ones(N,1));
    dep(k) = circ.getDepth();
    siz(k) = circ.getGateCount();
    fprintf('%3d %6d | %8d %10d | %8d %10d\n', n, N, dep(k), ...
            n^2+2*n+6, siz(k), n^2+8*n+3);
end
yn = {'NO', 'YES'};
fprintf('depth = n^2+2n+6 holds for every n = %d..%d: %s\n', nA(1), nA(end), ...
        yn{isequal(dep, nA.^2+2*nA+6) + 1});
fprintf('size  = n^2+8n+3 holds for every n = %d..%d: %s\n\n', nA(1), nA(end), ...
        yn{isequal(siz, nA.^2+8*nA+3) + 1});

%% ---- Part B: polylog exponent p of the transpiled exponential loader --
N = 2.^(2:10);
L = log2(N);

% Transpiled columns of the CSV
T = readtable(fullfile(root, 'results', 'fable_python', ...
                       'fable_comparison.csv'));
T = T(strcmp(T.method, 'exponential'), :);
T = sortrows(T, 'N');
assert(isequal(T.N(:)', N), 'CSV does not hold N = 4 ... 1024');
new_depth = T.depth_transpiled(:)';
new_size  = T.size_transpiled(:)';

fprintf('Part B: fit y = C*(log2 N)^p, N = 4 ... 1024\n');
[p_d, R2_d] = report('depth ', L, new_depth);
[p_s, R2_s] = report('size  ', L, new_size);

fprintf('\nNumbers for the paper (Sec. 8.3, polylogarithmic fit):\n');
fprintf('  p (depth) = %.2f   R^2 (depth) = %.4f\n', p_d, R2_d);
fprintf('  p (size)  = %.2f   R^2 (size)  = %.4f\n', p_s, R2_s);
fprintf('  LaTeX: with $p\\approx %.1f$ for the depth and $p\\approx %.2f$ for the size ($R^2>%.3f$)\n', ...
        p_d, p_s, floor(min(R2_d, R2_s)*1000)/1000);
fprintf('  (C and the local exponents above are diagnostics; they are not quoted.)\n');

%% ---- local function ----------------------------------------------------
function [p, R2] = report(name, L, y)
    x = log(L); z = log(y);
    c = polyfit(x, z, 1);
    p = c(1);
    zhat = polyval(c, x);
    R2 = 1 - sum((z - zhat).^2) / sum((z - mean(z)).^2);
    loc = diff(log(y)) ./ diff(log(L));       % local exponents
    fprintf('%s p = %.3f  C = %.3g  R^2 = %.4f\n', name, p, exp(c(2)), R2);
    fprintf('    local exponents: %s\n', sprintf('%.2f ', loc));
end
