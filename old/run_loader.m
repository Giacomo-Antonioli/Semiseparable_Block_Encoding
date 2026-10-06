% run_loader_example.m
%
% Example script for SemiseparableGeneratorLoader.
%
% This script loads the saved generators from:
%
%   semiseparable_single_generators.txt
%
% and shows how to:
%
%   1. select the results folder using a GUI;
%   2. create the loader;
%   3. print all available configurations;
%   4. retrieve one generator by index;
%   5. retrieve one generator by N and test_index;
%   6. iterate through all generators using a standard for loop;
%   7. iterate through all generators using forEach.

clear;
clc;
close all;

% -------------------------------------------------------------------------
% Select results folder using GUI
% -------------------------------------------------------------------------

results_folder = uigetdir(pwd, 'Select the results folder');

if isequal(results_folder, 0)
    error('No folder selected.');
end

generators_filename = fullfile(results_folder, ...
    'semiseparable_single_generators.txt');

if ~isfile(generators_filename)
    error('Generator file not found in selected folder:\n%s', ...
        generators_filename);
end

% -------------------------------------------------------------------------
% Create loader
% -------------------------------------------------------------------------

loader = SemiseparableGeneratorLoader(generators_filename);

fprintf('Loaded generators from:\n%s\n\n', loader.filename());

fprintf('Number of stored generators: %d\n\n', loader.numGenerators());

% -------------------------------------------------------------------------
% Show available configurations
% -------------------------------------------------------------------------

configs = loader.configurations();

disp('Available configurations:');
disp(configs);

% -------------------------------------------------------------------------
% Retrieve first generator by internal index
% -------------------------------------------------------------------------

fprintf('\nRetrieve first generator by internal index\n');

[x, y, N, test_index] = loader.get(1);

S = tril(x * y') + triu(y * x', 1);

fprintf('index = %d, N = %d, test_index = %d\n', ...
    1, N, test_index);

fprintf('norm(S, 2) = %.16e\n', norm(S, 2));

% -------------------------------------------------------------------------
% Retrieve generator by N and test_index
% -------------------------------------------------------------------------

fprintf('\nRetrieve generator by N and test_index\n');

target_N = configs.N(1);
target_test_index = configs.test_index(1);

[x, y, internal_index] = loader.getByNAndTestIndex( ...
    target_N, target_test_index);

S = tril(x * y') + triu(y * x', 1);

fprintf('N = %d, test_index = %d, internal_index = %d\n', ...
    target_N, target_test_index, internal_index);

fprintf('norm(S, 2) = %.16e\n', norm(S, 2));

% -------------------------------------------------------------------------
% Build matrix directly using class method
% -------------------------------------------------------------------------

fprintf('\nBuild matrix directly using loader.matrix(1)\n');

S_direct = loader.matrix(1);

fprintf('norm(S_direct, 2) = %.16e\n', norm(S_direct, 2));

% -------------------------------------------------------------------------
% Iterate using standard MATLAB for loop
% -------------------------------------------------------------------------

fprintf('\nIterating with standard for loop\n');

for q = 1:loader.numGenerators()

    [x, y, N, test_index] = loader.get(q);

    S = tril(x * y') + triu(y * x', 1);

    fprintf('index = %-3d  N = %-4d  test = %-3d  norm(S,2) = %.6e\n', ...
        q, N, test_index, norm(S, 2));

end

% -------------------------------------------------------------------------
% Iterate using forEach
% -------------------------------------------------------------------------

fprintf('\nIterating with forEach\n');

loader.forEach(@printGeneratorInfo);

% -------------------------------------------------------------------------
% Local helper function
% -------------------------------------------------------------------------

function printGeneratorInfo(x, y, N, test_index)

    S = tril(x * y') + triu(y * x', 1);

    fprintf('N = %-4d  test = %-3d  norm(S,2) = %.6e\n', ...
        N, test_index, norm(S, 2));

end