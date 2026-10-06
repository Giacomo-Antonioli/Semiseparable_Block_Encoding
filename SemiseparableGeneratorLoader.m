classdef SemiseparableGeneratorLoader < handle
%SEMISEPARABLEGENERATORLOADER Load saved semiseparable generators.
%
% This class loads generator pairs saved in a text file with columns:
%
%   N test_index entry_index x y
%
% The file is expected to contain one row per entry of the vectors x and y.
% For example, for N = 8 and test_index = 2, the file contains rows:
%
%   8 2 1 x_1 y_1
%   8 2 2 x_2 y_2
%   ...
%   8 2 8 x_8 y_8
%
% The class groups the rows by the pair:
%
%   (N, test_index)
%
% and reconstructs the corresponding generator vectors x and y.
%
% -------------------------------------------------------------------------
% Usage examples
% -------------------------------------------------------------------------
%
% Example 1: Load all saved generators
%
%   filename = fullfile( ...
%       'Results_singles_N2_to_N8_20260515_123456', ...
%       'semiseparable_single_generators.txt');
%
%   loader = SemiseparableGeneratorLoader(filename);
%
%   fprintf('Number of generators: %d\n', loader.numGenerators());
%
%
% Example 2: Retrieve a generator by internal index
%
%   [x, y, N, test_index] = loader.get(1);
%
%   S = tril(x * y') + triu(y * x', 1);
%
%
% Example 3: Retrieve a generator by N and test index
%
%   [x, y, internal_index] = loader.getByNAndTestIndex(8, 2);
%
%   S = tril(x * y') + triu(y * x', 1);
%
%
% Example 4: Iterate using a standard MATLAB for loop
%
%   for q = 1:loader.numGenerators()
%       [x, y, N, test_index] = loader.get(q);
%
%       S = tril(x * y') + triu(y * x', 1);
%
%       fprintf('N = %d, test = %d, norm(S) = %.6e\n', ...
%           N, test_index, norm(S, 2));
%   end
%
%
% Example 5: Iterate using forEach
%
%   loader.forEach(@(x, y, N, k) ...
%       fprintf('N = %d, test = %d\n', N, k));
%
%
% Example 6: Use forEach with an external function
%
%   loader.forEach(@processGenerator);
%
%   function processGenerator(x, y, N, k)
%       S = tril(x * y') + triu(y * x', 1);
%
%       fprintf('N = %d, test = %d, norm(S) = %.6e\n', ...
%           N, k, norm(S, 2));
%   end
%
%
% Example 7: Show available configurations
%
%   configs = loader.configurations();
%
%   disp(configs)
%
%
% Example 8: Build the semiseparable matrix directly
%
%   S = loader.matrix(1);
%
%   S = loader.matrixByNAndTestIndex(8, 2);

    properties (Access = private)
        filename_
        generators_
    end

    methods

        function obj = SemiseparableGeneratorLoader(filename)
        %SEMISEPARABLEGENERATORLOADER Construct loader from generator file.
        %
        % Input:
        %
        %   filename : path to semiseparable_single_generators.txt

            if nargin < 1
                error('You must provide the generators filename.');
            end

            if ~isfile(filename)
                error('File not found: %s', filename);
            end

            obj.filename_ = filename;
            obj.generators_ = obj.loadFile(filename);

        end

        function filename = filename(obj)
        %FILENAME Return the source filename.

            filename = obj.filename_;

        end

        function n = numGenerators(obj)
        %NUMGENERATORS Return the number of stored generator pairs.

            n = numel(obj.generators_);

        end

        function [x, y, N, test_index] = get(obj, index)
        %GET Retrieve one generator pair by internal index.
        %
        % Input:
        %
        %   index : integer between 1 and numGenerators()
        %
        % Output:
        %
        %   x          : first generator vector
        %   y          : second generator vector
        %   N          : matrix size
        %   test_index : configuration index

            if ~isscalar(index) || index ~= floor(index)
                error('Index must be an integer scalar.');
            end

            if index < 1 || index > obj.numGenerators()
                error('Index out of range. Valid range is 1 to %d.', ...
                    obj.numGenerators());
            end

            g = obj.generators_(index);

            x = g.x;
            y = g.y;
            N = g.N;
            test_index = g.test_index;

        end

        function [x, y, index] = getByNAndTestIndex(obj, N, test_index)
        %GETBYNANDTESTINDEX Retrieve generator pair by N and test index.
        %
        % Input:
        %
        %   N          : matrix size
        %   test_index : configuration index
        %
        % Output:
        %
        %   x     : first generator vector
        %   y     : second generator vector
        %   index : internal index of the matching generator

            index = [];

            for q = 1:obj.numGenerators()

                same_N = obj.generators_(q).N == N;
                same_test = obj.generators_(q).test_index == test_index;

                if same_N && same_test
                    index = q;
                    break;
                end

            end

            if isempty(index)
                error('No generator found for N = %d and test_index = %d.', ...
                    N, test_index);
            end

            x = obj.generators_(index).x;
            y = obj.generators_(index).y;

        end

        function configs = configurations(obj)
        %CONFIGURATIONS Return table of available configurations.
        %
        % Output:
        %
        %   configs : table with columns index, N, test_index

            num_generators = obj.numGenerators();

            index_col = zeros(num_generators, 1);
            N_col = zeros(num_generators, 1);
            test_index_col = zeros(num_generators, 1);

            for q = 1:num_generators

                index_col(q) = q;
                N_col(q) = obj.generators_(q).N;
                test_index_col(q) = obj.generators_(q).test_index;

            end

            configs = table(index_col, N_col, test_index_col, ...
                'VariableNames', {'index', 'N', 'test_index'});

        end

        function S = matrix(obj, index)
        %MATRIX Build the semiseparable matrix for one generator pair.
        %
        % Input:
        %
        %   index : internal generator index
        %
        % Output:
        %
        %   S = tril(x*y') + triu(y*x',1)

            [x, y] = obj.get(index);

            S = tril(x * y') + triu(y * x', 1);

        end

        function S = matrixByNAndTestIndex(obj, N, test_index)
        %MATRIXBYNANDTESTINDEX Build semiseparable matrix by N and test index.
        %
        % Input:
        %
        %   N          : matrix size
        %   test_index : configuration index
        %
        % Output:
        %
        %   S = tril(x*y') + triu(y*x',1)

            [x, y] = obj.getByNAndTestIndex(N, test_index);

            S = tril(x * y') + triu(y * x', 1);

        end

        function forEach(obj, fun)
        %FOREACH Apply a function to each stored generator.
        %
        % Input:
        %
        %   fun : function handle with signature
        %
        %       fun(x, y, N, test_index)
        %
        % Example:
        %
        %   loader.forEach(@(x, y, N, k) ...
        %       fprintf('N = %d, test = %d\n', N, k));

            if ~isa(fun, 'function_handle')
                error('Input must be a function handle.');
            end

            for q = 1:obj.numGenerators()

                [x, y, N, test_index] = obj.get(q);

                fun(x, y, N, test_index);

            end

        end

    end

    methods (Access = private)

        function generators = loadFile(~, filename)
        %LOADFILE Load and group all generator entries from file.

            T = readtable(filename, ...
                'FileType', 'text', ...
                'Delimiter', ' ', ...
                'MultipleDelimsAsOne', true);

            required_columns = {'N', 'test_index', 'entry_index', 'x', 'y'};

            for c = 1:numel(required_columns)

                if ~ismember(required_columns{c}, T.Properties.VariableNames)
                    error('Missing required column: %s', required_columns{c});
                end

            end

            pairs = unique([T.N, T.test_index], 'rows');

            generators = struct( ...
                'N', {}, ...
                'test_index', {}, ...
                'x', {}, ...
                'y', {});

            for q = 1:size(pairs, 1)

                N = pairs(q, 1);
                test_index = pairs(q, 2);

                idx = (T.N == N) & (T.test_index == test_index);

                Tq = T(idx, :);
                Tq = sortrows(Tq, 'entry_index');

                if height(Tq) ~= N
                    error(['Incomplete generator for N = %d, ', ...
                           'test_index = %d. Expected %d entries, found %d.'], ...
                           N, test_index, N, height(Tq));
                end

                expected_entries = (1:N)';

                if any(Tq.entry_index ~= expected_entries)
                    error(['Invalid entry indices for N = %d, ', ...
                           'test_index = %d. Expected entry_index = 1:N.'], ...
                           N, test_index);
                end

                generators(q).N = N;
                generators(q).test_index = test_index;
                generators(q).x = Tq.x;
                generators(q).y = Tq.y;

            end

        end

    end

end