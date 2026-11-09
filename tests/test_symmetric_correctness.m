classdef test_symmetric_correctness < matlab.unittest.TestCase
  % Verifies correctness of the symmetric semiseparable block-encoding
  % construction S = tril(u*v') + triu(v*u', 1).

  properties (TestParameter)
    N = {2, 4, 8};
  end

  methods (TestClassSetup)
    function addPaths(~)
      addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
      setup_paths();
    end
  end

  methods (Test)
    function unitarityAndBlockEncoding(test, N)
      rng(42);
      u = randn(N, 1);  u = u / norm(u);
      v = randn(N, 1);  v = v / norm(v);
      S = tril(u*v') + triu(v*u', 1);

      circ = build_semiseparable_circuit(u, v);
      U = circ.matrix;

      % (A) unitarity
      unitarity_error = norm(U' * U - eye(size(U)), 2);

      % (B) the top-left N x N block encodes S / (2 sqrt(N))
      S_tilde = 2 * sqrt(N) * U(1:N, 1:N);
      block_error = norm(S_tilde - S, 2);

      test.verifyLessThan(unitarity_error, 1e-10);
      test.verifyLessThan(block_error, 1e-10);
    end
  end
end
