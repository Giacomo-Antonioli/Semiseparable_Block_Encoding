classdef test_unsymmetric_correctness < matlab.unittest.TestCase
  % Verifies correctness of the unsymmetric semiseparable block-encoding
  % construction S = tril(u*v') + triu(x*y', 1).

  properties (TestParameter)
    N = {4, 8};
  end

  methods (TestClassSetup)
    function addPaths(~)
      addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
      setup_paths();
    end
  end

  methods (Test)
    function blockEncodingBySimulation(test, N)
      rng(42);
      [S, u, v, x, y, nu, nv, nx, ny] = generateUnsymm(N);

      circ = build_unsymmetric_semiseparable_circuit(u/nu, v/nv, x/nx, y/ny, ...
                                                     nu, nv, nx, ny);

      % extract the encoded block by simulating the first N basis states
      dim = 2^circ.nbQubits;
      S_tilde = zeros(N);
      for j = 1:N
        psi = zeros(dim, 1);
        psi(j) = 1;
        sim = circ.simulate(psi);
        out = sim.states;
        S_tilde(:, j) = out(1:N);
      end

      % S / alpha, alpha = sqrt(N) (||u|| ||v|| + ||x|| ||y||)
      alpha = sqrt(N) * (nu*nv + nx*ny);
      block_error = norm(S/alpha - S_tilde, 2);
      test.verifyLessThan(block_error, 1e-10);
    end
  end
end
