function circ = Zs_Gate(n, offset)
% ZS_GATE   Build a Zs gate circuit of size n with an optional offset.
%
%   circ = Zs_Gate(n)
%   circ = Zs_Gate(n, offset)
%
%   This function constructs a qclab.QCircuit implementing a sequence of
%   Z-shift (Zs) operations arranged over 2*n + 1 wires. Each iteration adds
%   a left-shift gate acting on decreasing target qubits, encoding the Zs
%   pattern.
%
%   Inputs:
%     n        - Number of qubits defining the Zs structure.
%                Must be a positive integer (n > 0).
%
%     offset   - (Optional) Wire offset passed to qclab.QCircuit.
%                Default = 0.
%                Must satisfy: offset < n.
%
%   Output:
%     circ     - Constructed qclab.QCircuit object representing the Zs gate.
%
%   The total number of wires in the circuit is 2*n + 1.
%

  % --- Default offset ---
  if nargin < 2
      offset = 0;
  end

  % --- Assertions ---
  assert(n > 0, 'n must be greater than zero.');


  % --- Build circuit ---
  tot_wires = 2*n + 1;
  circ = qclab.QCircuit(tot_wires, offset);

  for i = 0:n-1
      circ.push_back( ...
          leftshift( ...
              tot_wires, ...
              tot_wires - n : tot_wires - 1 - i, ...   % control wires
              [n-1-i, n], ...                           % target wires
              [1, 0] ) );                               % shift values
  end
    circ.asBlock("Zs")
end
