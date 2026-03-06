function circ = Zb_gate(n, offset)
% ZB_GATE   Build a Zb gate circuit of size n with an optional offset.
%
%   circ = Zb_gate(n)
%   circ = Zb_gate(n, offset)
%
%   This function constructs a qclab.QCircuit implementing a sequence of
%   Z-bottom (Zb) shift operations arranged over 2*n + 1 wires. Each
%   iteration appends a right-shift gate acting on increasing wire ranges,
%   producing the Zb pattern.
%
%   Inputs:
%     n        - Number of qubits defining the Zb structure.
%                Must be a positive integer (n > 0).
%
%     offset   - (Optional) Wire offset passed to qclab.QCircuit.
%                Default = 0.
%                Must satisfy: offset < n.
%
%   Output:
%     circ     - Constructed qclab.QCircuit object representing the Zb gate.
%
%   The total number of wires in the circuit is 2*n + 1.
%

  % --- Default offset ---
  if nargin < 2
      offset = 0;
  end

  % --- Assertions ---
  assert(n > 0, 'n must be greater than zero.');
  %assert(offset < n, 'offset must be smaller than n.');

  % --- Build circuit ---
  tot_wires = 2*n + 1;
  circ = qclab.QCircuit(tot_wires, offset);

  for i = 0:n-1
      circ.push_back( ...
          rightshift( ...
              tot_wires, ...
              tot_wires - n - 1 : tot_wires - n + i, ...   % active wires
              i ) );                                        % shift amount
  end

  circ.asBlock('Zb')

end
