function circ = Zus_gate(n, offset)
% ZBI_GATE   Build a Z-bi gate circuit of size n with an optional offset.
%
%   circ = Zbi_gate(n)
%   circ = Zbi_gate(n, offset)
%
%   Inputs:
%     n       - number of qubits (must be > 0)
%     offset  - optional wire offset (default = 0, must satisfy offset < n)
%
%   Output:
%     circ    - constructed qclab.QCircuit object

  % --- Handle default offset ---
  if nargin < 2
      offset = 0;
  end

  % --- Assertions ---
  assert(n > 0, 'n must be greater than zero.');
  assert(offset < n, 'offset must be smaller than n.');

  % --- Build circuit ---
  tot_wires = 2*n + 2;

  circ = qclab.QCircuit(tot_wires, n+1);

  circ.push_back(rightshift(tot_wires-n-1) );

  circ.asBlock("Zus")
end
