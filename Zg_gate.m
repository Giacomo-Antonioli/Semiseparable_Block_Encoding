function circ = Zg_gate(n, offset)
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

  circ = qclab.QCircuit(tot_wires, offset);
    tmp = qclab.QCircuit(n+1,0);
  tmp.push_back(rightshift(n+1));

  circ.push_back(tmp);
for i=0:n


  circ.push_back( ...
      leftshift(tot_wires, tot_wires-n-1:tot_wires-n+i-1, i, 0) );
end
 
    tmp = qclab.QCircuit(n+1,0);
  tmp.push_back(leftshift(n+1));
circ.push_back(tmp);
  
end
