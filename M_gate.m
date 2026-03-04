function totalCirc = M_gate(n,offset)
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


  % --- Build circuit ---
  tot_wires = 2*n + 2;
Zss_circ=Zs_Gate(n,1);
Zg_circ=Zg_gate(n);
totalCirc=qclab.QCircuit(2*n+2,offset);
totalCirc.push_back(Zg_circ);
totalCirc.push_back(Zss_circ);
totalCirc.push_back(qclab.qgates.Hadamard(0));
% totalCirc.asBlock('M')


  
end
