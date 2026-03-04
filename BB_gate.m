function circ = BB_gate(n)
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

  % --- Assertions ---
  assert(n > 0, 'n must be greater than zero.');


  % --- Build circuit ---
  tot_wires = 3*n + 2;

  ZZus_circ=Zb_gate(n);
  Mcirc=M_gate(n,n);
  circ = qclab.QCircuit(tot_wires);
  circ.push_back(Mcirc);
  circ.barrier(true);
  circ.push_back(ZZus_circ );


  
end
