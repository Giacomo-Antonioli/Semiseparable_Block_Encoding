


syms a [4 4] matrix
syms b [4 4] matrix
syms c [4 4] matrix
syms d [4 4] matrix

A = symmatrix2sym(a); % convert to array of symbolic scalars if you need elementwise sym objects
B = symmatrix2sym(b); % convert to array of symbolic scalars if you need elementwise sym objects
C = symmatrix2sym(c); % convert to array of symbolic scalars if you need elementwise sym objects
D = symmatrix2sym(d); % convert to array of symbolic scalars if you need elementwise sym objects

res=blkdiag(A,B,C,D);

circ=qclab.QCircuit(4);

circ.push_back(qclab.qgates.MCX([1,2],3,[1,0]));
circ.push_back(qclab.qgates.MCX([2,3],1,[0,1]));
circ.push_back(qclab.qgates.SWAP(0,2));
circ.push_back(qclab.qgates.SWAP(1,3));

% circ=qclab.QCircuit(4);
% circ.push_back(qclab.qgates.MCX([0],2,[0]));
% circ.push_back(qclab.qgates.MCX([1,2],3,[1,0]));
% circ.push_back(qclab.qgates.MCX([2,3],1,[0,1]));


circ.draw

tmp=circ.matrix;
finalResult = res * tmp % Multiply the block diagonal matrix with the circuit matrix