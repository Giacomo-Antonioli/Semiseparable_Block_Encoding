n = 8;            % number of matrices
m = 4;            % matrix size (m x m)
prefix = 'a':'h'; % prefixes for each matrix (length must be n)

M = cell(1,n);
for k = 1:n
    name = prefix(k);           % 'a', 'b', ... 'h'
    M{k} = sym(char(name), [m m]); % creates m-by-m array of symbols with that prefix
end



% Form block diagonal
res = blkdiag(M{:});

circ=qclab.QCircuit(6);

% circ.push_back(qclab.qgates.MCX([1,2],3,[1,0]));
% circ.push_back(qclab.qgates.MCX([2,3],1,[0,1]));
circ.push_back(qclab.qgates.SWAP(0,4));
circ.push_back(qclab.qgates.SWAP(1,5));
% circ.push_back(qclab.qgates.SWAP(1,3));
total=qclab.QCircuit(6);
% total.push_back(qclab.qgates.MCX([1,2],0,[0,1]))
% total.push_back(qclab.qgates.MCX([0,1],2,[1,0]))
total.push_back(circ)
subplot(1,2,1)
tmp=total.matrix;
tosee=kron(eye(2),res) *tmp;
spy(tosee)
subplot(1,2,2)


% total.push_back(qclab.qgates.MCX(1,0,[1]))
% total.push_back(qclab.qgates.MCX([0,1],2,[1,0]));
% total.push_back(qclab.qgates.MCX([1,2],0,[0,1]));
% circ=qclab.QCircuit(4);
% circ.push_back(qclab.qgates.MCX([0],2,[0]));
% circ.push_back(qclab.qgates.MCX([1,2],3,[1,0]));
% circ.push_back(qclab.qgates.MCX([2,3],1,[0,1]));


total.draw

tmp=total.matrix;
finalResult = kron(eye(2),res) * tmp; % Multiply the block diagonal matrix with the circuit matrix
spy(finalResult)