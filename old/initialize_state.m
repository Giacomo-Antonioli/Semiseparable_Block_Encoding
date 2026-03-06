function circ = initialize_state(n)
    circ=qclab.QCircuit(2*n+2,n);
    %SWAPGATE
    t=qclab.QCircuit(n+2);
    t.push_back(qclab.qgates.SWAP(0,n+1));
    SWAPGATE=t.matrix;
    size(SWAPGATE)
    circ.push_back(qclab.qgates.Hadamard(0))
    for i=0:n-1
       1+i:2*n+i
    circ.push_back(qclab.qgates.MCMatrixGate(0,1+i:n+2+i,SWAPGATE))

    end

end