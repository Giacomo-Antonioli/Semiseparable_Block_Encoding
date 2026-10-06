
function circ = state_prep(psi,offset)
% STATE_PREP_UNITARY  Builds a unitary whose first column is psi.
%
% Any unitary U satisfying U|0> = |psi> is a valid state preparation.
% We complete the basis via QR decomposition on [psi | I_{dim, dim-1}].
%
% Input:
%   psi - normalized complex column vector of length dim (dim = 2^n)
%
% Output:
%   U   - (dim x dim) unitary matrix with U(:,1) = psi
    
    circ=qclab.QCircuit(log2(size(psi,1)),offset);
    dim = numel(psi);
    psi = psi(:) / norm(psi);

    basis    = [psi, eye(dim, dim - 1)];
    [U, ~]   = qr(basis);

    if norm(U(:,1) - psi) > norm(-U(:,1) - psi)
        U(:,1) = -U(:,1);
    end

    circ.push_back(qclab.qgates.MatrixGate(0:circ.nbQubits-1,U));

end