function circ = InitializeStatevector(psi, offset)
    %INITIALIZESTATEVECTOR Builds a quantum circuit to prepare a statevector |ψ>
    %
    %   circ = InitializeStatevector(psi)
    %   circ = InitializeStatevector(psi, offset)
    %
    %   The quantum state is initialized using the iterative method of controlled 
    %   rotations. The algorithm systematically applies single and multi-controlled 
    %   Y, Z, and phase rotations, conditioned on basis states, to construct the 
    %   target superposition from the computational basis.
    %
    %   Inputs:
    %       psi    - Column vector of quantum state coefficients (length must be 2^n)
    %       offset - Qubit offset for circuit (default: 0)
    %
    %   Outputs:
    %       circ - qclab.QCircuit object representing the state preparation circuit
    
    if nargin == 1, offset = 0; end
    L = length(psi);
    n = log2(L);
    
    [thetas, phis, gammas] = get_all_angles(psi);
    
    basis_states = cell(1, n);
    for j1 = 1:n
        basis_states{j1} = make_list(2^(n-j1), n-j1);
    end
    basis_states = flip(basis_states);
    
    U_circ = struct();
    U_circ.n = n;
    U_circ.gates = {};
    
    circ = qclab.QCircuit(n, offset);
    circ.push_back(qclab.qgates.RotationY(0, thetas{1}(1)))
    
    for j1 = 2:length(thetas)
        curr_thetas = thetas{j1};
        curr_phis = phis{j1};
        curr_gammas = gammas{j1};
        
        for j2 = 1:length(curr_thetas)
            now_theta = curr_thetas(j2);
            now_phi = curr_phis(j2);
            now_gamma = curr_gammas(j2);
            
            for j3 = 1:length(basis_states{j1}{j2})
                if basis_states{j1}{j2}(j3) == '0'
                    circ.push_back(qclab.qgates.PauliX(j3-1))
                end
            end
            
            target_qubit = j1 - 1;
            control_qubits = (0:(j1-1-1));
            control_qubits_phase = (0:(j1-2));
            t = size(control_qubits);
            
            circ.push_back(qclab.qgates.MCRotationY(control_qubits, target_qubit, ones(1, t(2)), now_theta));
            circ.push_back(qclab.qgates.MCRotationZ(control_qubits, target_qubit, ones(1, t(2)), now_gamma));
            circ.push_back(qclab.qgates.MCPhase(control_qubits_phase, target_qubit, ones(1, t(2)), now_phi))
            circ.push_back(qclab.qgates.MCPhase(control_qubits_phase, target_qubit, ones(1, t(2)), -now_gamma))
            
            for j3 = 1:length(basis_states{j1}{j2})
                if basis_states{j1}{j2}(j3) == '0'
                    circ.push_back(qclab.qgates.PauliX(j3-1));
                end
            end
        end
    end
end