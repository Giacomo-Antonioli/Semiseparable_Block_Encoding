function circ = make_circuit(psi,offset)
    %MAKE_CIRCUIT Builds a circuit to prepare a quantum state |ψ>
    %
    %   U_circ = make_circuit(psi)
    %
    %   Inputs:
    %       psi - column vector of the quantum state coefficients
    %
    %   Outputs:
    %       U_circ - structure representing the quantum circuit
    %                (you need to implement gates yourself)
    
    if nargin == 1, offset = 0; end
    L = length(psi);
    n = log2(L);
    
    % Get the rotation angles and phases
    [thetas, phis, gammas] = get_all_angles(psi);

    % Build basis states for iterative control
    basis_states = cell(1,n);
    for j1 = 1:n
        basis_states{j1} = make_list(2^(n-j1), n-j1);
    end
    basis_states = flip(basis_states);  % Reverse order
    
    % Initialize a structure to hold the circuit
    U_circ = struct();
    U_circ.n = n;
    U_circ.gates = {};  % List of gates will go here
    
    % First Y rotation on the first qubit
    %%%%%%%%%%%%%TO MODIFY
    % Implement single-qubit Y rotation on qubit n
    circ=qclab.QCircuit(n,offset);
    circ.push_back(qclab.qgates.RotationY(0, thetas{1}(1)))
  
    % Iteratively apply controlled rotations
    for j1 = 2:length(thetas)
        curr_thetas = thetas{j1};
        curr_phis = phis{j1};
        curr_gammas = gammas{j1};
        
        for j2 = 1:length(curr_thetas)
            now_theta = curr_thetas(j2);
            now_phi = curr_phis(j2);
            now_gamma = curr_gammas(j2);
            
            % Apply X gates for zero qubits
            for j3 = 1:length(basis_states{j1}{j2})
                if basis_states{j1}{j2}(j3) == '0'
                    %%%%%%%%%%%%%TO MODIFY
                    % Apply X gate on qubit (n - j3)
                    circ.push_back(qclab.qgates.PauliX(j3-1))
                   % U_circ.gates{end+1} = struct('gate','X','targets',n-j3);
                end
            end
            
            % Multi-controlled rotations
            %%%%%%%%%%%%%TO MODIFY

% 0-based qubit indexing + opposite of Qiskit ordering
target_qubit =  j1-1;

control_qubits =  (0:(j1-1-1));

control_qubits_phase = (0:(j1-2));

t=size(control_qubits);
circ.push_back(qclab.qgates.MCRotationY(control_qubits,target_qubit,ones(1,t(2)),now_theta));
%% ---- Multi-controlled RY ----

circ.push_back(qclab.qgates.MCRotationZ(control_qubits,target_qubit,ones(1,t(2)),now_gamma));

%% ---- Multi-controlled RZ ----


%% ---- Multi-controlled PHASE = MCX + RZ(phi) + MCX ----
% circ.push_back(qclab.qgates.MCX(control_qubits_phase,target_qubit));
% % First MCX
% 
% circ.push_back(qclab.qgates.Phase(target_qubit,now_phi));
% % RZ(phi) on the target
% 
% 
% % Second MCX (uncompute)
% circ.push_back(qclab.qgates.MCX(control_qubits_phase,target_qubit));

circ.push_back(qclab.qgates.MCPhase(control_qubits_phase,target_qubit,ones(1,t(2)),now_phi))


%% ---- Multi-controlled PHASE(-gamma) = MCX + RZ(-gamma) + MCX ----
% circ.push_back(qclab.qgates.MCX(control_qubits_phase,target_qubit));
% 
% % First MCX
% 
% circ.push_back(qclab.qgates.Phase(target_qubit,-now_gamma));
% 
% % RZ(-gamma)
% 
% circ.push_back(qclab.qgates.MCX(control_qubits_phase,target_qubit));
circ.push_back(qclab.qgates.MCPhase(control_qubits_phase,target_qubit,ones(1,t(2)),-now_gamma))
% Second MCX (uncompute)


            
            % Uncompute X gates
            for j3 = 1:length(basis_states{j1}{j2})
                if basis_states{j1}{j2}(j3) == '0'
                    %%%%%%%%%%%%%TO MODIFY
                    % Apply X gate again to revert
                    circ.push_back(qclab.qgates.PauliX(j3-1));


                end
            end
           %  circ.barrier(true)
           % circ.draw
          
        end
    end
end
