from __future__ import annotations
from qiskit import transpile, QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit.circuit.library import MCXGate
import matplotlib.pyplot as plt
from collections import defaultdict
from fable import fable
from qiskit.quantum_info import Operator
import numpy as np
from qiskit.quantum_info import Statevector
import matplotlib
matplotlib.use('TkAgg')



# --- Verification Flag ---
RUN_VERIFICATION = False

# Plot the individual gate counts (True) or the total circuit size (False)
PLOT_GATE_COUNTS = False
# -------------------------



def qmap(i: int, nq: int) -> int:
    return nq - 1 - i


def mcx_with_states(qc: QuantumCircuit, ctrl_qubits: list[int],
                    target: int, ctrl_states: list[int]) -> None:
    if not ctrl_qubits:
        qc.x(target)
        return

    n_ctrl = len(ctrl_qubits)
    states = list(ctrl_states)
    if len(states) < n_ctrl:
        states += [1] * (n_ctrl - len(states))
    elif len(states) > n_ctrl:
        states = states[:n_ctrl]

    ctrl_state = "".join(str(int(s)) for s in reversed(states))

    gate = MCXGate(n_ctrl, ctrl_state=ctrl_state)
    qc.append(gate, ctrl_qubits + [target])


def _validate_shift_args(qc: QuantumCircuit, targets: list[int],
                         controls: list[int], control_states: list[int]) -> None:
    n = qc.num_qubits
    assert all(0 <= t < n for t in targets), "target qubit out of range"
    assert all(0 <= c < n for c in controls), "control qubit out of range"
    assert set(targets).isdisjoint(controls), "targets and controls overlap"
    assert len(controls) == len(control_states), "controls/controlStates length mismatch"


def rightshift(qc: QuantumCircuit, targets: list[int],
               controls: list[int] | None = None,
               control_states: list[int] | None = None) -> QuantumCircuit:
    controls = [] if controls is None else list(controls)
    control_states = [1] * len(controls) if control_states is None else list(control_states)
    targets = list(targets)

    _validate_shift_args(qc, targets, controls, control_states)
    if not targets:
        return qc

    for i in range(1, len(targets)):
        ctrl = controls + targets[i:]
        states = control_states + [0] * len(targets[i:])
        mcx_with_states(qc, ctrl, targets[i - 1], states)

    if not controls:
        qc.x(targets[-1])
    else:
        mcx_with_states(qc, controls, targets[-1], control_states)

    return qc


def leftshift(qc: QuantumCircuit, targets: list[int],
              controls: list[int] | None = None,
              control_states: list[int] | None = None) -> QuantumCircuit:
    controls = [] if controls is None else list(controls)
    control_states = [1] * len(controls) if control_states is None else list(control_states)
    targets = list(targets)

    _validate_shift_args(qc, targets, controls, control_states)
    if not targets:
        return qc

    for i in range(len(targets) - 1):
        ctrl = controls + targets[i + 1:]
        states = control_states + [1] * len(targets[i + 1:])
        mcx_with_states(qc, ctrl, targets[i], states)

    if not controls:
        qc.x(targets[-1])
    else:
        mcx_with_states(qc, controls, targets[-1], control_states)

    return qc


def init_circ(n: int) -> QuantumCircuit:
    nq = 3 * n + 2
    qc = QuantumCircuit(nq, name="Init")

    for i in range(n):
        qc.swap(i, 2 * n + 2 + i)

    qc.h(n)

    for i in range(n):
        qc.cswap(n, n + 1 + i, 2 * n + 2 + i)

    return qc


def Sel_gate(n: int, as_block: bool = False):
    assert n > 0, "n must be greater than zero."

    tot_wires = 2 * n + 2
    qc = QuantumCircuit(tot_wires, name="Z_sel")

    for i in range(n + 1):
        targets = list(range(tot_wires - n - 1, tot_wires - n + i))
        rightshift(qc, targets, controls=[i], control_states=[1])

    return qc.to_gate(label="Z_{sel}") if as_block else qc


def Fix_gate(n: int, as_block: bool = False):
    assert n > 0, "n must be greater than zero."

    nq = 2 * n + 2
    qc = QuantumCircuit(nq, name="Z_fix")

    for i in range(n):
        targets = list(range(nq - n, nq - i))
        controls = [n - i, n + 1]
        control_states = [1, 0]
        leftshift(qc, targets, controls, control_states)

    return qc.to_gate(label="Z_{fix}") if as_block else qc


def Flip_gate(n: int) -> QuantumCircuit:
    nq = 3 * n + 2
    qc = QuantumCircuit(nq, name="Flip")

    for i in range(n):
        qc.cx(i, n + i + 1)

    return qc


def full_circuit(n: int) -> QuantumCircuit:
    nq = 3 * n + 2
    qc = QuantumCircuit(nq, name="full_circuit")

    # Init acts on the full register.
    qc.compose(init_circ(n), qubits=[qmap(i, nq) for i in range(nq)], inplace=True)
    qc.barrier()

    # Sel acts on qclab qubits n .. 3n+1 (2n+2 wires).
    qc.compose(Sel_gate(n), qubits=[qmap(n + i, nq) for i in range(2 * n + 2)], inplace=True)
    qc.barrier()

    # Fix acts on the same 2n+2-wire block.
    qc.compose(Fix_gate(n), qubits=[qmap(n + i, nq) for i in range(2 * n + 2)], inplace=True)
    qc.h(2 * n + 1)
    qc.barrier()

    # Flip acts on the full register.
    qc.compose(Flip_gate(n), qubits=[qmap(i, nq) for i in range(nq)], inplace=True)
    for i in range(n):
        qc.h(3 * n + 1 - i)
    qc.barrier()

    return qc


def random_normalized_state(dim: int, rng: np.random.Generator) -> np.ndarray:
    """Return a random real-valued state vector of the given dimension, L2-normalized."""
    v = rng.random(dim)
    return v / np.linalg.norm(v)





def prepare_uv(n: int) -> tuple[QuantumCircuit, QuantumCircuit, np.ndarray, np.ndarray]:
    """
    Exact O(n)-gate, depth O(1) preparation of:
      u_i = exp(theta * i / N)   and   v_i = exp(-theta * i / N)
    for i = 0, ..., N-1 where N = 2^n.

    Returns:
        circ_u: QuantumCircuit for U
        circ_v: QuantumCircuit for V
        u: Statevector of circ_u (first column of the unitary)
        v: Statevector of circ_v (first column of the unitary)
    """
    theta = 3.0
    N = 2 ** n

    # Initialize separate Qiskit circuits
    circ_u = QuantumCircuit(n)
    circ_v = QuantumCircuit(n)

    for k in range(n):
        # Weight of qubit k (MSB-first mapping)
        w = 2 ** (n - 1 - k)

        # Calculate rotation angles
        theta_u = 2 * np.arctan(np.exp(theta * w / N))
        theta_v = np.pi - theta_u

        # Apply Y-rotations to both circuits
        circ_u.ry(theta_u, k)
        circ_v.ry(theta_v, k)

    # Get the state vectors (equivalent to circ.matrix(:, 1) in MATLAB)
    u = Statevector.from_instruction(circ_u).data
    v = Statevector.from_instruction(circ_v).data

    return circ_u, circ_v, u, v




def analyze_decomposed_circuit(qc: QuantumCircuit) -> tuple[dict[str, int], int, int]:
    """
    Filters out State Preparation, decomposes the entire circuit
    into the universal Clifford+T + SWAP + Ry basis.
    """
    # 1. Filter out StatePreparation and Initialize operations
    clean_qc = QuantumCircuit(*qc.qregs, *qc.cregs)
    for instruction in qc.data:
        op_name = instruction.operation.name
        if isinstance(instruction.operation, StatePreparation) or op_name == "initialize":
            continue
        clean_qc.append(instruction)

    # 2. Transpile using the universal basis set *including* 'ry'
    universal_basis = ['cx', 'h', 'x', 't', 'tdg', 's', 'sdg', 'swap', 'ry']

    #decomposed_qc = transpile(
    #    clean_qc,
    #    basis_gates=universal_basis,
    #    optimization_level=0
    #)

    decomposed_qc=clean_qc

    # 3. Collect gate counts (excluding non-physical barriers)
    counts = dict(decomposed_qc.count_ops())
    counts.pop('barrier', None)

    # 4. Get the depth of the decomposed circuit
    depth = decomposed_qc.depth()
    size = decomposed_qc.size()

    return counts, depth, size




def mottonen_state_prep_qiskit_old(psi: np.ndarray, offset: int = 0) -> QuantumCircuit:
    """
    Möttönen et al. state preparation as a Qiskit circuit using Ry/Rz
    uniformly-controlled rotations, CNOTs, and a global-phase correction.

    References:
        Möttönen et al., arXiv:quant-ph/0407010
    """
    # 1. Validation & Normalization
    psi = np.asarray(psi, dtype=complex).flatten()
    N = len(psi)
    n = int(np.round(np.log2(N)))

    if 2 ** n != N:
        raise ValueError(f"Length of psi must be a power of two (got {N}).")

    nrm = np.linalg.norm(psi)
    if np.abs(nrm - 1.0) > 1e-8:
        print(f"Warning: ||psi|| = {nrm:.6g} != 1; normalizing state vector.")
        psi = psi / nrm

    # 2. Magnitude/Phase Tree Recursion
    r_k = np.abs(psi)
    phi_k = np.angle(psi)

    ry_angles = [None] * n
    rz_angles = [None] * n

    for k in range(n - 1, -1, -1):  # equivalent to n downto 1 (0-indexed)
        # Reshape to (2, -1) similar to MATLAB's reshape(r_k, 2, [])
        R = r_k.reshape(2, -1)
        Phi = phi_k.reshape(2, -1)

        # Note: MATLAB index 2 matches Python index 1, and 1 matches 0
        ry_angles[k] = 2 * np.arctan2(R[1, :], R[0, :])
        rz_angles[k] = Phi[1, :] - Phi[0, :]

        r_k = np.hypot(R[0, :], R[1, :])
        phi_k = (Phi[0, :] + Phi[1, :]) / 2.0

    global_phase = phi_k[0]

    # 3. Create the Qiskit Circuit
    circuit = QuantumCircuit(n + offset)

    # Helper function for recursive uniformly-controlled rotations
    def add_multiplexed_rotation(gate_type: str, target: int, controls: list[int], angles: np.ndarray):
        s = len(controls)
        if s == 0:
            if gate_type == 'ry':
                circuit.ry(angles[0], target)
            elif gate_type == 'rz':
                circuit.rz(angles[0], target)
            return

        half = len(angles) // 2
        a = angles[:half]  # top control == 0
        b = angles[half:]  # top control == 1

        a_plus = (a + b) / 2.0
        a_minus = (a - b) / 2.0

        c1 = controls[0]
        rest = controls[1:]

        add_multiplexed_rotation(gate_type, target, rest, a_plus)
        circuit.cx(c1, target)
        add_multiplexed_rotation(gate_type, target, rest, a_minus)
        circuit.cx(c1, target)

    for k in range(n):
        target = k + offset
        controls = list(range(target - 1, offset - 1, -1))  # descending, not range(offset, target)
        add_multiplexed_rotation('ry', target, controls, ry_angles[k])

    for k in range(n):
        target = k + offset
        controls = list(range(target - 1, offset - 1, -1))  # same fix here
        add_multiplexed_rotation('rz', target, controls, rz_angles[k])

    # 4. Mandatory Global-Phase Correction
    # Evaluates to e^{i * global_phase} * I on the entire register
    circuit.p(global_phase, offset)
    circuit.x(offset)
    circuit.p(global_phase, offset)
    circuit.x(offset)

    return circuit


def mottonen_state_prep_qiskit(psi: np.ndarray, offset: int = 0) -> QuantumCircuit:
    """
    Möttönen et al. state preparation implemented using the explicit Gray code
    and transformation matrix multiplication architecture inspired by the Colab notebook.

    References:
        Möttönen et al., arXiv:quant-ph/0407010
    """
    # 1. Validation & Normalization
    psi = np.asarray(psi, dtype=complex).flatten()
    N = len(psi)
    n = int(np.round(np.log2(N)))

    if 2 ** n != N:
        raise ValueError(f"Length of psi must be a power of two (got {N}).")

    nrm = np.linalg.norm(psi)
    if np.abs(nrm - 1.0) > 1e-8:
        psi = psi / nrm

    # Create the base circuit to build the disentangling sequence
    qc = QuantumCircuit(n + offset)

    if n == 0:
        return qc

    # --- Core architectural sub-routines from the Colab implementation ---
    def get_gray_code(bits: int) -> list[str]:
        gray_list = []
        for i in range(1 << bits):
            gray_binary = bin(i ^ (i >> 1))[2:]
            gray_binary = gray_binary.zfill(bits)
            gray_list.append(gray_binary)
        return gray_list

    def g(i: int) -> int:
        return i ^ (i >> 1)

    def get_m_matrix_entry(i: int, j: int, k_num: int) -> float:
        bitwise_dot_product = bin(j & g(i)).count('1')
        return (2 ** (-k_num)) * ((-1) ** bitwise_dot_product)

    current_r = np.abs(psi)
    current_phi = np.angle(psi)

    all_ry_thetas = []
    all_rz_thetas = []

    # 2. Compute tree angles and map them using the explicit M matrix
    for k in range(n - 1, -1, -1):
        size = 1 << k
        alphas_y = np.zeros(size)
        alphas_z = np.zeros(size)

        R = current_r.reshape(2, -1)
        Phi = current_phi.reshape(2, -1)

        for j in range(size):
            alphas_y[j] = 2 * np.arctan2(R[1, j], R[0, j])
            alphas_z[j] = Phi[1, j] - Phi[0, j]

        current_r = np.hypot(R[0, :], R[1, :])
        current_phi = (Phi[0, :] + Phi[1, :]) / 2.0

        # Explicitly construct and multiply the conversion matrix M
        if size == 1:
            thetas_y = alphas_y
            thetas_z = alphas_z
        else:
            matrix_M = np.zeros((size, size))
            for i in range(size):
                for j in range(size):
                    matrix_M[i][j] = get_m_matrix_entry(i, j, k)

            thetas_y = np.dot(matrix_M, alphas_y)
            thetas_z = np.dot(matrix_M, alphas_z)

        all_ry_thetas.append((k, thetas_y))
        all_rz_thetas.append((k, thetas_z))

    global_phase = current_phi[0]

    # 3. Apply the Gray-code multiplexed gate sequence
    def apply_gate_sequence(q_circuit: QuantumCircuit, gate_type: str,
                            theta_angles: np.ndarray, num_controls: int, target_qubit: int):
        size = len(theta_angles)
        if num_controls == 0:
            if gate_type == 'ry':
                q_circuit.ry(theta_angles[0], target_qubit)
            elif gate_type == 'rz':
                q_circuit.rz(theta_angles[0], target_qubit)
            return

        gray_code = get_gray_code(num_controls)
        control_gates = []
        for i in range(len(gray_code)):
            c_val = int(gray_code[i], 2) ^ int(gray_code[(i + 1) % len(gray_code)], 2)
            control_gates.append(int(np.log2(c_val)))

        for i in range(size):
            if gate_type == 'ry':
                q_circuit.ry(theta_angles[i], target_qubit)
            elif gate_type == 'rz':
                q_circuit.rz(theta_angles[i], target_qubit)

            # Map control line dynamically using offset shifts
            ctrl_index = control_gates[i] + offset
            q_circuit.cx(ctrl_index, target_qubit)

    # Build the disentangling operations (Rz then Ry)
    for k, thetas_z in all_rz_thetas:
        target = k + offset
        apply_gate_sequence(qc, 'rz', thetas_z, k, target)

    for k, thetas_y in all_ry_thetas:
        target = k + offset
        apply_gate_sequence(qc, 'ry', thetas_y, k, target)

    # 4. Correctly invert the disentangling circuit to achieve State Preparation
    qc_inv = qc.inverse()
    qc_inv.global_phase = global_phase

    return qc_inv
if __name__ == "__main__":
    # Define the range of n we want to analyze
    # Keeping upper limit lower if RUN_VERIFICATION=True, as computing operators
    # for large registers is slow.
    n_values = list(range(2, 11))
    N_values = [2 ** n for n in n_values]
    #print("NVALS")
    #print(N_values)
    theta_param = 3.0

    # Track metrics per method
    depths = {"mottonen": [], "exponential": [], "fable": []}
    sizes = {"mottonen": [], "exponential": [], "fable": []}
    gate_data = {
        "mottonen": defaultdict(list),
        "exponential": defaultdict(list),
        "fable": defaultdict(list),
    }

    tracked_gates = ['cx', 't', 'tdg', 'h', 'swap', 'x', 's', 'sdg', 'ry']

    # 2. Sequential Data Gathering
    for n in n_values:
        N = 2 ** n

        # --------------------------------------------------------
        # Build the single reference matrix S from the exact
        # geometric u, v vectors (shared target across all methods)
        # --------------------------------------------------------
        circ_u_exp, circ_v_exp, u_vec, v_vec = prepare_uv(n)

        u_v_prime = np.outer(u_vec, np.conj(v_vec))
        v_u_prime = np.outer(v_vec, np.conj(u_vec))
        S_tril = np.tril(u_v_prime)
        S_triu = np.triu(v_u_prime, k=1)
        S = S_tril + S_triu

        circ = full_circuit(n)
        nq_circ = circ.num_qubits
        qubits_for_u = list(range(nq_circ - n, nq_circ))
        qubits_for_v = list(range(nq_circ - (2 * n + 1), nq_circ - n - 1))

        # --------------------------------------------------------
        # 1. Möttönen state prep + full_circuit
        # --------------------------------------------------------
        circ_u_mot = mottonen_state_prep_qiskit(u_vec, offset=0)
        circ_v_mot = mottonen_state_prep_qiskit(v_vec, offset=0)

        mot_counts, mot_depth, mot_size = analyze_decomposed_circuit(circ_u_mot)
        #print("MOT DEPTH:", mot_depth)
        #print("MOT SIZE:", mot_size)
        #print(circ_u_mot)


        combined_mot = QuantumCircuit(nq_circ)
        combined_mot.compose(circ_u_mot, qubits=qubits_for_u, inplace=True)
        combined_mot.compose(circ_v_mot, qubits=qubits_for_v, inplace=True)
        combined_mot.barrier()
        combined_mot.compose(circ, inplace=True)

        mot_counts, mot_depth, mot_size = analyze_decomposed_circuit(combined_mot)
        depths["mottonen"].append(mot_depth)
        sizes["mottonen"].append(mot_size)
        #print(combined_mot)
        # --------------------------------------------------------
        # 2. Exact exponential state prep (prepare_uv) + full_circuit
        # --------------------------------------------------------

        combined_exp = QuantumCircuit(nq_circ)
        combined_exp.compose(circ_u_exp, qubits=qubits_for_u, inplace=True)
        combined_exp.compose(circ_v_exp, qubits=qubits_for_v, inplace=True)
        combined_exp.barrier()
        combined_exp.compose(circ, inplace=True)

        exp_counts, exp_depth, exp_size = analyze_decomposed_circuit(circ_u_exp)


        #print(circ_u_exp)


        exp_counts, exp_depth, exp_size = analyze_decomposed_circuit(combined_exp)
        depths["exponential"].append(exp_depth)
        sizes["exponential"].append(exp_size)

        #print("EXP DETPHT: ", exp_depth)
        #print("EXP size: ", exp_size)

        #print(combined_exp)
        # --------------------------------------------------------
        # 3. FABLE block encoding of the same S
        # --------------------------------------------------------
        fable_circ, alpha = fable(S, 0)
        fable_counts, fable_depth, fable_size = analyze_decomposed_circuit(fable_circ)
        depths["fable"].append(fable_depth)
        sizes["fable"].append(fable_size)
        #print("FABLE DEPHT: ", fable_depth)
        #print("FABLE size: ", fable_size)
        #print(fable_circ)
        # --------------------------------------------------------
        # OPTIONAL VERIFICATION
        # --------------------------------------------------------
        if RUN_VERIFICATION:
            print(f"\n--- Running verification for n = {n} (N = {N}) ---")

            mot_unitary = Operator(combined_mot).data
            mot_extracted_S = mot_unitary[:N, :N] * 2 * np.sqrt(N)
            mot_error = np.linalg.norm(S - mot_extracted_S)

            exp_unitary = Operator(combined_exp).data
            exp_extracted_S = exp_unitary[:N, :N] * 2 * np.sqrt(N)
            exp_error = np.linalg.norm(S - exp_extracted_S)

            fable_unitary = Operator(fable_circ).data
            fable_extracted_S = fable_unitary[:N, :N] * alpha * N
            fable_error = np.linalg.norm(S - fable_extracted_S)

            print(f"  [Möttönen]    Difference from mathematical S: {mot_error:.2e}")
            print(f"  [Exponential] Difference from mathematical S: {exp_error:.2e}")
            print(f"  [FABLE]       Difference from mathematical S: {fable_error:.2e}")

            assert mot_error < 1e-8, f"Möttönen circuit block encoding failed for n={n}!"
            assert exp_error < 1e-8, f"Exponential circuit block encoding failed for n={n}!"
            assert fable_error < 1e-8, f"FABLE circuit block encoding failed for n={n}!"
            print("  -> Verification Successful! All circuits encode S correctly.\n")

        # Record gate counts for all three
        for gate in tracked_gates:
            gate_data["mottonen"][gate].append(mot_counts.get(gate, 0))
            gate_data["exponential"][gate].append(exp_counts.get(gate, 0))
            gate_data["fable"][gate].append(fable_counts.get(gate, 0))

        print(f"n = {n} (N = {N}) | Möttönen depth: {mot_depth} | "
              f"Exponential depth: {exp_depth} | FABLE depth: {fable_depth}")
        print(f"n = {n} (N = {N}) | Möttönen size: {mot_size} | "
              f"Exponential size: {exp_size} | FABLE size: {fable_size}")
        #print("MOTONNEN DEPTH: ", circ_u_mot.depth())
        #print("MOTONNEN SIZE: ", circ_u_mot.size())

    # --- 3. Plotting: log-log comparison of decomposed circuit depth and size vs N ---
    N_arr = np.asarray(N_values, dtype=float)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    method_styles = {
        "mottonen": dict(fmt='o-', label='This work (Möttönen)'),
        "exponential": dict(fmt='s-.', label='This work (Exponential)'),
        "fable": dict(fmt='x--', label='FABLE'),
    }

    # --- Left subplot: decomposed circuit depth ---
    for method, style in method_styles.items():
        ax1.loglog(N_values, depths[method], style['fmt'], linewidth=2, markersize=7,
                   label=style['label'])

    depth_guide_end = depths["fable"][0] * (N_arr[-1] / N_arr[0]) ** 2
    ax1.loglog([N_values[0], N_values[-1]], [depths["fable"][0], depth_guide_end],
               ':', color=(.6, .6, .6), label=r'$\propto N^2$')
    ax1.set_xlabel('$N$')
    ax1.set_ylabel('Circuit depth')
    ax1.set_title('Non decomposed circuit depth')
    ax1.grid(True, which='both', linestyle=':', alpha=0.6)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.set_xticks(N_values)
    ax1.set_xticklabels([str(v) for v in N_values])
    ax1.minorticks_off()

    # --- Right subplot: decomposed circuit size (total gate count) ---
    for method, style in method_styles.items():
        ax2.loglog(N_values, sizes[method], style['fmt'], linewidth=2, markersize=7,
                   label=style['label'])

    size_guide_end = sizes["fable"][0] * (N_arr[-1] / N_arr[0]) ** 2
    ax2.loglog([N_values[0], N_values[-1]], [sizes["fable"][0], size_guide_end],
               ':', color=(.6, .6, .6), label=r'$\propto N^2$')
    ax2.set_xlabel('$N$')
    ax2.set_ylabel('Circuit size (gate count)')
    ax2.set_title('Non decomposed circuit size')
    ax2.grid(True, which='both', linestyle=':', alpha=0.6)
    ax2.legend(loc='upper left', fontsize=9)
    ax2.set_xticks(N_values)
    ax2.set_xticklabels([str(v) for v in N_values])
    ax2.minorticks_off()

    plt.tight_layout()
    plt.show()



