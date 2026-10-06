from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import MCXGate
from qiskit.quantum_info import Operator


def qmap(i: int, nq: int) -> int:
    """
    Convert a qclab qubit index to a Qiskit qubit index.

    qclab convention: qubit 0 is the bottom wire.
    Qiskit convention: qubit 0 is the top wire.
    """
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

from qiskit.circuit.library import StatePreparation
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


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
def demo(n: int = 2, seed: int = 0) -> None:
    """
    Build the full circuit for a given n, prepare two random normalized
    input vectors on separate qubit blocks, compose them with the circuit,
    and print the leading N columns of the resulting unitary.
    """
    rng = np.random.default_rng(seed)

    circ = full_circuit(n)
    print(f"Full circuit for n={n} ({circ.num_qubits} qubits):")
    #print(circ.draw(output="text"))

    N = 2 ** n
    nq_circ = circ.num_qubits

    v1 = random_normalized_state(N, rng)
    v2 = random_normalized_state(N, rng)
    u = v1
    v = v2

    # Calculate u*v' (outer product)
    u_v_prime = np.outer(u, np.conj(v)) # Using conj for complex conjugate transpose

    # Calculate v*u' (outer product)
    v_u_prime = np.outer(v, np.conj(u))

    # Extract lower triangular part of u*v'
    S_tril = np.tril(u_v_prime)

    # Extract upper triangular part of v*u', excluding the diagonal (k=1)
    S_triu = np.triu(v_u_prime, k=1)

    # Combine to form S
    S = S_tril + S_triu
    # Highest-indexed n qubits for v1, then skip one qubit, then the next
    # n qubits (going down) for v2.
    qubits_for_v1 = list(range(nq_circ - n, nq_circ))
    skipped_qubit = nq_circ - n - 1
    qubits_for_v2 = list(range(nq_circ - (2 * n + 1), nq_circ - n - 1))

    init_state = QuantumCircuit(nq_circ, name="StatePreparation")
    init_state.append(StatePreparation(v1), qubits_for_v1)
    init_state.append(StatePreparation(v2), qubits_for_v2)

    init_state.barrier()

    print("\nState-preparation circuit:")
    #print(init_state.draw(output="text"))

    combined = init_state.compose(circ, qubits=list(range(nq_circ)))
    print("\nCombined circuit (state prep + full_circuit):")
    #print(combined.draw(output="text"))

    # Leading N columns of the unitary implemented by the (un-initialized) circuit.
    U = Operator(combined).data*2*np.sqrt(N)
    myS=U[:N, :N]
    print(f"\nFirst {N} columns of the unitary matrix for full_circuit(n={n}):")
    print(np.round(myS, 4))
    print(np.round(S,4))
    print(np.linalg.norm(S-myS))






import numpy as np
import matplotlib

# Force Matplotlib to use a standard interactive pop-up window
# This bypasses PyCharm's "tostring_rgb" tool window bug
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from collections import defaultdict
from fable import fable  # Import FABLE block-encoding
from qiskit import transpile, QuantumCircuit
from qiskit.circuit.library import StatePreparation



def analyze_decomposed_circuit(qc: QuantumCircuit) -> tuple[dict[str, int], int,int]:
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

    decomposed_qc = transpile(
        clean_qc,
        basis_gates=universal_basis,
        optimization_level=1
    )

    # 3. Collect gate counts (excluding non-physical barriers)
    counts = dict(decomposed_qc.count_ops())
    counts.pop('barrier', None)

    # 4. Get the depth of the decomposed circuit
    depth = decomposed_qc.depth()
    size = decomposed_qc.size()

    return counts, depth,size





def mottonen_state_prep_qiskit(psi: np.ndarray, offset: int = 0) -> QuantumCircuit:
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



"""
if __name__ == "__main__":
    # Define the range of n we want to analyze (e.g., 2 to 6)
    # Note: If running verification, keep the upper limit low (n <= 5)
    # because calculating the full Operator of a 3n+2 register gets extremely slow.
    n_values = list(range(2, 8))
    N_values = [2 ** n for n in n_values]

    rng = np.random.default_rng()

    # Track metrics for your original circuit (solid lines)
    orig_depths = []
    orig_gate_data = defaultdict(list)

    # Track metrics for FABLE circuit (dashed lines)
    fable_depths = []
    fable_gate_data = defaultdict(list)

    tracked_gates = ['cx', 't', 'tdg', 'h', 'swap', 'x', 's', 'sdg', 'ry']

    # 2. Sequential Data Gathering
    for n in n_values:
        N = 2 ** n

        # --- 1. Prepare exact U and V state preparation circuits and vectors ---
        circ_u, circ_v, u_vec, v_vec = prepare_uv(n)

        # --- 2. Generate the semiseparable matrix S using u_vec and v_vec ---
        u = u_vec
        v = v_vec

        # Calculate outer products
        u_v_prime = np.outer(u, np.conj(v))
        v_u_prime = np.outer(v, np.conj(u))

        # Combine lower-triangular and strictly upper-triangular parts
        S_tril = np.tril(u_v_prime)
        S_triu = np.triu(v_u_prime, k=1)
        S = S_tril + S_triu

        # --- 3. Build FABLE Circuit using the constructed matrix S ---
        fable_circ, alpha = fable(S, 0)
        fable_counts, fable_depth = analyze_decomposed_circuit(fable_circ)
        fable_depths.append(fable_depth)

        # --- 4. Build Your Original Circuit with U and V State Prep Appended ---
        circ = full_circuit(n)
        nq_circ = circ.num_qubits

        # Define target qubit registers
        qubits_for_u = list(range(nq_circ - n, nq_circ))
        qubits_for_v = list(range(nq_circ - (2 * n + 1), nq_circ - n - 1))

        # Assemble combined circuit: [State Prep] -> [Barrier] -> [Full Circuit]
        combined_circ = QuantumCircuit(nq_circ)
        combined_circ.compose(circ_u, qubits=qubits_for_u, inplace=True)
        combined_circ.compose(circ_v, qubits=qubits_for_v, inplace=True)
        combined_circ.barrier()  # Visual boundary
        combined_circ.compose(circ, inplace=True)

        # Analyze original + state-prep circuit
        orig_counts, orig_depth = analyze_decomposed_circuit(combined_circ)
        orig_depths.append(orig_depth)

        # --- OPTIONAL VERIFICATION STEP ---
        if RUN_VERIFICATION:
            print(f"\n--- Running verification for n = {n} (N = {N}) ---")

            # Extract Original Block matrix (scaled up by 2 * sqrt(N))
            orig_unitary = Operator(combined_circ).data
            orig_extracted_S = orig_unitary[:N, :N] * 2 * np.sqrt(N)
            orig_error = np.linalg.norm(S - orig_extracted_S)

            # Extract FABLE Block matrix (scaled up by alpha)
            fable_unitary = Operator(fable_circ).data
            fable_extracted_S = fable_unitary[:N, :N] * alpha *N
            fable_error = np.linalg.norm(S - fable_extracted_S)

            print(f"  [Original] Difference from mathematical S: {orig_error:.2e}")
            print(f"  [FABLE]    Difference from mathematical S: {fable_error:.2e}")

            # Standard tolerance check
            assert orig_error < 1e-8, f"Original circuit block encoding failed for n={n}!"
            assert fable_error < 1e-8, f"FABLE circuit block encoding failed for n={n}!"
            print("  -> Verification Successful! Both circuits encode S correctly.\n")

        # Record gate counts for both
        for gate in tracked_gates:
            orig_gate_data[gate].append(orig_counts.get(gate, 0))
            fable_gate_data[gate].append(fable_counts.get(gate, 0))

        print(f"n = {n} (N = {N}) processed | Orig Depth: {orig_depth} | FABLE Depth: {fable_depth}")

    # --- 3. Plotting ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # --- Left Plot: Depth comparison ---
    ax1.plot(N_values, orig_depths, marker='o', color='royalblue', linewidth=2.5,
             linestyle='-', label="Original + Prep (Depth)")
    ax1.plot(N_values, fable_depths, marker='x', color='crimson', linewidth=2,
             linestyle='--', label="FABLE Circuit (Depth)")

    ax1.set_title("Decomposed Circuit Depth vs. $N$ (Log Scale)", fontsize=13, fontweight='bold')
    ax1.set_xlabel("$N$ ($2^n$)", fontsize=12)
    ax1.set_ylabel("Depth", fontsize=12)
    ax1.set_xscale('log', base=2)
    ax1.set_xticks(N_values)
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(fontsize=10)

    # --- Right Plot: Gate Counts comparison ---
    colors = {
        'cx': '#1f77b4', 't': '#ff7f0e', 'tdg': '#2ca02c', 'h': '#d62728',
        'swap': '#9467bd', 'x': '#8c564b', 's': '#e377c2', 'sdg': '#7f7f7f',
        'ry': '#bcbd22'
    }

    for gate in tracked_gates:
        orig_list = orig_gate_data[gate]
        fable_list = fable_gate_data[gate]

        if sum(orig_list) > 0 or sum(fable_list) > 0:
            color = colors.get(gate, None)

            if sum(orig_list) > 0:
                ax2.plot(N_values, orig_list, marker='s', color=color, linestyle='-',
                         label=f"{gate} (Original + Prep)")
            if sum(fable_list) > 0:
                ax2.plot(N_values, fable_list, marker='^', color=color, linestyle='--',
                         alpha=0.8, label=f"{gate} (FABLE)")

    ax2.set_title("Decomposed Gate Counts vs. $N$ (Log Scale)", fontsize=13, fontweight='bold')
    ax2.set_xlabel("$N$ ($2^n$)", fontsize=12)
    ax2.set_ylabel("Gate Count", fontsize=12)
    ax2.set_xscale('log', base=2)
    ax2.set_xticks(N_values)
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.grid(True, linestyle=':', alpha=0.6)

    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0, fontsize=9)

    plt.tight_layout()
    plt.show()
    
    """

if __name__ == "__main__":
    # Define the range of n we want to analyze (e.g., 2 to 5)
    # Keeping upper limit to 5 if RUN_VERIFICATION=True, as computing operators for large registers is slow.
    n_values = list(range(2, 10))
    N_values = [2 ** n for n in n_values]

    theta_param = 3.0

    # Track metrics for your original circuit (solid lines)
    orig_depths = []
    orig_gate_data = defaultdict(list)

    # Track metrics for FABLE circuit (dashed lines)
    fable_depths = []
    fable_gate_data = defaultdict(list)

    # Track total circuit size
    orig_sizes = []
    fable_sizes = []

    tracked_gates = ['cx', 't', 'tdg', 'h', 'swap', 'x', 's', 'sdg', 'ry']

    # 2. Sequential Data Gathering
    for n in n_values:
        N = 2 ** n

        # --- 1. Mathematically generate normalized geometric statevectors u and v ---
        indices = np.arange(N)
        #u_raw = np.exp(theta_param * indices / N)
        #v_raw = np.exp(-theta_param * indices / N)

        u_raw=np.random.rand(N)
        v_raw = np.random.rand(N)
        # Normalize vectors
        u_vec = u_raw / np.linalg.norm(u_raw)
        v_vec = v_raw / np.linalg.norm(v_raw)

        # --- 2. Build Möttönen state preparation circuits (n-qubit subcircuits, offset=0) ---
        circ_u = mottonen_state_prep_qiskit(u_vec, offset=0)
        circ_v = mottonen_state_prep_qiskit(v_vec, offset=0)

        # --- 3. Generate the semiseparable matrix S using u_vec and v_vec ---
        u_v_prime = np.outer(u_vec, np.conj(v_vec))
        v_u_prime = np.outer(v_vec, np.conj(u_vec))

        # Combine lower-triangular and strictly upper-triangular parts
        S_tril = np.tril(u_v_prime)
        S_triu = np.triu(v_u_prime, k=1)
        S = S_tril + S_triu

        # --- 4. Build FABLE Circuit using the constructed matrix S ---
        fable_circ, alpha = fable(S, 0)
        fable_counts, fable_depth, fable_size= analyze_decomposed_circuit(fable_circ)
        fable_depths.append(fable_depth)
        fable_sizes.append(fable_size)
        # --- 5. Build Your Original Circuit and append the Möttönen initializations ---
        circ = full_circuit(n)
        nq_circ = circ.num_qubits

        # Define target qubit registers in the main circuit
        qubits_for_u = list(range(nq_circ - n, nq_circ))
        qubits_for_v = list(range(nq_circ - (2 * n + 1), nq_circ - n - 1))

        # Assemble combined circuit: [Möttönen State Prep] -> [Barrier] -> [Full Circuit]
        combined_circ = QuantumCircuit(nq_circ)
        combined_circ.compose(circ_u, qubits=qubits_for_u, inplace=True)
        combined_circ.compose(circ_v, qubits=qubits_for_v, inplace=True)
        combined_circ.barrier()  # Visual boundary
        combined_circ.compose(circ, inplace=True)
        #print(combined_circ)
        # Analyze original + state-prep circuit
        orig_counts, orig_depth, orig_size = analyze_decomposed_circuit(combined_circ)
        orig_depths.append(orig_depth)
        orig_sizes.append(orig_size)

        # --- OPTIONAL VERIFICATION STEP ---
        if RUN_VERIFICATION:
            print(f"\n--- Running verification for n = {n} (N = {N}) ---")

            # Extract Original Block matrix (scaled up by 2 * sqrt(N))
            orig_unitary = Operator(combined_circ).data
            orig_extracted_S = orig_unitary[:N, :N] * 2 * np.sqrt(N)
            orig_error = np.linalg.norm(S - orig_extracted_S)

            # Extract FABLE Block matrix (scaled up by alpha)
            fable_unitary = Operator(fable_circ).data
            fable_extracted_S = fable_unitary[:N, :N] * alpha * N
            fable_error = np.linalg.norm(S - fable_extracted_S)

            print(f"  [Original (Möttönen)] Difference from mathematical S: {orig_error:.2e}")
            print(f"  [FABLE]              Difference from mathematical S: {fable_error:.2e}")


            # Print matrices (real part only, 2 decimals)
            def print_matrix(name, M):
                print(f"\n{name}:")
                print(np.real(M).round(2))


            print_matrix("Mathematical S", S)
            print_matrix("Original extracted S", orig_extracted_S)
            print_matrix("FABLE extracted S", fable_extracted_S)

            # Standard tolerance check
            assert orig_error < 1e-8, f"Original circuit block encoding failed for n={n}!"
            assert fable_error < 1e-8, f"FABLE circuit block encoding failed for n={n}!"
            print("  -> Verification Successful! Both circuits encode S correctly.\n")

        # Record gate counts for both
        for gate in tracked_gates:
            orig_gate_data[gate].append(orig_counts.get(gate, 0))
            fable_gate_data[gate].append(fable_counts.get(gate, 0))

        print(f"n = {n} (N = {N}) processed | Orig Depth (with Möttönen): {orig_depth} | FABLE Depth: {fable_depth}")

    # --- 3. Plotting ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # --- Left Plot: Depth comparison ---
    ax1.plot(N_values, orig_depths, marker='o', color='royalblue', linewidth=2.5,
             linestyle='-', label="Original + Möttönen Prep (Depth)")
    ax1.plot(N_values, fable_depths, marker='x', color='crimson', linewidth=2,
             linestyle='--', label="FABLE Circuit (Depth)")

    ax1.set_title("Decomposed Circuit Depth vs. $N$ (Log Scale)", fontsize=13, fontweight='bold')
    ax1.set_xlabel("$N$ ($2^n$)", fontsize=12)
    ax1.set_ylabel("Depth", fontsize=12)
    ax1.set_xscale('log', base=2)
    ax1.set_xticks(N_values)
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(fontsize=10)

    # --- Right Plot ---
    if PLOT_GATE_COUNTS:

        colors = {
            'cx': '#1f77b4', 't': '#ff7f0e', 'tdg': '#2ca02c',
            'h': '#d62728', 'swap': '#9467bd', 'x': '#8c564b',
            's': '#e377c2', 'sdg': '#7f7f7f', 'ry': '#bcbd22'
        }

        for gate in tracked_gates:
            orig_list = orig_gate_data[gate]
            fable_list = fable_gate_data[gate]

            if sum(orig_list) > 0 or sum(fable_list) > 0:
                color = colors.get(gate)

                if sum(orig_list) > 0:
                    ax2.plot(
                        N_values,
                        orig_list,
                        marker='s',
                        color=color,
                        linestyle='-',
                        label=f"{gate} (Original + Möttönen)"
                    )

                if sum(fable_list) > 0:
                    ax2.plot(
                        N_values,
                        fable_list,
                        marker='^',
                        color=color,
                        linestyle='--',
                        alpha=0.8,
                        label=f"{gate} (FABLE)"
                    )

        ax2.set_title("Decomposed Gate Counts vs. $N$", fontsize=13, fontweight='bold')
        ax2.set_ylabel("Gate Count", fontsize=12)
        ax2.legend(
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
            borderaxespad=0,
            fontsize=9
        )

    else:

        ax2.plot(
            N_values,
            orig_sizes,
            marker='o',
            linewidth=2.5,
            color='royalblue',
            label="Original + Möttönen"
        )

        ax2.plot(
            N_values,
            fable_sizes,
            marker='x',
            linewidth=2,
            linestyle='--',
            color='crimson',
            label="FABLE"
        )

        ax2.set_title("Decomposed Circuit Size vs. $N$", fontsize=13, fontweight='bold')
        ax2.set_ylabel("Circuit Size", fontsize=12)
        ax2.legend(fontsize=10)

    ax2.set_xlabel("$N$ ($2^n$)", fontsize=12)
    ax2.set_xscale('log', base=2)
    ax2.set_xticks(N_values)
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.grid(True, linestyle=':', alpha=0.6)

    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0, fontsize=9)

    plt.tight_layout()
    plt.show()