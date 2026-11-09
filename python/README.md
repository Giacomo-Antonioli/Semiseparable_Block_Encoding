# Semiseparable block encoding: Qiskit

Qiskit implementation of the block encoding of rank-1 semiseparable matrices.
This branch holds only the working code. Tests, experiments and the FABLE
comparison live on the `tests` branch of the same repository.

```
pip install -r requirements.txt
python semiseparable_qiskit.py        # small correctness demo
```

```python
from semiseparable_qiskit import full_circuit, mottonen_state_prep_qiskit
circ = full_circuit(n)                 # 3n + 2 qubits
```

The module contains the circuit (`full_circuit` and its building blocks),
Mottonen state preparation and `analyze_decomposed_circuit` (gate counts in a
universal basis). It does not depend on FABLE or matplotlib.
