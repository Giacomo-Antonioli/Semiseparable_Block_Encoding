# Semiseparable block encoding: library

Circuit builders for the block encoding of rank-1 semiseparable matrices,
written for [QCLAB](https://github.com/QuantumComputingLab/qclab). This branch
holds only the working code, to be imported as a module: MATLAB/QCLAB in
`circuits/` and `state_prep/`, Qiskit in `python/`. Tests, experiments, the
FABLE comparison and the paper data live on the `tests` branch of the same
repository.

## Requirements

A QCLAB that contains `qgates.CSWAP` and `QCircuit.InitializeStateVector`,
`getGateCount`, `getDepth`, `getQubitCount`, `getControlledStats`. These are in
the `feature/qsemisep-extras` branch of the fork
`Giacomo-Antonioli/qclab` (QCLAB v1.2.0 plus these additions).

## Use

```matlab
addpath(genpath('path/to/this/folder'));
circ = build_semiseparable_circuit(u, v);   % qclab.QCircuit
```

## Contents

| Folder | Content |
|---|---|
| `circuits/` | `build_semiseparable_circuit`, `_mottonen`, `_exp`, `build_unsymmetric_semiseparable_circuit`, and the building blocks `init_circ`, `Sel_gate`, `Fix_gate`, `Flip_gate`, `leftshift`, `rightshift` |
| `state_prep/` | state preparation of the generators (Mottonen state preparation is `qclab.mottonen_state_prep` in the fork) |

Test generators (`generateUnsymm`, random and OU generators) are on the `tests`
branch.

## Qiskit

The Qiskit implementation is in `python/` (`semiseparable_qiskit.py`, see
`python/README.md`).
