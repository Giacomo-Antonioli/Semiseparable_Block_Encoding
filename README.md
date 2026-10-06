# Semiseparable Block Encoding

Code accompanying the paper **TODO: title, authors, arXiv link**.

It builds quantum circuits that block-encode rank-1 semiseparable matrices

- symmetric: `S = tril(u*v') + triu(v*u', 1)`, subnormalization `2*sqrt(N)`
- unsymmetric: `S = tril(u*v') + triu(x*y', 1)`, subnormalization `sqrt(N)*(‖u‖‖v‖ + ‖x‖‖y‖)`

using `3n + 2` qubits for `N = 2^n`. The scripts in `experiments/` reproduce
every numerical result in the paper, including the comparison against
[FABLE](https://github.com/QuantumComputingLab/fable).

## Requirements

- MATLAB (tested with R2025b)
- [QCLAB](https://github.com/QuantumComputingLab/qclab) on the MATLAB path.
  **TODO:** this code uses QCLAB features that are not in the upstream release
  (`qgates.CSWAP`, `QCircuit.InitializeStateVector`, `QCircuit.getGateCount`,
  `getDepth`, `getQubitCount`, `getControlledStats`). Point this to the fork or
  release that contains them.
- [FABLE](https://github.com/QuantumComputingLab/fable), included as a git
  submodule (only needed for the FABLE comparisons)
- Optional, for the Qiskit port in `python/`: Python 3.10+ and the packages in
  `python/requirements.txt`

```bash
git clone --recurse-submodules <this-repo-url>
```

## Quickstart

From the repository root in MATLAB:

```matlab
setup_paths                         % adds src/, experiments/utils, plotting/, FABLE

N = 8;
u = randn(N,1); u = u/norm(u);
v = randn(N,1); v = v/norm(v);
S = tril(u*v') + triu(v*u', 1);

circ = build_semiseparable_circuit(u, v);   % qclab.QCircuit on 3*log2(N)+2 qubits
U = circ.matrix;
norm(2*sqrt(N)*U(1:N,1:N) - S)              % ~1e-16
```

See also `examples/` and `tests/`.

## Repository layout

```
setup_paths.m        add everything to the MATLAB path
src/
  circuits/          block-encoding circuits and their building blocks
  state_prep/        state preparation of the generators
  generators/        test generators (random, tridiagonal-inverse, exponential/OU)
experiments/         one script per experiment in the paper (+ utils/)
plotting/            PGFPlots/TikZ export of the saved results
examples/            small demos
tests/               correctness checks
python/              Qiskit implementation
results/             data and figures used in the paper
external/fable/      FABLE (git submodule)
```

### Circuits (`src/circuits`)

| Function | Description |
|---|---|
| `build_semiseparable_circuit(v,u)` | Symmetric case; generators prepared with a dense unitary (`state_prep`) |
| `build_semiseparable_circuit_mottonen(v,u)` | Same, with Möttönen state preparation |
| `build_semiseparable_circuit_exp(v,u)` | Same, for exponential generators prepared with `n` RY rotations (`prepareUV`) |
| `build_unsymmetric_semiseparable_circuit(v,u,x,y,nu,nv,nx,ny)` | Unsymmetric case |
| `init_circ`, `Sel_gate`, `Fix_gate`, `Flip_gate`, `leftshift`, `rightshift` | Building blocks |

## Reproducing the paper

The scripts below can be run from any folder (each one calls `setup_paths`). Each one writes to its own
folder under `results/`. New runs of the correctness and perturbation
experiments go to timestamped subfolders, which git ignores.

| Paper | Script | Output in `results/` | LaTeX figure |
|---|---|---|---|
| TODO Fig. – symmetric correctness | `experiments/symmetric_correctness.m` | `symmetric_correctness/` | `plotting/latex_correctness_plot.m` on `exp_1/runs` |
| TODO Fig. – symmetric perturbation bound | `experiments/symmetric_perturbation.m` | `symmetric_perturbation/` | `plotting/latex_perturbation_plot.m` on `runs` |
| TODO Fig. – unsymmetric correctness | `experiments/unsymmetric_correctness.m` | `unsymmetric_correctness/` | printed by the script |
| TODO Fig. – unsymmetric perturbation bound | `experiments/unsymmetric_perturbation.m` | `unsymmetric_perturbation/` | written by the script |
| TODO Fig. – subnormalization vs FABLE | `experiments/subnormalization_vs_fable.m` | `subnormalization_vs_fable/` | written by the script |
| TODO Tab. – gate count / depth vs FABLE | `experiments/resources_vs_fable.m` | `resources_vs_fable/` | – |

Problem sizes, the number of trials and the random seeds are set in the
`cfg` block at the top of each script.

## Citation

TODO: BibTeX entry.

## License

TODO. FABLE is distributed under its own license (see `external/fable`).
