# Block encoding of semiseparable matrices via their generators

Code and data for the paper

> G. Antonioli, G. M. Del Corso, R. Van Beeumen,
> **Block encoding of semiseparable matrices via their generators**,
> [arXiv:2610.13228](https://arxiv.org/abs/2610.13228) (submitted to SIAM Journal on Matrix Analysis and Applications).

This repository exists to reproduce the numerical results of that paper. It
contains the quantum circuits that block-encode generator-representable
(one-pair) semiseparable matrices, in QCLAB (MATLAB) and in Qiskit (Python),
and the scripts and saved data behind every figure and numerical claim of
Section 8 (Numerical experiments).

- symmetric case: `S = tril(u*v') + triu(v*u', 1)`, subnormalization `2*sqrt(N)`
- unsymmetric case: `S = tril(u*v') + triu(x*y', 1)`, subnormalization `sqrt(N)*(‖u‖‖v‖ + ‖x‖‖y‖)`

with `3n + 2` qubits (`2n + 2` ancillas) for `N = 2^n`.

**Version of the submitted paper:** the git tag `TODO-tag` (see *Releases* below).

## Requirements

| Component | Version | Needed for |
|---|---|---|
| MATLAB | R2025b | everything in `experiments/`, `tests/`, `examples/` |
| [QCLAB](https://github.com/QuantumComputingLab/qclab) fork `Giacomo-Antonioli/qclab`, branch `feature/qsemisep-extras`: QCLAB v1.2.0 plus `CSWAP`, `QGate3`, `InitializeStateVector`, gate-count and depth helpers | submodule `external/qclab`, commit `2e253b2` | MATLAB code |
| [FABLE](https://github.com/QuantumComputingLab/fable) (`fable-circuits` 1.0.2) | submodule `external/fable`, commit `6c1766d` | FABLE comparison (MATLAB and Python) |
| Python | 3.12.3 | `python/` |
| Qiskit, NumPy, Matplotlib | 2.5.1, 2.5.1, 3.11.1 (pinned in `python/requirements.txt`) | `python/` |

```bash
git clone --recurse-submodules https://github.com/Giacomo-Antonioli/Semiseparable_Block_Encoding
cd Semiseparable_Block_Encoding
pip install -r python/requirements.txt          # Python part only
pip install external/fable/py                   # FABLE for Python, from the submodule
```

In MATLAB, from the repository root, run `setup_paths` once per session. It
puts `src/`, `experiments/utils/`, `plotting/`, FABLE and the QCLAB fork on the
path. The stock QCLAB v1.2.0 is **not** enough (it lacks `qgates.CSWAP` and the
other additions listed above).

## Quickstart

```matlab
setup_paths
N = 8;
u = randn(N,1); u = u/norm(u);
v = randn(N,1); v = v/norm(v);
S = tril(u*v') + triu(v*u', 1);
circ = build_semiseparable_circuit(u, v);   % qclab.QCircuit on 3*log2(N)+2 qubits
U = circ.matrix;
norm(2*sqrt(N)*U(1:N,1:N) - S)              % ~1e-16
```

The automated checks (`tests/`, also run by CI with `runtests`) are:

```matlab
setup_paths; results = runtests('tests'); table(results)
```

## File-by-file description

### Circuits and building blocks: `src/circuits/`

| File | Description |
|---|---|
| `build_semiseparable_circuit.m` | symmetric case; generators loaded by a dense unitary (`state_prep`) |
| `build_semiseparable_circuit_mottonen.m` | symmetric case; generators loaded by Möttönen state preparation |
| `build_semiseparable_circuit_exp.m` | symmetric case; dedicated loader for the Exponential generators (`n` R_y rotations, `prepareUV`) |
| `build_unsymmetric_semiseparable_circuit.m` | unsymmetric case (Section 6 of the paper) |
| `init_circ.m`, `Sel_gate.m`, `Fix_gate.m`, `Flip_gate.m`, `leftshift.m`, `rightshift.m` | shift, selective-shift and swap building blocks (Section 4) |

### State preparation: `src/state_prep/`

`state_prep.m`, `generate_sp_block.m` (dense unitary whose first column is the state), `prepareUV.m` (R_y loader for the exponential generators).

### Test generators: `src/generators/`

| File | Description |
|---|---|
| `local_generate_uv.m` | the three families of the paper: type 1 Random, type 2 Tridiagonal, type 3 Exponential |
| `ss_generators_from_tridiag.m` | generators from the first and last columns of the inverse of an irreducible, diagonally dominant symmetric tridiagonal matrix |
| `gen_ornstein_uhlenbeck.m` | exponential generators (Ornstein–Uhlenbeck covariance) |
| `generateUnsymm.m` | random unsymmetric generators `u, v, x, y` |

### Experiments: `experiments/`

One script per experiment (each calls `setup_paths`, fixes its seed with `rng(...)` and has a `cfg` block at the top with sizes, ε values and number of trials). Helpers are in `experiments/utils/`.

| File | Experiment |
|---|---|
| `symmetric_correctness.m` | error `‖S − S̃‖₂` for the three families (Section 8.1) |
| `unsymmetric_correctness.m` | same for random unsymmetric matrices (Section 8.1) |
| `symmetric_perturbation.m` | error versus ε and versus N against the bound `2√2 ε` (Section 8.2, Fig. `fig:err_combined`) |
| `unsymmetric_perturbation.m` | same for the unsymmetric case, bound `2√(w_lo²+w_up²) ε` |
| `subnormalization_vs_fable.m` | subnormalization `α/‖S‖₂` of our circuits and of FABLE, N = 4 … 1024 |
| `resources_vs_fable.m` | gate count and depth in QCLAB, ours versus FABLE (small N, needs the full unitary) |

### Plotting: `plotting/`

`latex_correctness_plot.m`, `latex_perturbation_plot.m`, `generate_latex_figure.m`, `coords_str.m` write the PGFPlots coordinates used in the paper from the saved results.

### Python / Qiskit: `python/`

| File | Description |
|---|---|
| `semiseparable_qiskit.py` | Qiskit implementation of the circuits (Möttönen and exponential loaders) |
| `semiseparable_vs_fable_qiskit.py` | depth and size of our circuits and of FABLE for N = 4 … 1024 (Section 8.3, Fig. `fig:fable`) |
| `requirements.txt` | Python dependencies |

### Data: `results/`

Output of the experiment scripts that was used in the paper (`.mat`, `.csv`, `.tex`, figures). New runs go to timestamped subfolders, which git ignores.

### Other

`examples/` small demos, `tests/` function-based tests, `external/` submodules (QCLAB fork, FABLE), `scripts/make_supplement_zip.sh` builds the zip for the SIMAX Supplementary Materials, `LICENSE` (BSD 3-Clause).

## Reproducing the paper

Random seeds: every test, example and experiment script sets a fixed seed (`rng(42)`, except `unsymmetric_correctness.m` with `rng(46)`). The Python FABLE script is deterministic (no random numbers).

Generators (Section 8): Random is a normalized Gaussian vector pair; Tridiagonal uses the first and last columns of the inverse of a random irreducible, diagonally dominant symmetric tridiagonal matrix; Exponential is `u_i = exp(3 i/N)`, `v_i = exp(-3 i/N)` normalized.

| Paper result | Command | Parameters in the paper | Output |
|---|---|---|---|
| Exactness: error below 1e-15 for N = 4 … 128 (Sec. 8.1) | `experiments/symmetric_correctness.m` | `n = 2:7`, three families | `results/symmetric_correctness/` |
| Exactness, unsymmetric (Sec. 8.1) | `experiments/unsymmetric_correctness.m` | `n = 2:7` | `results/unsymmetric_correctness/` |
| Error versus ε, N = 64, ε from 1e-14 to 1e-6 (Fig. `fig:first`) | `experiments/symmetric_perturbation.m` | `n_fixed = 6` | `results/symmetric_perturbation/epsilon_summary.csv`, `plots/error_vs_epsilon.png` |
| Error versus N, ε = 1e-6, N = 4 … 128 (Fig. `fig:second`) | same script | `n_range = 2:7`, `eps_fixed = 1e-6` | `results/symmetric_perturbation/N_summary.csv`, `plots/error_vs_N.png` |
| Unsymmetric tightness, N up to 32 (Sec. 8.2) | `experiments/unsymmetric_perturbation.m` | `n_range = 2:5` | `results/unsymmetric_perturbation/` |
| FABLE comparison, depth and size, N = 4 … 1024, fitted exponents, crossover (Fig. `fig:fable`, Sec. 8.3) | `python python/semiseparable_vs_fable_qiskit.py` (options `--nmax`, `--out-dir`, `--show`) | `n = 2:10`, θ = 3, transpiled with `optimization_level=0` to `{cx, h, x, t, tdg, s, sdg, swap, ry}` | `results/fable_python/fable_comparison.csv` (raw and transpiled depth/size), `fits.txt`, `fable_comparison.png` |
| Subnormalization versus FABLE | `experiments/subnormalization_vs_fable.m` | `dims = 2:10` | `results/subnormalization_vs_fable/` |

LaTeX/PGFPlots coordinates for the figures: `plotting/latex_correctness_plot.m` (on `results/symmetric_correctness/exp_1/runs`) and `plotting/latex_perturbation_plot.m` (on `results/symmetric_perturbation/runs`).

### Known gaps (to be closed before the submitted version is tagged)

- The `cfg` blocks of the MATLAB scripts currently hold reduced debugging values (for example `symmetric_correctness.m` has `n_range = 1:3`, `symmetric_perturbation.m` has `n_fixed = 3`). The paper values are those in the table above and must be set before re-running. The committed files in `results/` should be checked against the paper figures.
- Section 8.3 states the gate set as {CNOT, X, T, T†, S, S†, SWAP, R_y}. The script also keeps `h` (FABLE circuits contain Hadamards, which cannot be translated without it). The paper numbers for our circuits are reproduced with `h` in the set.
- FABLE counts: the script reproduces the paper values for N = 4, 8, 16 and differs from them for N ≥ 32 (see the checklist delivered with this README). Expected runtime and memory per script: TODO (measure on the machine used for the paper).

## Citation

```bibtex
@article{antonioli2026semiseparable,
  author  = {Antonioli, Giacomo and Del Corso, Gianna M. and Van Beeumen, Roel},
  title   = {Block encoding of semiseparable matrices via their generators},
  journal = {arXiv preprint arXiv:2610.13228},
  year    = {2026}
}
```

## License

BSD 3-Clause, see `LICENSE`. QCLAB and FABLE (submodules) are distributed under their own licenses.
