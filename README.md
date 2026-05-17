# Orbital Coherence Function (OCF)
### A Keplerian Pre-Filter for Ensemble Orbit Determination

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Reference:** Maia, K. (2026). *The Orbital Coherence Function: A Keplerian Pre-Filter for Ensemble Orbit Determination* Independent Researcher, Marietta, Georgia, USA.

---

## Overview

This repository provides a reference implementation of the **Orbital Coherence Function (OCF)**, a modular filter for accelerating orbit determination of near-Earth objects.

The method exploits a previously underexploited three-way identity in classical Keplerian mechanics:

```
K = G'/G  =  p/r  =  1 + e·cos(θ)
```

where `G' ≡ ω²r³/M` is the *instantaneous gravitational coupling* — not a new physical constant, but a reparametrisation of angular momentum conservation in gravitational units, expressed **locally at every orbital point**.

### Claimed performance (Apophis synthetic benchmark)

| Metric | Value |
|---|---|
| Candidate reduction | 95.6% (500 → 22) |
| True orbit rank | #1 (OCF = 0) |
| Observational arc | 30 days |
| Computation time | 38 s (500 candidates) |
| Scaling | O(N_c × N_obs) |
| Timeline compression | 3× vs traditional (30 vs 90 days) |

> **Scope.** This is a proof-of-concept implementation validating the method in the ideal two-body regime. Perturbation modelling, real observational data (MPC), and statistical benchmarking across asteroid populations are required for operational deployment — see §10 of the paper.

---

## Quick Start

```bash
git clone https://github.com/Maia-KarlosMarden/ocf-keplerian.git
cd ocf-keplerian
pip install -e ".[dev]"

# Reproduce Table 1 (Mercury verification)
python scripts/verify_mercury.py

# Reproduce Table 2 (Apophis benchmark)
python scripts/benchmark_apophis.py

# Arc-length sensitivity study
python scripts/arc_sensitivity.py

# Full test suite
pytest -v
```

---

## Installation

```bash
pip install -e .        # standard
pip install -e ".[dev]" # with pytest
pip install -e ".[notebooks]"  # with Jupyter + matplotlib
```

**Requirements:** Python ≥ 3.11, NumPy ≥ 1.26, SciPy ≥ 1.12.

---

## Repository Structure

```
ocf-keplerian/
├── ocf/
│   ├── __init__.py         # public API
│   ├── keplerian.py        # G', K identities, relational invariant
│   ├── propagator.py       # two-body Kepler propagator
│   ├── ocf.py              # OCF evaluation (Eq. 9)
│   ├── candidates.py       # ensemble generation
│   └── filters.py          # χ² threshold, filtering (Eq. 10)
├── tests/
│   ├── test_keplerian.py   # unit tests: foundational identities
│   └── test_ocf.py         # unit tests: OCF pipeline
├── scripts/
│   ├── verify_mercury.py   # Table 1 reproduction
│   ├── benchmark_apophis.py# Table 2 reproduction
│   └── arc_sensitivity.py  # arc-length sensitivity
├── docs/
│   └── theory.md           # derivation summary
├── pyproject.toml
└── README.md
```

---

## Core API

### Keplerian Index

```python
from ocf import KeplerianIndex
from ocf.keplerian import AU

ki = KeplerianIndex(e=0.1912, a=0.9224 * AU)

# All three forms of K — must agree
K_kin  = ki.at_theta(theta)                      # 1 + e·cos(θ)
K_geom = ki.p / ki.radius_at_theta(theta)        # p/r
K_dyn  = ki.G_prime_at_theta(theta) / 6.6743e-11 # G'/G

# Verify relational invariant G'·r = G·p
result = ki.verify_invariant(theta)
# {'G_prime_r': ..., 'G_p': ..., 'rel_error': ..., 'passes': True}
```

### Running the OCF on a custom asteroid

```python
import numpy as np
from ocf import compute_ocf_ensemble, generate_ensemble, filter_candidates
from ocf.filters import reduction_stats
from ocf.keplerian import AU

# Define your asteroid
a_true = 1.5 * AU     # semi-major axis [m]
e_true = 0.25         # eccentricity
i_true = np.radians(10.0)

# Observation epochs (e.g. 30 daily observations)
obs_times = np.linspace(0.0, 30 * 86400.0, 30)  # [s]

# Generate 500 candidate orbits
candidates, true_idx = generate_ensemble(
    a_true, e_true, i_true, n_candidates=500, seed=42
)

# Compute OCF for all candidates
results = compute_ocf_ensemble(
    candidates, obs_times,
    ref_propagator=candidates[true_idx],
    noise_sigma=0.0,   # 0 = ideal two-body; set >0 for noisy observations
)

# Filter
survivors, epsilon = filter_candidates(results, sigma=1e-6, n_obs=30)
stats = reduction_stats(results, survivors, true_id=true_idx)

print(f"Reduction: {stats['reduction_pct']:.1f}%")
print(f"True orbit rank: #{stats['true_rank']}")
```

### Using a custom propagator (BYOP)

You can supply any propagator that satisfies the `KeplerPropagator` interface
(i.e. exposes `.e`, `.p`, and `.state_at(t) → dict` with keys `K` and `theta`).
This allows integration with existing IOD pipelines.

---

## Theoretical Background

### The Relational Invariant

For any elliptical Keplerian orbit, the product `G'·r` is conserved:

```
G'(r) · r = G · p   [for all r along the orbit]
```

This is mathematically equivalent to angular momentum conservation (`h = const`),
but expressed as a **local, pointwise identity** rather than a global one.

### The Keplerian Index

```
K ≡ G'/G = p/r = 1 + e·cos(θ)
```

Range: `[1−e, 1+e]`, varying periodically with true anomaly.

At perihelion: `K = 1 + e`
At aphelion: `K = 1 − e`
At semi-latus rectum crossings: `K = 1`

### OCF Definition (Eq. 9)

For candidate orbit `O_k` and observations at epochs `t_i`:

```
OCF(O_k) = Σ_i [ p_k/r_{k,i} − (1 + e_k·cos θ_{k,i}) − δ_ref(t_i) ]²
```

- True orbit: OCF = 0 exactly
- False candidates: OCF > 0, growing with arc length
- `δ_ref` = common-mode subtraction from nominal Gauss orbit

### Chi-squared threshold (Eq. 10)

```
ε = σ² · χ²_α(N_obs)
```

where `σ` is the noise standard deviation on K and `α` is the confidence level (default 3σ).

---

## Tests

```bash
pytest -v                    # all tests
pytest tests/test_keplerian.py  # foundational identity tests
pytest tests/test_ocf.py        # pipeline tests
pytest --cov=ocf             # coverage report
```

The test suite verifies:
- `G'·r = G·p` to < 0.01% relative error at all orbital positions (Mercury, Apophis)
- Three-way identity `K_dynamics = K_geometry = K_kinematics`
- True orbit ranks #1 with OCF = 0
- Candidate reduction > 90% at 30 days
- True orbit survives the chi-squared filter

---

## Extending to Perturbations

Section 10.2 of the paper defines the perturbed OCF:

```
OCF_pert(O_k) = Σ_i [ ... − δ_ref(t_i) − δ_pert(t_i) ]²
```

To incorporate perturbations, subclass `KeplerPropagator` and override `state_at(t)`
to include planetary, SRP, Yarkovsky, and relativistic corrections.

---

## Contributing

Third-party validation of OCF against real MPC data is actively welcomed.
See the recommended development pathway in §11.2 of the paper.

Issues, pull requests, and benchmark results from real observational data
are all encouraged.

---

## Citation

```bibtex
@unpublished{maia2026ocf,
  author = {Maia, Karlos},
  title  = {The Orbital Coherence Function: A Keplerian Pre-Filter for Ensemble Orbit Determination},
  year   = {2026},
  note   = {Independent Researcher, Marietta, Georgia, USA}
}
```

---

## License

MIT License — see [LICENSE](LICENSE).
