# SNOPE

**S**pherical-symmetric spacetime **N**umerical **O**rbit and **P**arameter **E**stimation.

One integrator, one MCMC engine, one plotting module, shared across every
spacetime you want to test against the S2 star's orbit around Sgr A*.
Schwarzschild-de Sitter, quadratic gravity, and Brans-Dicke ship as
worked examples. Adding a new metric means writing one short file, not
another 800-line script.

## Why this exists

The three original scripts (SdS, quadratic gravity, Brans-Dicke) each
carried their own copy of the same ~800 lines: data loading, geodesic
integration, Roemer delay, relativistic Doppler + gravitational
redshift, sky projection, MCMC bookkeeping, corner/convergence/orbit
plotting. Only the metric functions and priors actually differed
between them. SNOPE pulls all of that shared machinery into one tested
core (`snope/`) and reduces each metric to a single config file
(`snope/metrics/*.py`) defining g_tt(r), g_rr(r), the new physics
parameter(s), and their priors.

Nothing about accuracy was sacrificed to get here: every metric still
integrates with `DOP853` at `rtol=atol=1e-12`, at whatever `n_points`
you ask for -- exactly as the original scripts did.

## Install

```bash
git clone <this repo>
cd SNOPE
pip install -r requirements.txt
```

Needs Python 3.9+. The MCMC step uses `emcee`'s multiprocessing `Pool`,
so wrap standalone scripts in `if __name__ == "__main__":` on platforms
where that matters. Notebooks are fine as they are.

## Quickstart

Drop your S2 astrometry/RV tables in `data/` (format described in
`data/README.md`), then:

```python
from snope.metrics.sds import main

sampler, flat_samples, best_fit = main(
    data_path_pos="data/tab_gillessen_pos.csv",
    data_path_rv="data/tab_gillessen_vr.csv",
)
```

That's the whole notebook. It optimizes a good starting point with
L-BFGS-B, runs the `emcee` MCMC, prints a chi2/AIC/BIC summary, and
saves a corner plot, a walker-trace plot, and a best-fit orbit/RV plot.
Swap `sds` for `quadratic_gravity` or `brans_dicke` to fit a different
spacetime -- nothing else changes.

`main()` takes overrides for anything in the pipeline:

```python
from snope.priors import PriorMode

sampler, flat_samples, best_fit = main(
    prior_mode=PriorMode.FLAT,      # see "Prior modes" below
    n_walkers=64, n_steps=50000, burn_in=5000, n_points=1000,
    n_cores=8,
)
```

See `notebooks/` for one ready-to-run notebook per metric, including a
quick effective-potential check (see below) before you commit to a
full MCMC run.

## Layout

```
snope/
  constants.py        physical constants (G, c, yr->s, arcsec->rad)
  data.py               S2 astrometry/RV CSV loader
  orbit_model.py        S2OrbitModel -- the shared physics engine
                         (geodesic integration, Roemer delay, redshift,
                         sky projection). Metric-agnostic.
  priors.py             ParamPrior / PriorSet -- flat / mixed / Gaussian
  mcmc.py                S2OrbitMCMC -- the emcee wrapper
  plotting.py            chi2/AIC/BIC summary, corner, convergence, orbit plots
  pipeline.py             run_fit() / make_main() -- wires it all together
  metrics/
    __init__.py           MetricSpec interface (start here to add a metric)
    sds.py                 Schwarzschild-de Sitter    (parameter: log10_lambda)
    quadratic_gravity.py   Quadratic gravity            (parameter: k)
    brans_dicke.py          Brans-Dicke (PPN form)       (parameter: omega_bd)
notebooks/               one notebook per metric
data/                     put your S2 tables here
tests/                    fast sanity tests, no MCMC (tests/test_metrics.py)
```

## The physics engine, briefly

Every metric here is static, diagonal, and equatorial:
`ds^2 = g_tt(r) dt^2 + g_rr(r) dr^2 + r^2 dphi^2`, with `r` in units of
`R_s = 2GM_bh/c^2` (so the horizon sits near `r=1`). `S2OrbitModel`
integrates the geodesic equations using the general diagonal-metric
Christoffel symbols:

```
Gamma^t_tr     =  g_tt'/(2 g_tt)
Gamma^r_tt     = -g_tt'/(2 g_rr)
Gamma^r_rr     =  g_rr'/(2 g_rr)
Gamma^r_phiphi = -r/g_rr
Gamma^phi_rphi =  1/r
```

This is a standard identity for any diagonal metric of this form -- it
doesn't assume `g_rr = -1/g_tt`, so it's correct whether your metric is
Schwarzschild-like (`g_rr = 1/(-g_tt)`, e.g. SdS) or has independent
`g_tt` and `g_rr` (quadratic gravity, Brans-Dicke). One integrator,
correct for anything you plug in.

The turning-point condition used to find the orbit's conserved energy
and angular momentum (`find_EL`), and the combined gravitational +
Doppler redshift used in `compute_observables`, both depend on `g_tt`
alone -- so they're identical across metrics too, and also live in the
shared core rather than in each metric file.

## Checking a metric before you fit it: the effective potential

Before spending hours on a full MCMC run, it's worth checking that your
chosen metric + parameter values actually produce a bound orbit --
i.e. that the effective potential `V_eff(r) = -g_tt(r) * (1 + L^2/r^2)`
has the periapsis/apoapsis turning points you expect at `E^2`. Each
notebook has a cell for this: it plots `V_eff(r)` across the orbit's
radial range, marks `E^2`, and shows where the two curves cross. If
they don't cross twice, the orbit isn't bound for those parameters, and
the integrator will fail or return garbage -- better to catch that here
than after a long MCMC run.

## Adding a new metric

Create `snope/metrics/my_metric.py`:

```python
from . import MetricSpec
from ..priors import ParamPrior
from ..pipeline import make_main
from ..priors import PriorMode

M_VAL = 0.5  # Schwarzschild-convention geometric mass (r in units of R_s)

def gtt(r, R_s, my_param):
    ...   # return g_tt(r)

def grr(r, R_s, my_param):
    ...

def dgtt_dr(r, R_s, my_param):
    ...

def dgrr_dr(r, R_s, my_param):
    ...

metric = MetricSpec(
    name="my_metric",
    extra_param_names=("my_param",),
    gtt=gtt, grr=grr, dgtt_dr=dgtt_dr, dgrr_dr=dgrr_dr,
    horizon_radius=lambda R_s, **extra: 1.0,
    keplerian_mass_guess=lambda R_s, **extra: M_VAL,
    description="One-line description of the metric",
)

param_priors = {
    'M_bh':     ParamPrior(bounds=(3.7e6, 4.9e6), gaussian=(4.3e6, 0.2e6)),
    'distance': ParamPrior(bounds=(5.33, 11.33),  gaussian=(8.33, 0.5)),
    'a':        ParamPrior(bounds=(112.0, 139.0), gaussian=(125.5, 2.0)),
    'e':        ParamPrior(bounds=(0.855, 0.912), gaussian=(0.884, 0.005)),
    'i':        ParamPrior(bounds=(128.17, 140.20), gaussian=(134.18, 1.5)),
    'omega':    ParamPrior(bounds=(56.95, 74.03), gaussian=(65.49, 2.0)),
    'Omega':    ParamPrior(bounds=(217.95, 235.94), gaussian=(226.95, 2.0)),
    't_peri':   ParamPrior(bounds=(2002.07, 2002.57), gaussian=(2002.32, 0.05)),
    'x0':       ParamPrior(bounds=(-0.001, 0.001), gaussian=(0.0, 0.0002)),
    'y0':       ParamPrior(bounds=(-0.001, 0.001), gaussian=(0.0, 0.0002)),
    'vx0':      ParamPrior(bounds=(-0.0005, 0.0005), gaussian=(0.0, 0.0001)),
    'vy0':      ParamPrior(bounds=(-0.0005, 0.0005), gaussian=(0.0, 0.0001)),
    'vz0':      ParamPrior(bounds=(-25.0, 25.0), gaussian=(0.0, 5.0)),
    'my_param': ParamPrior(bounds=(lo, hi), gaussian=(mean, std)),
}

main = make_main(metric, param_priors, default_prior_mode=PriorMode.MIXED,
                  n_walkers=40, n_steps=30000, burn_in=3000, n_points=1000)
```

That's it -- `from snope.metrics.my_metric import main` now works just
like the three built-in metrics. `R_s` (metres) is passed to every
metric function in case your parameter has physical units that need
converting into the dimensionless `r` used internally (see `sds.py`'s
`log10_lambda`, in m^-2); ignore it otherwise (see
`quadratic_gravity.py`, `brans_dicke.py`).

If your metric needs a parameter-dependent horizon location (like
quadratic gravity's `rh(k)`), implement that in `horizon_radius` -- it
sets the integration-clamp radius and, by default, seeds the Keplerian
guess for the root-solver.

## Prior modes

Every parameter in a metric's `param_priors` dict carries both a flat
range and a Gaussian `(mean, std)` at once. `prior_mode` (an argument to
`main()`, or to `PriorSet` directly) picks which one is actually used,
for every parameter, with a single switch:

| mode | orbital elements + new-physics parameter | t_peri, x0, y0, vx0, vy0, vz0 |
|---|:---:|:---:|
| `PriorMode.FLAT`             | flat | flat |
| `PriorMode.MIXED` (default)  | flat | Gaussian |
| `PriorMode.GAUSSIAN`         | Gaussian | Gaussian |

```python
from snope.priors import PriorMode
main(prior_mode=PriorMode.FLAT)
main(prior_mode=PriorMode.MIXED)      # default
main(prior_mode=PriorMode.GAUSSIAN)
```

To change which parameters count as "offsets" in mixed mode, pass
`gaussian_params=(...)` to `main()`. A Gaussian parameter's hard sampler
bounds are always mean +/- 5*std (`priors.N_SIGMA_BOUND`), so switching
modes never changes where the sampler is allowed to go -- only whether
it's penalized for straying from the center.

## Outputs

Each run writes files prefixed with the metric's name (e.g. `sds_*`):

```
<prefix>_flat_samples.npy       flattened post-burn-in chain
<prefix>_samples.npy             full chain (n_steps, n_walkers, n_params)
<prefix>_log_likelihoods.npy     post-burn-in log-likelihoods
<prefix>_statistical_info.txt    chi2/AIC/BIC + best-fit params + autocorr times
<prefix>_corner.png
<prefix>_convergence.png
<prefix>_orbit.png               4-panel: sky orbit, RV(t), RA(t), Dec(t)
```

## A note on runtime

The defaults ported from the original scripts (`n_steps=30000,
n_walkers=40, n_points=1000`) are a cluster-scale run: every likelihood
evaluation re-integrates the geodesic at `rtol=atol=1e-12`, and
L-BFGS-B's numerical gradient means the initial optimization alone can
take on the order of `100 x n_params` evaluations. For local iteration,
turn everything down first -- e.g. `n_points=200, n_steps=200` -- to
sanity-check a new metric or prior choice quickly, then scale back up
for the real run. `n_cores` defaults to all available cores via
`multiprocessing.Pool`.

## Tests

```bash
pip install pytest
pytest tests/
```

`tests/test_metrics.py` builds small synthetic S2-like data on the fly
and checks that every built-in metric integrates a closed orbit and
returns finite, physically plausible observables. It's a fast smoke
test, not a scientific validation of any particular metric's physics.
