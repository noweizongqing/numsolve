# ODE design and scope

## Why this layer exists

`numsolve` is a pure MoonBit numerical toolkit. The ODE layer follows the
usual explicit Runge-Kutta formulas so that the implementation is easy to
audit and can run in MoonBit's Wasm target without a Python, C, or Julia
runtime.

That makes the contribution an integration and API contribution rather than
a new numerical method. Mature ecosystems such as SciPy (`solve_ivp`),
Sundials, and DifferentialEquations.jl cover a much larger design space,
including implicit and stiff methods, event detection, dense output, and
sensitivity analysis. This package does not claim to replace them. It aims to
provide a coherent subset for MoonBit programs and to make its boundaries
explicit.

## Public API

The original scalar functions remain available:

| Function | Method | State |
| --- | --- | --- |
| `euler` | explicit Euler | `Double` |
| `midpoint` | explicit midpoint | `Double` |
| `rk4` | classic fourth-order Runge-Kutta | `Double` |
| `rk45` | Dormand-Prince 5(4) with step-size control | `Double` |

The system layer uses `Array[Double]` for the state and adds:

| Function | Method | State |
| --- | --- | --- |
| `euler_system` | explicit Euler | `Array[Double]` |
| `midpoint_system` | explicit midpoint | `Array[Double]` |
| `rk4_system` | classic fourth-order Runge-Kutta | `Array[Double]` |
| `rk45_system` | Dormand-Prince 5(4) with component tolerances | `Array[Double]` |
| `solve_ivp` | dispatch through `OdeMethod` | `Array[Double]` |

`VectorTrajectory` exposes `times`, `states`, `end_state`, `length`, and
`dimension`. The state returned by `solve_ivp` therefore has a shape familiar
to users of array-oriented ODE libraries.

## Numerical design

For a problem

```text
y'(t) = f(t, y),    y(t0) = y0,
```

the fixed-step methods divide `[t0, t1]` into equal steps. The adaptive
Dormand-Prince method keeps its own accepted steps. It computes a fifth-order
and a fourth-order candidate from the same seven stages, accepts the step when
the estimated error is within tolerance, and sizes the next step from that
estimate.

The scalar `rk45` uses

```text
|y5 - y4| / (1 + |y5|)
```

as its relative error measure. `rk45_system` scales each component by both an
absolute and a relative tolerance:

```text
scale_i = atol + rtol * max(|y_old_i|, |y_new_i|)
error   = sqrt(mean(((y5_i - y4_i) / scale_i)^2))
```

The step is accepted when `error <= 1`. This prevents a component near zero
from being ignored or allowed to dominate a component with a much larger
magnitude.

Dimension mismatches, empty system states, non-positive tolerances, invalid
step counts, and exhausted evaluation budgets are reported as `OdeError`
values instead of being converted into a plausible but wrong trajectory.

## Verification

The tests include:

- order checks for Euler, midpoint, and RK4 on a known exponential system;
- a two-dimensional harmonic oscillator over a full period;
- scalar and vector adaptive runs against analytic solutions;
- backward integration;
- `solve_ivp` dispatch;
- dimension and empty-state validation.

The core scalar ODE tests also retain their numerical goldens. These checks
measure the advertised convergence behaviour and the actual error control,
not only a final number copied from a reference.

## Not yet covered

The following are deliberate boundaries of this version:

- no event detection or root-of-solution callbacks;
- no dense output or continuous interpolation between accepted steps;
- no implicit, Rosenbrock, or BDF methods for stiff systems;
- no automatic Jacobian, sensitivities, or complex states;
- no sparse matrix or large-scale state support.

## Roadmap

The next useful additions, in dependency order, would be:

1. optional evaluation at user-provided times using dense output;
2. event detection with terminal and directional options;
3. implicit solvers for stiff systems, sharing a linear-system interface;
4. multiple trajectories or batched systems for common parameter sweeps;
5. stronger cross-validation against SciPy on a broader problem set.

Each addition should keep the existing scalar APIs stable and continue to
return explicit errors for unsupported or ill-posed inputs.
