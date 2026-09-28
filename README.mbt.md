# numsolve

Numerical solvers in pure MoonBit: root finding, ordinary differential
equations, linear systems, and interpolation — written against published
algorithms, and checked against a reference implementation of each.

```moonbit nocheck
///|
let root = @numsolve.brent(fn(x) { x * x - 2.0 }, 1.0, 2.0) // 1.4142135623731364

///|
let faster = @numsolve.newton(fn(x) { x * x - 2.0 }, fn(x) { 2.0 * x }, 1.5) // 1.4142135623730951, in four steps
```

## Root finding

`roots.mbt` finds where a function of one variable is zero. Three of the five
searches keep a bracket — two points whose function values have opposite
signs, which promises a root between them — and two are unbracketed, starting
from a point instead.

| Function | Method | Steps to 1e-10 on `x² − 2` from `[1, 2]` |
| --- | --- | --- |
| `bisection` | halve the interval, keep the half the sign change is in | 32 |
| `false_position` | the line through the two ends of the bracket, with the Illinois halving | 8 |
| `secant` | the line through the last two points, unbracketed | 6 |
| `newton` | the tangent at the current point, derivative needed | 4 |
| `brent` | bisection, the secant and inverse quadratic interpolation, chosen between | 6 |

Bisection is the slow one and the only one whose accuracy can be written down
before it runs. False position is bisection's bracket with a much better point
than the midpoint, and needs the Illinois modification to be usable at all: on
`x¹⁰ − 1`, which curves the same way across the whole bracket, the plain method
is two thousandths out after a hundred steps while the halving version is
exact in seventeen. The secant gets speed without a derivative, at the cost of
the bracket. Newton is the fastest and the most demanding. Brent is the one to
reach for by default: it keeps bisection's bracket, converges as fast as the
secant on smooth functions, and falls back to a bisection step whenever a
fitted one looks like a bad idea.

```moonbit nocheck
///|
let cubic = @numsolve.brent(fn(x) { x * x * x - x - 2.0 }, 1.0, 2.0)

///|
let sine = @numsolve.false_position(fn(x) { @math.sin(x) - 0.5 }, 0.0, 1.0)
```

The tolerance is measured in the variable being solved for, and each method
reads it in the way its own search can honour:

| Method | Stops when |
| --- | --- |
| `bisection`, `false_position`, `brent` | the interval left is no wider than `tolerance · (1 + \|x\|)` |
| `secant`, `newton` | the step is that small, or the function value is no larger than `tolerance` |

The `(1 + |x|)` factor is what makes one tolerance mean the same thing for a
root near zero and a root a long way from it. A tolerance of exactly zero
turns the tolerance test off, and the search then runs the full
`max_iterations` steps and returns what it has — that is how a fixed number of
steps is asked for, and what the convergence tests are built on. A function
value that comes out exactly zero ends a search at once, at any tolerance, and
an endpoint that is already a root is not searched for at all.

What a search cannot do, or cannot finish, is raised rather than answered
around:

| Variant | Payload | Raised when |
| --- | --- | --- |
| `NotBracketed(Double, Double)` | `f(a)` and `f(b)` | the two points given as a bracket have the same sign |
| `ZeroDerivative(Double)` | the point where it happened | a step would divide by a slope of exactly zero |
| `NonPositiveTolerance(Double)` | the tolerance | the tolerance is negative |
| `MaxIterations(Int)` | the limit | the search took `max_iterations` steps without meeting its tolerance |

```moonbit nocheck
///|
let root = @numsolve.bisection(fn(x) { x * x - 2.0 }, 1.0, 2.0) catch {
  _ => 0.0
}
```

A bracket may be given either way round, the two points of a bracket that hold
opposite signs promise a root and anything else is refused with both values,
and an unbracketed search — which promises nothing — reports a stall instead:
the secant on `x¹⁰ − 1` from 0 and 1.5 walks off to 1.9e11 and back, lands on a
point it repeats exactly, and reports the zero slope it is left with. Newton
on `x³ − 2x + 2` from zero steps to one and back to zero for ever, and reports
running out of iterations.

The tests measure the convergence each algorithm claims, not only its final
answer:

| Claim | Measured |
| --- | --- |
| Bisection after `k` steps is within `(b − a) / 2^k` of the root | held for k = 1, 2, 4, 8, 12, 16, 20 |
| Newton's error is about the square of the error before it | 8.6e-2, 2.5e-3, 2.1e-6, 1.6e-12, 0 |
| The secant's error falls as the golden ratio's power, 1.618 | log-ratios of 1.83, 1.68 and 1.67 for the last steps |
| The Illinois halving is what keeps false position off an end | exact in 17 steps where the plain method is 2.2e-3 out at 100 |

## Ordinary differential equations

`ode.mbt` walks a single equation `y'(t) = f(t, y)` from a starting value, and
returns a `Trajectory`: the times it stopped at and the value at each, the
first of them the starting point. A trajectory is read with `times()`,
`values()`, `end_value()` and `length()`. Integration backwards is allowed —
`t1` may be below `t0` — and the steps then have the sign of `t1 - t0`.

```moonbit nocheck
///|
let doubling = @numsolve.rk4(fn(_t, y) { y }, 1.0, 0.0, 1.0, 10) // 2.718279744135166
```

Four steppers, the first three of them told how many steps to take and the
fourth told how accurate to be:

| Stepper | Order | Evaluations per step | Error when the step halves |
| --- | --- | --- | --- |
| `euler` | 1 | 1 | ÷ 2 |
| `midpoint` | 2 | 2 | ÷ 4 |
| `rk4` | 4 | 4 | ÷ 16 |
| `rk45` | 5, adaptive | 7 | set by the tolerance |

On `y' = y` from 1 over `[0, 1]`, where the answer is e, ten steps of each
stepper give:

| Stepper | Ten steps | Error |
| --- | --- | --- |
| `euler` | 2.5937424601 | 1.2e-1 |
| `midpoint` | 2.714080846608224 | 4.2e-3 |
| `rk4` | 2.718279744135166 | 2.1e-6 |
| `rk45` at 1e-8 | 2.7182818362088534 | 7.8e-9, in 11 steps and 77 evaluations |

`euler` is here to be measured against rather than used. `rk4` is the default
for an interval that can be stepped evenly: order four from four evaluations
per step. `rk45` is Dormand-Prince 5(4)7M — seven evaluations produce both a
fifth-order answer and a fourth-order one, their difference estimates the
error of the step, and the step size of the next step is set from that
estimate. The trajectory it returns is sampled where it chose to stop, which
is neither even nor as many points: for the same error of about 8e-9 on
`y' = y` it spends 77 evaluations where `rk4` needs forty steps and 160.

The tolerance of `rk45` is not the tolerance of `roots.mbt`. It is the error a
step is allowed, measured relative to the value as `|y5 - y4| / (1 + |y5|)`,
and it must be positive — a tolerance of zero can never be met, and is refused
here rather than read as "take a fixed number of steps":

| Variant | Payload | Raised when |
| --- | --- | --- |
| `NonPositiveSteps(Int)` | the step count | a fixed-step stepper was given zero steps or fewer |
| `NonPositiveTolerance(Double)` | the tolerance | `rk45` was given a tolerance that is not positive |
| `MaxEvaluations(Int)` | the budget | `rk45` spent its evaluations of `f` without reaching `t1` |

`max_evaluations` is a budget on calls to `f` rather than on steps, and each
step costs seven of them; running out of it raises rather than returning a
trajectory that stops short of `t1`.

The tests measure the order each stepper claims — halving the step divides the
error by 2, 4 and 16, measured at 1.92, 3.85 and 15.35 — and that the adaptive
stepper meets the tolerance it is given rather than a fixed number of steps:
the same run at 1e-4, 1e-6, 1e-8 and 1e-10 lands at errors of 1.9e-5, 5.6e-7,
7.8e-9 and 8.7e-11, each under the tolerance asked for. A stiff equation,
`y' = −1000y`, is walked to the end of the interval by spending steps where
the solution moves, and a right-hand side that does not depend on `y` is
integrated exactly, because there a Runge-Kutta step is a quadrature rule.

## Linear systems

`linear.mbt` solves `A · x = b` for a dense square matrix, given as an
`Array[Array[Double]]` of rows.

```moonbit nocheck
///|
let a = [[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]

///|
let x = @numsolve.solve_linear_system(a, [1.0, 2.0, 4.0]) // [0.5, 0.0, 0.0]
```

| Function | Method | Use it when |
| --- | --- | --- |
| `LU::decompose` and `solve` | Gaussian elimination with partial pivoting, `P · A = L · U` | the matrix is general, and one factorization is reused for several right-hand sides |
| `Cholesky::decompose` and `solve` | `A = L · Lᵀ`, no pivoting | the matrix is symmetric positive definite: half the work of an LU |
| `solve_tridiagonal` | the Thomas recurrence over three diagonals | the matrix is tridiagonal — `O(n)` instead of `O(n³)` |

A decomposition is worth keeping when there is more than one right-hand side:
`LU::solve_many` answers with a list of them after one factorization, and
`LU::inverse` is that call with the identity. `LU::determinant` and
`Cholesky::determinant` come from the diagonal of the factorization at no
extra cost.

Pivoting is not an optimization here but the difference between an answer and
a wrong one. On `[[0, 2], [3, 4]]` there is no answer at all without a row
exchange — the pivot is zero — and on `[[1e-18, 1], [1, 1]]`, whose exact
solution is `[1, 1]`, elimination in the order the rows come in divides by
1e-18 and answers `[0, 1]`. Both cases are in the tests, the second one
against an elimination written without the pivot search.

What the arithmetic costs shows up in the condition number. The Hilbert
matrix of order 8 — `H[i][j] = 1 / (i + j + 1)`, condition number 1.5e10 —
solved for the right-hand side made of its row sums has the exact solution of
all ones, and comes back with a largest error of 1.4e-7: about ten of the
sixteen digits a double has are left, and that is the matrix rather than the
method. The test asserts 1e-5 rather than that figure, so that a different
order of summation inside the same algorithm cannot fail it.

| Variant | Payloads | Raised when |
| --- | --- | --- |
| `NotSquare(Int, Int)` | rows, columns | the matrix is not square, or its rows are not all the same length |
| `SizeMismatch(Int, Int)` | needed, given | a right-hand side or a diagonal is the wrong length |
| `Singular(Int)` | the column | a pivot is exactly zero and no row below it can be exchanged for it |
| `NotPositiveDefinite(Int)` | the index | Cholesky was given a matrix whose pivot is not positive |

## Interpolation

`interp.mbt` draws a curve through a table of points — two arrays, `xs` and
`ys`, of the same length, with the `xs` in increasing order — and reads the
value it takes between them.

```moonbit nocheck
///|
let table_x = [0.0, 1.0, 2.0, 3.0]

///|
let table_y = [1.0, 2.0, 0.0, 4.0]

///|
let line = @numsolve.LinearInterp::new(table_x, table_y)

///|
let spline = @numsolve.CubicSpline::new(
  table_x,
  table_y,
  @numsolve.SplineBoundary::Natural,
)
```

| Type | The curve | Use it when |
| --- | --- | --- |
| `LinearInterp` | the straight line between the two knots around `x` | the table is dense, or a kink at every knot does not matter |
| `CubicSpline` | one cubic per interval, joined so that value, slope and second derivative agree at every knot | the curve has to be smooth: it is `C²`, and its slope is continuous across the knots |

A spline is the piecewise cubic whose second derivatives `M[i]` at the knots
solve `h[i-1]·M[i-1] + 2(h[i-1] + h[i])·M[i] + h[i]·M[i+1] = 6·(s[i] - s[i-1])`,
where `h` is an interval's width and `s` its slope — `n` unknowns and `n - 2`
equations, so two more conditions are needed and `SplineBoundary` is which
ones:

| Boundary | The two extra conditions | End behaviour |
| --- | --- | --- |
| `Natural` | `M[0] = M[n-1] = 0` | free to bend at both ends: the choice when nothing is known about them |
| `Clamped(d0, dn)` | `S'(x0) = d0` and `S'(xn) = dn` | prescribed slopes at the two ends |
| `NotAKnot` | the third derivative is continuous at the second and second-to-last knot | no assumption about the ends: the two outermost cubics are one polynomial |

The system is solved by `solve_tridiagonal` from `linear.mbt`, so a spline
through `n` points costs `O(n)` — one pass down and one back — rather than the
`O(n³)` of a dense elimination. The not-a-knot condition is the one that does
not fit that system as a row: it relates three unknowns where row 0 has room
for two, so the end curvature is eliminated with it rather than solved for, and
read back once the interior ones are known.

On the table above, the three boundaries give:

| Boundary | `S(0.5)` | `S(1.5)` | `S(2.5)` |
| --- | --- | --- | --- |
| `Natural` | 1.95 | 0.775 | 1.325 |
| `Clamped(0, 0)` | 1.65 | 0.625 | 2.225 |
| `NotAKnot` | 2.4375 | 0.8125 | 0.6875 |

`moments()` is the second derivatives `M` the cubics are built from — `[0,
-7.2, 10.8, 0]` for the natural spline above — and `coefficients()` is the
cubics themselves, one row `[a, b, c, d]` per interval for the polynomial
`a·x³ + b·x² + c·x + d` in the global `x` rather than in the distance from the
knot at the left of the interval. Two neighbouring rows of it describe the
same value, the same slope and the same second derivative at the knot they
share, which is what the tests assert to check the `C²` claim at every knot
under all three boundaries.

A table of two points is one interval, and under `Natural` or `Clamped` that
is a straight line; a not-a-knot spline needs a knot that is not one, so it
needs three.

A point outside the range of the table has no answer — a spline describes the
interval it was built from, not the line beyond it — so `evaluate` returns
`None` there, and what cannot be built at all is raised:

| Variant | Payload | Raised when |
| --- | --- | --- |
| `TooFewKnots(Int)` | the number of points | fewer than two points were given, or fewer than three for a not-a-knot spline |
| `UnevenLengths(Int, Int)` | the two lengths | `xs` and `ys` are not the same length |
| `NotSorted(Double)` | the first knot out of order | the knots do not increase |
| `DuplicateKnots(Double)` | the repeated knot | two knots are the same point |

```moonbit nocheck
///|
let value = @numsolve.CubicSpline::new(table_x, table_y, @numsolve.SplineBoundary::Natural) catch {
  _ => abort("the table is a valid one")
}.evaluate(1.5) // Some(0.775)
```

The tests hold the fit to what SciPy's `CubicSpline` returns for the same
table under the same boundary, and to the properties a spline is chosen for:
every knot is reproduced exactly, the three boundaries agree on a table that
lies on a straight line, and the cubics of neighbouring intervals meet in
value, slope and second derivative. Sampling `e^(x)` at 0, 0.25, 0.5 and 1 and
fitting with its true end slopes 1 and e puts the curve within 2.2e-5 of it at
0.125 and within 4.4e-4 at 0.75 — the price of a cubic per interval.

## Command line

The solvers are also wired to a command-line demo under `cmd/main`, which runs
one solver on one problem and prints its answer next to the number the problem
is known to have:

```
moon run cmd/main -- root sqrt2
moon run cmd/main -- ode growth --steps 100
moon run cmd/main -- spline natural
moon run cmd/main -- solve hilbert
```

| Command | Subjects | Solver | Prints |
| --- | --- | --- | --- |
| `root` | `sqrt2`, `cubic`, `trig` | `brent`, `false_position` | the root, against √2 or π/6 |
| `ode` | `growth`, `decay`, optionally `--steps N` | `rk45`, or `rk4` with N steps | the end value, e or 1/e, and the error |
| `spline` | `natural`, `clamped`, `notaknot` | `CubicSpline` through (0,1), (1,2), (2,0), (3,4) | every piece's coefficients and three samples |
| `solve` | `demo`, `hilbert` | `solve_linear_system` | the 3×3 solution, or how close Hilbert 8's is to all ones |

`--help` prints the same list, and anything that is not a command is refused
with a message rather than answered with a guess — the parser is tested on
both. The demo lives in a package of its own, so its tests also hold the
package boundary to its word: they call `brent` and the other solvers from
outside, with and without the optional arguments, and assert the same values
the library's own tests do.

## Testing

```
moon test
```

66 tests, and every golden value in them was computed first with SciPy 1.18
and NumPy 2.5 and then written down in full: SciPy's `brentq` on `x² − 2`
returns 1.4142135623731364, which is what the test holds Brent to, and the
convergence rates above are read off the iterates rather than assumed. The
ODE goldens come from a model of the same arithmetic — the same step sizes,
the same order of operations — so the adaptive stepper's eleven steps and 77
evaluations are what the test asserts, not a range it is allowed to fall in.

The package also carries benchmarks, run against the same functions:

```
moon bench --release
```
