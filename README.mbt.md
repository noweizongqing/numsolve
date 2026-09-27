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

## Testing

```
moon test
```

18 tests, and every golden value in them was computed first with SciPy 1.18
and NumPy 2.5 and then written down in full: SciPy's `brentq` on `x² − 2`
returns 1.4142135623731364, which is what the test holds Brent to, and the
convergence rates above are read off the iterates rather than assumed.

The package also carries benchmarks, run against the same functions:

```
moon bench --release
```
