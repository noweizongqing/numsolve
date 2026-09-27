///|
/// Benchmarks for the searches a root is found with, run against the
/// quadratic x² − 2 and the cubic x³ − x − 2:
///
///     moon bench --release benchmarks
///
/// Each benchmark is named for one search, so a number reads as what that
/// search costs on a smooth function — and so that a change to bisection's
/// halving, Brent's curve fitting or Newton's tangent can be seen in the
/// numbers rather than argued about.
///
/// The timed closures do not raise: a bracket the tests already assert is
/// assumed to be good, and a failure inside a timed closure would only
/// measure how fast the benchmark fails.

///|
/// The function the slower benchmarks solve, and its derivative.
fn quadratic(x : Double) -> Double {
  x * x - 2.0
}

///|
/// The derivative of the quadratic, which Newton's method needs.
fn quadratic_derivative(x : Double) -> Double {
  2.0 * x
}

///|
/// The cubic x³ − x − 2, whose root is 1.5213797068045676 — a function that is
/// not symmetric about its root, which is what the curve-fitting searches are
/// for.
fn cubic(x : Double) -> Double {
  x * x * x - x - 2.0
}

///|
/// A root found by halving an interval: the slow search, and the one whose
/// cost is known in advance — thirty-two steps for a tolerance of 1e-10,
/// whatever the function does.
test "a bisection root" (it : @bench.T) {
  let fallback = @numsolve.bisection(quadratic, 1.0, 2.0) catch {
    _ => abort("x^2 - 2 is bracketed by 1 and 2")
  }
  it.bench(fn() {
    it.keep(@numsolve.bisection(quadratic, 1.0, 2.0) catch { _ => fallback })
  })
}

///|
/// A root found by false position with the Illinois halving: the same
/// function in a quarter of bisection's steps, each of them a line drawn
/// across the remaining interval.
test "a false position root" (it : @bench.T) {
  let fallback = @numsolve.false_position(quadratic, 1.0, 2.0) catch {
    _ => abort("x^2 - 2 is bracketed by 1 and 2")
  }
  it.bench(fn() {
    it.keep(
      @numsolve.false_position(quadratic, 1.0, 2.0) catch {
        _ => fallback
      },
    )
  })
}

///|
/// A root found by Brent's method: a few steps, each of them a curve fitted
/// through the last three points and checked against a bisection step before
/// it is taken.
test "a brent root" (it : @bench.T) {
  let fallback = @numsolve.brent(cubic, 1.0, 2.0) catch {
    _ => abort("x^3 - x - 2 is bracketed by 1 and 2")
  }
  it.bench(fn() {
    it.keep(@numsolve.brent(cubic, 1.0, 2.0) catch { _ => fallback })
  })
}

///|
/// A root found by Newton's method: four steps, each of them one function
/// value and one derivative, which is the cheapest search here when the
/// derivative is at hand.
test "a newton root" (it : @bench.T) {
  let fallback = @numsolve.newton(quadratic, quadratic_derivative, 1.5) catch {
    _ => abort("x^2 - 2 has a tangent at 1.5")
  }
  it.bench(fn() {
    it.keep(
      @numsolve.newton(quadratic, quadratic_derivative, 1.5) catch {
        _ => fallback
      },
    )
  })
}
