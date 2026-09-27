///|
/// Root finding: where a function of one variable is zero.
///
/// Five searches, each taking the function and either a bracket or a starting
/// point, and each asked for a `tolerance` and a `max_iterations`.
///
/// The three bracketing methods — `bisection`, `false_position` and `brent` —
/// are given two points whose function values have opposite signs, which is a
/// promise that a root lies between them. The two open methods — `secant` and
/// `newton` — are given starting points and promise nothing: they settle on
/// whichever root their start leads them to, and can stop short or run away
/// when the start is poor.
///
/// The tolerance is measured in the variable being solved for, and each
/// method reads it in the way its own search can honour:
///
/// | Method | Stops when |
/// | --- | --- |
/// | `bisection`, `false_position`, `brent` | the interval left is no wider than `tolerance · (1 + \|x\|)` |
/// | `secant`, `newton` | the step is that small, or the function value is no larger than `tolerance` |
///
/// The `(1 + |x|)` factor is what makes one tolerance mean the same thing for
/// a root near zero and a root a long way from it: the answer comes out right
/// in the digits that matter rather than to an absolute width.
///
/// A `tolerance` of exactly zero turns the tolerance test off, and the search
/// then runs the full `max_iterations` steps and returns what it has: that is
/// how a fixed number of steps is asked for. A negative tolerance is refused.
/// A function value that comes out exactly zero ends a search at once, at any
/// tolerance, and an endpoint that is already a root is not searched for at
/// all.

///|
/// A root search that could not be run, or that ran out of iterations.
pub suberror RootError {
  /// The two points offered as a bracket have function values of the same
  /// sign, so no root is promised between them. The payloads are those two
  /// values, `f(a)` and `f(b)`.
  NotBracketed(Double, Double)
  /// A step would have had to divide by a slope that came out exactly zero —
  /// the derivative in Newton's method, or the difference between the last
  /// two function values in the secant method. The payload is the point at
  /// which the search stopped.
  ZeroDerivative(Double)
  /// The tolerance is negative, and no search can be run with one. The
  /// payload is the tolerance. Zero is not refused: it asks for a fixed
  /// number of iterations.
  NonPositiveTolerance(Double)
  /// The search took `max_iterations` steps without meeting its tolerance.
  /// The payload is that limit.
  MaxIterations(Int)
}

///|
/// The gap between 1.0 and the double above it, which is the smallest
/// relative step a double can still tell apart from no step at all. Brent's
/// method uses it to keep its steps meaningful once the tolerance is down at
/// the rounding of the numbers themselves.
let double_epsilon : Double = 2.220446049250313e-16

///|
/// Bisection: halve the interval and keep the half the sign change is in.
///
/// It is the slowest of the five — one bit of the answer per step, three or
/// four steps per decimal digit — and the only one whose accuracy is known
/// before it runs: whatever the function does, after `k` steps the root is
/// within `(b - a) / 2^k` of the point returned, and the half-width of the
/// interval is what the tolerance is compared against.
///
/// The two points may be given in either order.
fn bisection_run(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance : Double,
  max_iterations : Int,
) -> (Double, Int) raise RootError {
  if tolerance < 0.0 {
    raise RootError::NonPositiveTolerance(tolerance)
  }
  let f_a = f(a)
  if f_a == 0.0 {
    return (a, 0)
  }
  let f_b = f(b)
  if f_b == 0.0 {
    return (b, 0)
  }
  if f_a * f_b > 0.0 {
    raise RootError::NotBracketed(f_a, f_b)
  }
  let mut left = a
  let mut right = b
  let mut f_left = f_a
  let mut middle = 0.5 * (a + b)
  let mut iterations = 0
  let mut converged = false
  while iterations < max_iterations {
    middle = 0.5 * (left + right)
    let f_middle = f(middle)
    iterations = iterations + 1
    if f_middle == 0.0 {
      converged = true
      break
    }
    if f_left * f_middle < 0.0 {
      right = middle
    } else {
      left = middle
      f_left = f_middle
    }
    if tolerance > 0.0 &&
      (right - left).abs() <= tolerance * (1.0 + middle.abs()) {
      converged = true
      break
    }
  }
  if tolerance > 0.0 && !converged {
    raise RootError::MaxIterations(max_iterations)
  }
  (middle, iterations)
}

///|
/// The root of `f` between `a` and `b`, by bisection.
///
/// The two points have to bracket a root — their function values of opposite
/// signs — and either of them may be the root itself. The answer is the
/// midpoint of the last interval, so it is within `tolerance · (1 + |x|)` of
/// the root, and within `(b - a) / 2^k` after `k` steps.
///
/// It is the method to reach for when nothing is known about the function
/// beyond its sign at two points, and it never fails on a promise that holds.
pub fn bisection(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 100,
) -> Double raise RootError {
  bisection_run(f, a, b, tolerance, max_iterations).0
}

///|
/// False position: take the point where the straight line through the two
/// ends of the bracket crosses zero, and keep the half the sign change is in,
/// as bisection does but with a much better point than the midpoint.
///
/// The line alone is not enough on a function that curves the same way across
/// the whole interval: it holds one end of the bracket in place and each step
/// then buys almost nothing. This is the Illinois modification of the method,
/// which is what gives the search its speed back — when an end is kept two
/// steps in a row, the function value it is kept for is halved, and the next
/// line is drawn across the interval at a new angle.
///
/// The two points may be given in either order.
fn false_position_run(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance : Double,
  max_iterations : Int,
) -> (Double, Int) raise RootError {
  if tolerance < 0.0 {
    raise RootError::NonPositiveTolerance(tolerance)
  }
  let f_a = f(a)
  if f_a == 0.0 {
    return (a, 0)
  }
  let f_b = f(b)
  if f_b == 0.0 {
    return (b, 0)
  }
  if f_a * f_b > 0.0 {
    raise RootError::NotBracketed(f_a, f_b)
  }
  let mut left = a
  let mut right = b
  let mut f_left = f_a
  let mut f_right = f_b
  // Which end the last step kept: 0 for neither, 1 for the left, 2 for the
  // right. An end that is kept twice is the one the Illinois halving is
  // applied to.
  let mut kept = 0
  let mut x = a
  let mut iterations = 0
  let mut converged = false
  while iterations < max_iterations {
    x = (left * f_right - right * f_left) / (f_right - f_left)
    let f_x = f(x)
    iterations = iterations + 1
    if f_x == 0.0 {
      converged = true
      break
    }
    if f_left * f_x < 0.0 {
      right = x
      f_right = f_x
      if kept == 1 {
        f_left = 0.5 * f_left
      }
      kept = 1
    } else {
      left = x
      f_left = f_x
      if kept == 2 {
        f_right = 0.5 * f_right
      }
      kept = 2
    }
    if tolerance > 0.0 && (right - left).abs() <= tolerance * (1.0 + x.abs()) {
      converged = true
      break
    }
  }
  if tolerance > 0.0 && !converged {
    raise RootError::MaxIterations(max_iterations)
  }
  (x, iterations)
}

///|
/// The root of `f` between `a` and `b`, by false position with the Illinois
/// modification.
///
/// The bracket is read the same way as `bisection`'s, and on most functions
/// the answer arrives in a fraction of the steps: six to nine rather than
/// thirty-odd for a tolerance of 1e-10. On a function that curves one way
/// across the whole bracket the modification is what keeps it from stalling —
/// the tests hold `x^10 - 1` as the case in point, where the same search
/// without the halving is still two thousandths out after a hundred steps.
pub fn false_position(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 100,
) -> Double raise RootError {
  false_position_run(f, a, b, tolerance, max_iterations).0
}

///|
/// The secant method: draw the line through the last two points and see where
/// it crosses zero, over and over. The derivative is never asked for, and
/// neither is a bracket: the two points are just a start.
///
/// Being unbracketed, it can leave the interval it started in and settle on
/// another root, or on nothing at all — a poor start is answered with a poor
/// root rather than with a failure, so the result is worth reading against the
/// function. What it buys is speed: its error falls as the golden ratio's
/// power, about 1.6, against bisection's one bit per step.
fn secant_run(
  f : (Double) -> Double,
  x0 : Double,
  x1 : Double,
  tolerance : Double,
  max_iterations : Int,
) -> (Double, Int) raise RootError {
  if tolerance < 0.0 {
    raise RootError::NonPositiveTolerance(tolerance)
  }
  let mut x_previous = x0
  let mut f_previous = f(x0)
  if f_previous == 0.0 {
    return (x0, 0)
  }
  let mut x = x1
  let mut f_x = f(x1)
  if f_x == 0.0 {
    return (x1, 0)
  }
  let mut iterations = 0
  let mut converged = false
  while iterations < max_iterations {
    if f_x - f_previous == 0.0 {
      raise RootError::ZeroDerivative(x)
    }
    let next = x - f_x * (x - x_previous) / (f_x - f_previous)
    iterations = iterations + 1
    x_previous = x
    f_previous = f_x
    x = next
    f_x = f(next)
    if f_x == 0.0 {
      converged = true
      break
    }
    let moved = x - x_previous
    // A point that repeats exactly is a point the iteration has stuck on, and
    // being stuck is not converging unless the function is zero there: the
    // secant can land on a fixed point of its own that is not a root, and
    // reporting that as an answer would be reporting a point the function is
    // nowhere near zero at. Such a search keeps going and ends by running out
    // of iterations, which is what happened.
    if tolerance > 0.0 &&
      (
        f_x.abs() <= tolerance ||
        (moved != 0.0 && moved.abs() <= tolerance * (1.0 + x.abs()))
      ) {
      converged = true
      break
    }
  }
  if tolerance > 0.0 && !converged {
    raise RootError::MaxIterations(max_iterations)
  }
  (x, iterations)
}

///|
/// A root of `f`, by the secant method, from the two starting points `x0` and
/// `x1`.
///
/// The two points say which part of the function to work in; a root is not
/// promised between them, and a start far from a root can return a point that
/// is not one. Either starting point may be a root itself.
pub fn secant(
  f : (Double) -> Double,
  x0 : Double,
  x1 : Double,
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 100,
) -> Double raise RootError {
  secant_run(f, x0, x1, tolerance, max_iterations).0
}

///|
/// Newton's method: from a point, step along the tangent to where it crosses
/// zero, over and over. It needs the derivative as well as the function, and
/// it needs a starting point close enough to the root for the tangents to
/// keep pointing at it — but when it works it works faster than anything else
/// here, doubling the number of correct digits with every step.
///
/// A starting point at which the derivative is exactly zero has no tangent,
/// and is refused.
fn newton_run(
  f : (Double) -> Double,
  df : (Double) -> Double,
  x0 : Double,
  tolerance : Double,
  max_iterations : Int,
) -> (Double, Int) raise RootError {
  if tolerance < 0.0 {
    raise RootError::NonPositiveTolerance(tolerance)
  }
  let mut x = x0
  let mut iterations = 0
  let mut converged = false
  while iterations < max_iterations {
    let f_x = f(x)
    if f_x == 0.0 {
      converged = true
      break
    }
    let derivative = df(x)
    if derivative == 0.0 {
      raise RootError::ZeroDerivative(x)
    }
    let step = f_x / derivative
    x = x - step
    iterations = iterations + 1
    if tolerance > 0.0 &&
      (step.abs() <= tolerance * (1.0 + x.abs()) || f_x.abs() <= tolerance) {
      converged = true
      break
    }
  }
  if tolerance > 0.0 && !converged {
    raise RootError::MaxIterations(max_iterations)
  }
  (x, iterations)
}

///|
/// A root of `f`, by Newton's method, from `x0`, with the derivative `df`.
///
/// The starting point may be a root itself, in which case it is returned
/// unchanged. There is no bracket, so the root it finds is the one its
/// starting point leads to; on a function with several roots, start where the
/// root is wanted.
pub fn newton(
  f : (Double) -> Double,
  df : (Double) -> Double,
  x0 : Double,
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 100,
) -> Double raise RootError {
  newton_run(f, df, x0, tolerance, max_iterations).0
}

///|
/// Brent's method: bisection, the secant line and inverse quadratic
/// interpolation, chosen between at every step.
///
/// It keeps a bracket, so it cannot lose the root; it fits a curve through
/// its last three points, so it converges on smooth functions as fast as the
/// secant does; and it falls back to a bisection step whenever the fitted one
/// looks like a bad idea, so it cannot stall or run away. That is why it is
/// the default root finder of most numerical libraries, and why it is the one
/// to use when the function is smooth but a derivative is not at hand.
///
/// The two points may be given in either order.
fn brent_run(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance : Double,
  max_iterations : Int,
) -> (Double, Int) raise RootError {
  if tolerance < 0.0 {
    raise RootError::NonPositiveTolerance(tolerance)
  }
  let f_a = f(a)
  if f_a == 0.0 {
    return (a, 0)
  }
  let f_b = f(b)
  if f_b == 0.0 {
    return (b, 0)
  }
  if f_a * f_b > 0.0 {
    raise RootError::NotBracketed(f_a, f_b)
  }
  // The bracket is carried in increasing order, whichever order it came in.
  let mut left = a
  let mut right = b
  let mut f_left = f_a
  let mut f_right = f_b
  if left > right {
    let swap_x = left
    let swap_f = f_left
    left = right
    f_left = f_right
    right = swap_x
    f_right = swap_f
  }
  // `x` is the best point so far, `previous` the point before it, and
  // `holder` the far end of the bracket; the root is always between `x` and
  // `holder`.
  let mut x = right
  let mut f_x = f_right
  let mut previous = left
  let mut f_previous = f_left
  let mut holder = right
  let mut f_holder = f_right
  let mut step = right - left
  let mut step_before = step
  let mut iterations = 0
  let mut converged = false
  while iterations < max_iterations {
    if f_holder * f_x > 0.0 {
      holder = previous
      f_holder = f_previous
      step = x - previous
      step_before = step
    }
    if f_holder.abs() < f_x.abs() {
      let old_x = x
      let old_f_x = f_x
      previous = old_x
      f_previous = old_f_x
      x = holder
      f_x = f_holder
      holder = old_x
      f_holder = old_f_x
    }
    if tolerance > 0.0 && (holder - x).abs() <= tolerance * (1.0 + x.abs()) {
      converged = true
      break
    }
    if f_x == 0.0 {
      converged = true
      break
    }
    let bracket = 0.5 * (holder - x)
    // Below the rounding of the numbers themselves a step is no step, so
    // every step is measured against the double epsilon as well as against
    // the tolerance the caller asked for.
    let half_tolerance = 2.0 * double_epsilon * x.abs() + 0.5 * tolerance
    if step_before.abs() >= half_tolerance && f_previous.abs() > f_x.abs() {
      // The last step was a long one and the older point is the further from
      // zero, so the two of them are worth fitting a curve through: a
      // straight line when they and the far end leave no room for a curve, an
      // inverse parabola through the three otherwise.
      let slope = f_x / f_previous
      let mut p = 0.0
      let mut q = 0.0
      if previous == holder {
        p = 2.0 * bracket * slope
        q = 1.0 - slope
      } else {
        let a = f_previous / f_holder
        let b = f_x / f_holder
        p = slope * (2.0 * bracket * a * (a - b) - (x - previous) * (b - 1.0))
        q = (a - 1.0) * (b - 1.0) * (slope - 1.0)
      }
      if p > 0.0 {
        q = -q
      }
      p = p.abs()
      if 2.0 * p <
        (3.0 * bracket * q - (half_tolerance * q).abs()).min(
          (step_before * q).abs(),
        ) {
        // The fitted point is inside the bracket and closer to zero by a
        // useful margin, so take it.
        step_before = step
        step = p / q
      } else {
        // It is not, so bisect.
        step = bracket
        step_before = step
      }
    } else {
      step = bracket
      step_before = step
    }
    previous = x
    f_previous = f_x
    if step.abs() > half_tolerance {
      x = x + step
    } else if bracket >= 0.0 {
      x = x + half_tolerance
    } else {
      x = x - half_tolerance
    }
    f_x = f(x)
    iterations = iterations + 1
  }
  if tolerance > 0.0 && !converged {
    raise RootError::MaxIterations(max_iterations)
  }
  (x, iterations)
}

///|
/// The root of `f` between `a` and `b`, by Brent's method.
///
/// The two points have to bracket a root. The answer is within
/// `tolerance · (1 + |x|)` of the root, and usually at the cost of far fewer
/// steps than bisection needs for the same tolerance — six against thirty-two
/// on the tests' functions — with no derivative to supply.
pub fn brent(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 100,
) -> Double raise RootError {
  brent_run(f, a, b, tolerance, max_iterations).0
}
