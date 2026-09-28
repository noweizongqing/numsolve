///|
/// The command line the demo reads, and one thing read across the package
/// boundary: an optional argument, which travels to this package through the
/// generated interface.

///|
/// Whether two numbers agree to within `tolerance`.
fn cli_close(actual : Double, expected : Double, tolerance : Double) -> Bool {
  (actual - expected).abs() <= tolerance
}

///|
/// The subject of a command that is expected to be a root demo.
fn rooted(arguments : Array[String]) -> String {
  match command(arguments) {
    Root(subject) => subject
    Bad(message) => abort("expected a root demo, not \{message}")
    Help => abort("expected a root demo, not help")
    Ode(_, _) => abort("expected a root demo, not an ode")
    Spline(_) => abort("expected a root demo, not a spline")
    Solve(_) => abort("expected a root demo, not a solve")
  }
}

///|
/// The equation and step count of a command that is expected to be an ode
/// demo.
fn ode_command(arguments : Array[String]) -> (String, Int) {
  match command(arguments) {
    Ode(subject, steps) => (subject, steps)
    Bad(message) => abort("expected an ode, not \{message}")
    Help => abort("expected an ode, not help")
    Root(_) => abort("expected an ode, not a root")
    Spline(_) => abort("expected an ode, not a spline")
    Solve(_) => abort("expected an ode, not a solve")
  }
}

///|
/// The message of a command that is expected to be a mistake.
fn rejected(arguments : Array[String]) -> String {
  match command(arguments) {
    Bad(message) => message
    Help => abort("expected a mistake, not help")
    Root(_) => abort("expected a mistake, not a root")
    Ode(_, _) => abort("expected a mistake, not an ode")
    Spline(_) => abort("expected a mistake, not a spline")
    Solve(_) => abort("expected a mistake, not a solve")
  }
}

///|
/// A command and its subject are read from the arguments: the four demos with
/// their names, the step count of an ode, and help with no arguments at all.
test "the command line is read into commands" {
  assert_eq(rooted(["root", "sqrt2"]), "sqrt2")
  assert_eq(rooted(["root", "cubic"]), "cubic")
  assert_eq(rooted(["root", "trig"]), "trig")
  assert_eq(ode_command(["ode", "growth"]).0, "growth")
  assert_eq(ode_command(["ode", "growth"]).1, 0)
  assert_eq(ode_command(["ode", "decay", "--steps", "100"]).0, "decay")
  assert_eq(ode_command(["ode", "decay", "--steps", "100"]).1, 100)
  assert_eq(
    match command(["spline", "notaknot"]) {
      Spline(subject) => subject
      _ => abort("expected a spline")
    },
    "notaknot",
  )
  assert_eq(
    match command(["solve", "hilbert"]) {
      Solve(subject) => subject
      _ => abort("expected a solve")
    },
    "hilbert",
  )
  assert_true(command([]) is Help)
  assert_true(command(["--help"]) is Help)
  assert_true(command(["-h"]) is Help)
}

///|
/// What the parser cannot read is refused with a message rather than answered
/// with a guess.
test "mistakes are not commands" {
  assert_true(rejected(["root"]).contains("which function"))
  assert_true(rejected(["root", "bogus"]).contains("is not sqrt2"))
  assert_true(rejected(["ode"]).contains("which equation"))
  assert_true(rejected(["ode", "wave"]).contains("is not growth"))
  assert_true(rejected(["ode", "growth", "--steps"]).contains("how many"))
  assert_true(
    rejected(["ode", "growth", "--steps", "x"]).contains("is not a number"),
  )
  assert_true(rejected(["spline"]).contains("which boundary"))
  assert_true(rejected(["spline", "loose"]).contains("is not natural"))
  assert_true(rejected(["solve"]).contains("which system"))
  assert_true(rejected(["solve", "easy"]).contains("is not demo"))
  assert_true(rejected(["frobnicate"]).contains("is not a command"))
}

///|
/// The optional arguments of the solvers travel across the package boundary:
/// this package reads them from the generated interface, so the same function
/// is called with and without the tolerance it offers, and the answers hold
/// the reference values the root tests assert.
test "an optional argument crosses the package boundary" {
  let f = fn(x : Double) { x * x - 2.0 }
  let loose = @numsolve.brent(f, 1.0, 2.0) catch {
    _ => abort("1 and 2 bracket √2")
  }
  let tight = @numsolve.brent(f, 1.0, 2.0, tolerance=1.0e-12) catch {
    _ => abort("1 and 2 bracket √2")
  }
  assert_true(cli_close(loose, 1.4142135623731364, 1.0e-9))
  assert_true(cli_close(tight, 1.4142135623731364, 1.0e-9))
  let adaptive = @numsolve.rk45(fn(_t, y) { y }, 1.0, 0.0, 1.0) catch {
    _ => abort("y' = y is not hard")
  }
  let fixed = @numsolve.rk4(fn(_t, y) { y }, 1.0, 0.0, 1.0, 10) catch {
    _ => abort("y' = y is not hard")
  }
  assert_true(cli_close(adaptive.end_value(), 2.7182818362088534, 1.0e-9))
  assert_true(cli_close(fixed.end_value(), 2.718279744135166, 1.0e-9))
  let spline = @numsolve.CubicSpline::new(
    [0.0, 1.0, 2.0, 3.0],
    [1.0, 2.0, 0.0, 4.0],
    @numsolve.SplineBoundary::Natural,
  ) catch {
    _ => abort("the demonstration knots are valid")
  }
  assert_true(cli_close(spline.evaluate(0.5).unwrap(), 1.95, 1.0e-9))
}
