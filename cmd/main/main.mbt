///|
/// The solvers, one command at a time.
///
///     moon run cmd/main -- root sqrt2
///     moon run cmd/main -- ode growth --steps 100
///     moon run cmd/main -- spline natural
///     moon run cmd/main -- solve hilbert
///
/// Each command runs one solver on one problem and prints its answer next to
/// the number the problem is known to have, so the two can be compared.

///|
/// What the command line asked for.
priv enum Command {
  Root(String)
  Ode(String, Int)
  Spline(String)
  Solve(String)
  Help
  Bad(String)
}

///|
fn main {
  match command(@env.args()) {
    Help => println(usage())
    Bad(message) => {
      println(message)
      println("")
      println(usage())
    }
    Root(name) => root_demo(name)
    Ode(name, steps) => ode_demo(name, steps)
    Spline(name) => spline_demo(name)
    Solve(name) => solve_demo(name)
  }
}

///|
/// How the demo is used.
fn usage() -> String {
  [
    "numsolve — the solvers, one command at a time", "", "  moon run cmd/main -- <command> <subject> [--steps N]",
    "", "  root sqrt2|cubic|trig      a root by Brent or false position", "  ode growth|decay           y' = y or y' = -y over [0, 1] by rk45,",
    "                             or by rk4 when --steps N is given", "  spline natural|clamped|notaknot",
    "                             the three cubic splines through (0,1),", "                             (1,2), (2,0), (3,4): every piece's",
    "                             coefficients and three samples", "  solve demo|hilbert         the 3x3 demo system, or Hilbert 8 whose",
    "                             solution is known to be all ones", "  --help                     this",
  ].join("\n")
}

///|
/// The next argument after `index`, if there is one.
fn following(arguments : Array[String], index : Int) -> String? {
  if index + 1 < arguments.length() {
    Some(arguments[index + 1])
  } else {
    None
  }
}

///|
/// The command the arguments describe.
fn command(arguments : Array[String]) -> Command {
  if arguments.length() == 0 {
    return Help
  }
  let name = arguments[0]
  if name == "--help" || name == "-h" {
    return Help
  }
  if name == "root" {
    let subject = match following(arguments, 0) {
      None => return Bad("root: which function? sqrt2, cubic or trig")
      Some(value) => value
    }
    if subject == "sqrt2" || subject == "cubic" || subject == "trig" {
      return Root(subject)
    }
    return Bad("root: \{subject} is not sqrt2, cubic or trig")
  }
  if name == "ode" {
    let subject = match following(arguments, 0) {
      None => return Bad("ode: which equation? growth or decay")
      Some(value) => value
    }
    if subject != "growth" && subject != "decay" {
      return Bad("ode: \{subject} is not growth or decay")
    }
    let mut steps = 0
    if arguments.length() >= 3 && arguments[2] == "--steps" {
      let count = match following(arguments, 2) {
        None => return Bad("--steps: how many?")
        Some(value) => value
      }
      steps = @string.parse_int(count) catch {
        _ => return Bad("--steps: \{count} is not a number")
      }
    }
    return Ode(subject, steps)
  }
  if name == "spline" {
    let subject = match following(arguments, 0) {
      None => return Bad("spline: which boundary? natural, clamped or notaknot")
      Some(value) => value
    }
    if subject == "natural" || subject == "clamped" || subject == "notaknot" {
      return Spline(subject)
    }
    return Bad("spline: \{subject} is not natural, clamped or notaknot")
  }
  if name == "solve" {
    let subject = match following(arguments, 0) {
      None => return Bad("solve: which system? demo or hilbert")
      Some(value) => value
    }
    if subject == "demo" || subject == "hilbert" {
      return Solve(subject)
    }
    return Bad("solve: \{subject} is not demo or hilbert")
  }
  Bad("\{name} is not a command. Try --help.")
}

///|
/// A root of one of the three demonstration functions, against the value it
/// is known to have.
fn root_demo(name : String) -> Unit {
  if name == "sqrt2" {
    let root = @numsolve.brent(fn(x) { x * x - 2.0 }, 1.0, 2.0) catch {
      _ => abort("1 and 2 bracket the root of x² − 2")
    }
    println("x² − 2 = 0 at \{root}")
  } else if name == "cubic" {
    let root = @numsolve.brent(fn(x) { x * x * x - x - 2.0 }, 1.0, 2.0) catch {
      _ => abort("1 and 2 bracket the root of x³ − x − 2")
    }
    println("x³ − x − 2 = 0 at \{root}")
  } else {
    let root = @numsolve.false_position(fn(x) { @math.sin(x) - 0.5 }, 0.0, 1.0) catch {
      _ => abort("0 and 1 bracket sin x = 1/2")
    }
    println("sin x = 1/2 at \{root}, and π/6 = \{@math.PI / 6.0}")
  }
}

///|
/// A trajectory of `y' = y` or `y' = -y` over `[0, 1]` from 1, by rk45 — or by
/// rk4 with the number of steps given — against e or 1/e.
fn ode_demo(name : String, steps : Int) -> Unit {
  let rhs : (Double, Double) -> Double = match name {
    "growth" => fn(_t, y) { y }
    _ => fn(_t, y) { -y }
  }
  let truth = if name == "growth" {
    @math.exp(1.0)
  } else {
    1.0 / @math.exp(1.0)
  }
  if steps > 0 {
    let trajectory = @numsolve.rk4(rhs, 1.0, 0.0, 1.0, steps) catch {
      _ => abort("y' = y is not a hard equation")
    }
    let value = trajectory.end_value()
    println(
      "\{steps} steps of rk4: \{value} (truth \{truth}, error \{(value - truth).abs()})",
    )
  } else {
    let trajectory = @numsolve.rk45(rhs, 1.0, 0.0, 1.0) catch {
      _ => abort("y' = y is not a hard equation")
    }
    let value = trajectory.end_value()
    println(
      "rk45: \{value} (truth \{truth}, error \{(value - truth).abs()}, \{trajectory.length()} samples)",
    )
  }
}

///|
/// The three splines through (0,1), (1,2), (2,0), (3,4): the coefficients of
/// every piece and the three midpoints sampled.
fn spline_demo(name : String) -> Unit {
  let xs = [0.0, 1.0, 2.0, 3.0]
  let ys = [1.0, 2.0, 0.0, 4.0]
  let boundary = match name {
    "natural" => @numsolve.SplineBoundary::Natural
    "clamped" => @numsolve.SplineBoundary::Clamped(0.0, 0.0)
    _ => @numsolve.SplineBoundary::NotAKnot
  }
  let spline = @numsolve.CubicSpline::new(xs, ys, boundary) catch {
    _ => abort("the demonstration knots are valid")
  }
  let rows = spline.coefficients()
  let mut i = 0
  while i < rows.length() {
    println(
      "piece \{i}: [\{rows[i][0]}, \{rows[i][1]}, \{rows[i][2]}, \{rows[i][3]}]",
    )
    i = i + 1
  }
  for x in [0.5, 1.5, 2.5] {
    println("S(\{x}) = \{spline.evaluate(x).unwrap()}")
  }
}

///|
/// The 3×3 demo system, or the Hilbert matrix of order 8 solved against the
/// row sums, whose solution is known to be all ones.
fn solve_demo(name : String) -> Unit {
  if name == "demo" {
    let a = [[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]
    let x = @numsolve.solve_linear_system(a, [1.0, 2.0, 4.0]) catch {
      _ => abort("the demo system is not singular")
    }
    println("A·x = b at x = [\{x[0]}, \{x[1]}, \{x[2]}]")
  } else {
    let size = 8
    let hilbert : Array[Array[Double]] = Array(capacity=size)
    let mut i = 0
    while i < size {
      let row : Array[Double] = Array(capacity=size)
      let mut j = 0
      while j < size {
        row.push(1.0 / (i + j + 1).to_double())
        j = j + 1
      }
      hilbert.push(row)
      i = i + 1
    }
    let rhs : Array[Double] = Array(capacity=size)
    let mut row = 0
    while row < size {
      let mut sum = 0.0
      let mut j = 0
      while j < size {
        sum = sum + hilbert[row][j]
        j = j + 1
      }
      rhs.push(sum)
      row = row + 1
    }
    let x = @numsolve.solve_linear_system(hilbert, rhs) catch {
      _ => abort("Hilbert matrices are not singular")
    }
    let mut error = 0.0
    for value in x {
      error = error.max((value - 1.0).abs())
    }
    println("Hilbert 8 (condition 1.5e10): x is all ones to \{error}")
  }
}
