///|
fn vector_close(
  actual : Array[Double],
  expected : Array[Double],
  tolerance : Double,
) -> Bool {
  if actual.length() != expected.length() {
    return false
  }
  let mut i = 0
  while i < actual.length() {
    if (actual[i] - expected[i]).abs() > tolerance {
      return false
    }
    i = i + 1
  }
  true
}

///|
fn harmonic_oscillator(_t : Double, y : Array[Double]) -> Array[Double] {
  [y[1], -y[0]]
}

///|
fn coupled_decay(_t : Double, y : Array[Double]) -> Array[Double] {
  [y[0], -2.0 * y[1]]
}

///|
test "rk4_system integrates a two-dimensional oscillator" {
  let period = 6.283185307179586
  let trajectory = rk4_system(
    harmonic_oscillator,
    [1.0, 0.0],
    0.0,
    period,
    2000,
  ) catch {
    _ => abort("the harmonic oscillator is a valid system")
  }
  assert_eq(trajectory.dimension(), 2)
  assert_eq(trajectory.length(), 2001)
  assert_true(vector_close(trajectory.end_state(), [1.0, 0.0], 1.0e-9))
}

///|
test "the vector fixed-step methods converge at their advertised order" {
  let exact = @math.exp(1.0)
  let euler_error = {
    let result = euler_system(coupled_decay, [1.0, 1.0], 0.0, 1.0, 40) catch {
      _ => abort("the exponential system is valid")
    }
    (result.end_state()[0] - exact).abs()
  }
  let midpoint_error = {
    let result = midpoint_system(coupled_decay, [1.0, 1.0], 0.0, 1.0, 40) catch {
      _ => abort("the exponential system is valid")
    }
    (result.end_state()[0] - exact).abs()
  }
  let rk4_error = {
    let result = rk4_system(coupled_decay, [1.0, 1.0], 0.0, 1.0, 40) catch {
      _ => abort("the exponential system is valid")
    }
    (result.end_state()[0] - exact).abs()
  }
  assert_true(rk4_error < midpoint_error)
  assert_true(midpoint_error < euler_error)
  assert_true(rk4_error < 1.0e-7)
}

///|
test "rk45_system uses absolute and relative component tolerances" {
  let trajectory = rk45_system(
    coupled_decay,
    [1.0, 1.0],
    0.0,
    1.0,
    rtol=1.0e-9,
    atol=1.0e-12,
  ) catch {
    _ => abort("the exponential system is valid")
  }
  let expected = [@math.exp(1.0), @math.exp(-2.0)]
  assert_true(vector_close(trajectory.end_state(), expected, 1.0e-7))
  assert_true(trajectory.length() > 1)
}

///|
test "solve_ivp selects the fixed and adaptive methods" {
  let fixed = solve_ivp(
    coupled_decay,
    [1.0, 1.0],
    0.0,
    1.0,
    solver=OdeMethod::RK4,
    steps=100,
  ) catch {
    _ => abort("the exponential system is valid")
  }
  let adaptive = solve_ivp(
    coupled_decay,
    [1.0, 1.0],
    0.0,
    1.0,
    solver=OdeMethod::RK45,
    rtol=1.0e-8,
    atol=1.0e-10,
  ) catch {
    _ => abort("the exponential system is valid")
  }
  let expected = [@math.exp(1.0), @math.exp(-2.0)]
  assert_true(vector_close(fixed.end_state(), expected, 1.0e-7))
  assert_true(vector_close(adaptive.end_state(), expected, 1.0e-6))
}

///|
test "vector integration can run backwards" {
  let trajectory = rk4_system(
    harmonic_oscillator,
    [0.0, -1.0],
    1.5707963267948966,
    0.0,
    500,
  ) catch {
    _ => abort("the harmonic oscillator is a valid system")
  }
  assert_true(vector_close(trajectory.end_state(), [1.0, 0.0], 1.0e-9))
}

///|
test "a derivative with the wrong dimension is refused" {
  let failed = try {
    let _ = rk4_system(
      fn(_t : Double, _y : Array[Double]) { [1.0] },
      [1.0, 2.0],
      0.0,
      1.0,
      1,
    )
    false
  } catch {
    OdeError::DimensionMismatch(expected, actual) =>
      expected == 2 && actual == 1
    _ => false
  }
  assert_true(failed)
}

///|
test "an empty vector state is refused" {
  let failed = try {
    let _ = rk45_system(harmonic_oscillator, [], 0.0, 1.0)
    false
  } catch {
    OdeError::EmptyState => true
    _ => false
  }
  assert_true(failed)
}
