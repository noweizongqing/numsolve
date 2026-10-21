# Contributing to numsolve

## Development setup

Install the current MoonBit toolchain, then clone the repository and run the
checks from the project root:

```bash
moon version --all
moon fmt --check
moon check --deny-warn
moon test
moon bench --release --target wasm
```

Use `moon fmt` before committing source changes. The generated interface is
checked in; run `moon info` when a public API changes.

## Numerical changes

- Keep algorithms tied to a published method or a mathematical definition.
- Add tests for the claimed property, not only for one final value.
- Test boundary cases: empty inputs, wrong dimensions, singular systems,
  invalid tolerances, and non-convergence where applicable.
- Keep scalar convenience APIs stable when adding a system API.
- Document the complexity and the limits of the method.

## Documentation

Update `README.mbt.md` for user-facing behavior and `docs/ode.md` or
`docs/linear.md` for design details. Keep examples runnable and avoid claiming
compatibility that is not tested.

## Commits and releases

Use a descriptive prefix such as `feat(ode):`, `fix(linear):`, `docs:`, or
`chore(release):`. Include tests and documentation with the implementation.
Before publishing, verify the CI commands above and update `CHANGELOG.md`.

## License

By contributing, you agree that your contribution is licensed under the
Apache License 2.0 in `LICENSE`.
