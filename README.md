# Wing Structural Optimization

Wing structural optimization tools originally developed in MATLAB for the SAE
Aerodesign Brasil competition at the UFMG Aerodesign Team (Uai, Sô! Fly!). The
legacy MATLAB implementation remains at the repository root; the Julia package
is under `WingOpt/`.

## Julia package

The package separates geometry, loads, structural analysis, aerodynamics,
optimization, input/output, and visualization under the `WingOpt` module.
Install dependencies and run the test suite from a clean checkout with:

```sh
julia --project=WingOpt -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=WingOpt WingOpt/examples/smoke_test.jl
```

Inputs and generated results can be represented as TOML-compatible dictionaries
using `WingOpt.InputOutput.read_data` and `WingOpt.InputOutput.write_data`.
`WingOpt/Manifest.toml` records the resolved Julia dependency versions.

Direct dependencies are limited to the package's numerical and development
needs: `LinearAlgebra` for matrix operations, `Interpolations` for airfoil and
station interpolation, `Metaheuristics` and `Optim` for optimization methods,
`Plots` for visualization, and Julia's standard-library `TOML` for portable
input and result files. Formatting follows `WingOpt/.JuliaFormatter.toml` and
is checked in CI.

## Improvements
 - [ ] Implement others beam cross-sections.
 - [ ] Implement foil geometric torsion.
 - [ ] Implement wing sweep.
 - [ ] Implement foil aerodynamic torsion.
 - [ ] Implement composite ply analysis.

## References

 - [Greco, M. - Resistência dos Materiais: Uma Abordagem Sintética](https://www.amazon.com.br/Resist%C3%AAncia-dos-Materiais-Marcelo-Greco/dp/8535274588)
 - [Megson, T.H.G. - Aircraft Structures for Engineering Students](https://www.amazon.com.br/Aircraft-Structures-Engineering-Students-T-H-G/dp/0081009143)
 - Airframe and Equipment Engineering Report No. 45, AD-A955 270, 1955
 - [Mathworks Matlab Optimization Docs](https://www.mathworks.com/help/gads/index.html?s_tid=CRUX_lftnav)
 
 
 
