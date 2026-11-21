# Wing Structural Optimization
A set of tools for wing structural optimization in Matlab.
This code was orginally developed for the SAE Aerodesign Brasil competition at
the UFMG's Aerodesign Team (Uai, Sô! Fly!).

## TODO
 - [ ] Refactor the code to be more modular.
 - [ ] Full documentation corverge.

## Improvements
 - [ ] Implement others beam cross-sections.
 - [ ] Implement composite ply analysis.

### Follow-up work (not part of MATLAB parity)

The Julia port preserves the current MATLAB safety/uncertainty factors,
half-wing sign conventions, and torsion calculations. Wing sweep, geometric
airfoil torsion, and airfoil-based aerodynamic torsion remain unsupported
follow-up improvements; the existing geometry fields do not apply these
effects to the load model.

The thin-wall MATLAB helper currently computes `sum(thickness × length)`
instead of the dimensionally standard `sum(length / thickness)`. The Julia
port preserves this numerical behavior for parity; correcting it is a
follow-up because it changes torsion results.

## Next Version
 - [ ] Migration to Python.
 - [ ] Migration to OpenMDAO.

## References

 - [Greco, M. - Resistência dos Materiais: Uma Abordagem Sintética](https://www.amazon.com.br/Resist%C3%AAncia-dos-Materiais-Marcelo-Greco/dp/8535274588)
 - [Megson, T.H.G. - Aircraft Structures for Engineering Students](https://www.amazon.com.br/Aircraft-Structures-Engineering-Students-T-H-G/dp/0081009143)
 - Airframe and Equipment Engineering Report No. 45, AD-A955 270, 1955
 - [Mathworks Matlab Optimization Docs](https://www.mathworks.com/help/gads/index.html?s_tid=CRUX_lftnav)
 
 
 
