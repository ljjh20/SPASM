# Luke's Special Arbitrary Spring Moddler (SPASM)

_SPASM_ is a set of MATLAB classes and functions to stitch together an arbitrary sequence of curved and straight 'beam' springs and model their collective deflection under loading. It models the strain energy at sample points throughout the spring from all types of external loads, and then models the deflection at those points using Castigliano's second theorem. This provides a useful tool to predict the behaviour of isolated isotropic springs for strain energy storage, and was shown to be a reasonable approximation for laminate composite leaf springs loaded in bending through experimental verification.

## Installation

Place the files in your MATLAB path and run the _stitch_ function found in the `main/models/solver_template.m` example with your spring geometry and material properties stored in the _sections_ hashmap. Make sure you specify your material properties in the _constants_ hashmap and your external forces in the _forces_ hashmap. If you wish to apply a point moment, you can do so by specifying the _couple_ property of a beam section as non-zero.
