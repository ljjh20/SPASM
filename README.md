# Luke's Special Arbitrary Spring Modeller (SPASM)

![alt text](https://github.com/ljjh20/SPASM/blob/master/main/media/spasm_interface.png)
_SPASM_ is a set of MATLAB classes and functions to stitch together an arbitrary sequence of curved and straight beam sections model their collective deflection, stiffness, and energetics under point loading. It models the strain energy at sample points throughout the spring, and then models the deflection at those points using Castigliano's second theorem. This provides a useful tool to predict the behaviour of isolated isotropic springs for strain energy storage, and was shown to be a reasonable approximation for laminate composite leaf springs loaded in bending through experimental verification.

## Installation

Place the files in your MATLAB path and run the `template.mlx` script found in `main`  with your spring geometry and material properties stored in the `sections` hashmap. Make sure you specify your material properties in the `constants` hashmap and your external forces in the `forces` hashmap. If you wish to apply a point moment, you can do so by specifying the `couple` property of a beam section as non-zero.

## Sweep Functionality

If you wish to sweep across a variable, pass your hashmap of beam sections into the `stitch_sweep` function instead of the `stitch` function, and define your variable name as a string, number of sections as an int, and the range of values as a linspace input. 

![alt_text](https://github.com/ljjh20/SPASM/blob/master/main/media/spasm_sweep.png)
