# General information
This project is a modification of [IMFIT](http://www.mpe.mpg.de/~erwin/code/imfit/) ([Erwin, 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...799..226E/abstract)) with added models of spiral arms (firstly introduced in [Chugunov et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9605C/abstract)). The only difference from basic IMFIT 1.9.0 is that some models of spiral arms are implemented, and all other features remain the same.

In particular, this package has exactly the same dependencies as basic IMFIT, and the installation process is also has no differences (the instructions from the basic IMFIT can be found at docs/imfit_howto.pdf). This package do not require basic IMFIT to be installed.

# Spiral arms model
The most comprehensive documentation for the added spiral arm models is provided at [docs_spirals/model_description.pdf](docs_spirals/model_description.pdf).

The basic model is `SpiralArm0b`. This function produces 2D light distribution in the individual spiral arm. Produced spiral arms may have variable pitch angle, variable width and asymmetric perpendicular profile. Here, we provide a list of parameters of this function.
* `X0`, `Y0`: image coordinates of a center of the spiral structure (similar to other IMFIT components, generally should match the disc coordinates);
* `PA`, `ell`: parameters describing the orientation of the galactic plane (similar to other IMFIT components, generally should match the disc coordinates);
* `r_0`, `phi_0`: spiral arm beginning position in polar coordinates. (`phi_0` is counted counterclockwise from the PA direction). [`r_0` in pixels, `phi_0` in degrees].
* `r_end`, `phi_end`: spiral arm ending position in polar coordinates. (`phi_end` is counted counterclockwise from the PA direction). [`r_end` in pixels, `phi_end` in degrees].
* `mu_a_2`, `mu_a_3`, `mu_a_4`: coefficients describing the deviation of spiral arm from pure logarighmic spiral shape. [Dimensionless; reasonable values are positive or negative of the order of unity or severel units; all zero indicates logarithmic spiral].
* `I_0`: spiral arm surface brightness projected to the center. [In image counts].
* `part_growth`: part of azimuthal length of the arm where growth from zero brightness to the beginning of exponentially decreasing part occurs. [Dimensionless; reasonable values are positive and of the order of 1/10].
* `h_s`: radial exponential scale of spiral arm (defined in a similar way as for exponential disc). [In pixels].
* `part_cutoff`: part of the azimuthal length of the arm where decrease from the ending of exponential part to zero brightness occurs. [Dimensionless; reasonable values are positive and of the order of 1/10].
* `width`: radial arm width at the middle point (at the middle azimuthal length). [In pixels].
* `w_asymm`: spiral arm asymmetry at the middle point. [Dimensionless; reasonable values are from -1 to 1, and 0 indicates symmetry].
*	`n_out`, `n_in`: outer and inner Sersic index of the perpendicular profile of the spiral arm. [Dimensionless; reasonable values are positive and less or of the order of unity].
*	`gamma_out`, `gamma_in` --- outer and inner widening coefficients (inner/outer half-widths at the ending of arm are by a factor of exp(`gamma_in`) and exp(`gamma_out`) greater than inner/outer half-widths at the beginning, respectively). [Dimensionless; reasonable values are less or of the order of unity; both 0 indicate constant width, negative indicates decreasing width].
