# General information
This project is a modification of [IMFIT](http://www.mpe.mpg.de/~erwin/code/imfit/) ([Erwin, 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...799..226E/abstract)) with added models of spiral arms (see [Chugunov et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9605C/abstract), [Chugunov et al., 2025](https://ui.adsabs.harvard.edu/abs/2025Galax..13...44C/abstract) and related works). The only difference from basic IMFIT 1.9.0 is that some models of spiral arms are implemented, and all other features remain the same.

In particular, this package has exactly the same dependencies as basic IMFIT, and the installation process is also has no differences (the instructions from the basic IMFIT can be found at docs/imfit_howto.pdf). This package do not require basic IMFIT to be installed.

# Showcase examples

WIP

# Spiral arms model
Our model produces 2D light distribution in the individual spiral arm. Produced spiral arms may have variable pitch angle, variable width and asymmetric perpendicular profile. The most comprehensive documentation for the added spiral arm models is provided at [docs_spirals/model_description.pdf](docs_spirals/model_description.pdf). The properties of model are also described in [Chugunov et al., 2025b](https://ui.adsabs.harvard.edu/abs/2025Galax..13...44C/abstract), but the documentation contains more technical details and describes all varieties of the function implementation for convenient use.

Below, we provide a schematic illustration of some properties of spiral arms. (WIP)

## List of parameters

Here is the brief description of all parameters of our baseline model, `SpiralArm`.
The basic model is `SpiralArm0b`.  Here, we provide a list of parameters of this function.
* `X0`, `Y0`: image coordinates of a center of the spiral structure (similar to other IMFIT components, generally should match the disc coordinates);
* `PA`, `ell`: parameters describing the orientation of the galactic plane (similar to other IMFIT components, generally should match the disc coordinates);
* `r_0`, `phi_0`: spiral arm beginning position in polar coordinates.
* `r_end`, `phi_end`: spiral arm ending position in polar coordinates.
* `mu_a_2`, `mu_a_3`, `mu_a_4`: coefficients defining the deviation of spiral arm from pure logarighmic spiral shape.
* `I_0`: spiral arm surface brightness projected to the center.
* `part_growth`, `part_cutoff`: parts of azimuthal length of the arm where growth from zero brightness to exponential part, and decline from exponential to zero occur.
* `ih_s`: inverse radial exponential scale of spiral arm.
* `w_zp`, `w_i`: linear coefficients for arm width as a function of galactocentric radius.

# Version history
The history of development of this package is complicated. The first model of spiral arms by our group was introduced in [Chugunov et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9605C/abstract). This model (Alpha) was constructed without thorough validation, being essentially heuristic (however, it worked well anyway). Later, some of its properties were reconsidered based on workflow experience and for some time we used different model (Beta). Only after that we conducted a proper investigation of spiral arms surface brightness profile ([Chugunov et al., 2025b](https://ui.adsabs.harvard.edu/abs/2025Galax..13...44C/abstract)) which allowed us to construct more reliable model of spiral arms. Nevertheless, most of the core properties and reasoning behind them do not differ much between models.

Thus, we consider models from our first works to be obsolete and do not plan (and do not recommend) using them anymore. However, some of our published papers rely on these obsolete models, and we decided to keep them in separate branches for this reason. Below, we present a list of our papers which used obsolete versions of IMFIT_spirals.

* [Alpha](https://github.com/IVChugunov/IMFIT_spirals/tree/Alpha) branch:
** [Chugunov et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9605C/abstract)
** [Marchuk et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.528.1276M/abstract)
* [Beta](https://github.com/IVChugunov/IMFIT_spirals/tree/Beta) branch:
** [Chugunov et al., 2025a](https://ui.adsabs.harvard.edu/abs/2025PASA...42...29C/abstract) (in particular, models in [the corresponding repository](https://github.com/IVChugunov/Distant_spirals_decomposition) need this version)
** [Kostiuk et al., 2025](https://ui.adsabs.harvard.edu/abs/2025Galax..13...27K/abstract)
** [Marchuk et al., 2025](https://ui.adsabs.harvard.edu/abs/2025Galax..13...39M/abstract)
* [main](https://github.com/IVChugunov/IMFIT_spirals/tree/main) branch contains our final version.