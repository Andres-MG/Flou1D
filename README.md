# Flou1D

Flou1D is a one-dimensional solver for the Navier-Stokes equations. It was conceived as a research
tool and thus, ease of use has been favored over performance.

## Compilation instructions

Compilation is handled by `CMake`:

1. Create a separate folder for the build:

    * `mkdir build`

2. Move to the new folder:

    * `cd build`

3. Generate the compilation files:

    * `ccmake ..`

4. Compile the code:

    * `make [-j #nprocs]`

Other compilation tools can be used in step 4 instead of `make`. More information in this
[link](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html).

<!-- ### Forbidden keywords

* PDE
* TimeIntegrator
* Sensor
* TruncError
* TwoPointFlux
* DissipativeFlux
* RiemannSolver
* ArtViscousFlux

### Working on

* Sensor for WENO: tau or not??
* WENO epsilon =? dx
* WENO only with Dirichlet BDs
* Sensor with integeral?? -->
