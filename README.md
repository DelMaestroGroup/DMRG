# DMRG

DMRG for studying accessible entanglement in 1d systems.  Based on the amazing [ITensor](http://itensor.org/) software library.

## Installation

Requires:

 * [cxxopts](https://github.com/jarro2783/cxxopts) 
 * [ITensor](http://itensor.org/)
 * [Boost Format](https://www.boost.org/doc/libs/1_68_0/libs/format/)

You will need to modify the `CCFLAGS` variable in `Makefile` to point to where
you have installed the header files.

## Usage

For a list of all options try:

    ./dmrg --help

    Accessible Entanglement via DMRG
    Usage:
      ./dmrgSacc [OPTION...]

          --pbc                    periodic boundary conditions
          --random                 use a random initial state
          --log_scale              use a logarithmic spacing for the interaction
                                   strength
      -n, --num_interaction n      the number of interaction points
      -N, --number_particles N     number of particles
      -V, --interaction_stength V  interaction strength
      -i, --Vi Vi                  initial interaction strength
      -f, --Vf Vf                  final interaction strength
          --dV dV                  interaction strength step
          --help                   Print help

To run with periodic boundary conditions and a log scale try:


    ./dmrgSacc --pbc -N 13 --Vi=1E-5 --Vf=1 -n 10 --log_scale
