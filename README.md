# DMRG

DMRG for studying entanglement in 1d systems

## Installation

Requires:

 * (https://github.com/jarro2783/cxxopts)[cxxopts] 
 * (http://itensor.org/)[ITensor]

## Usage

For a list of all options try:

    ./dmrg --help

To run with periodic boundary conditions and a log scale try:


    ./dmrgSacc --pbc -N 13 --Vi=1E-5 --Vf=1 -n 10 --log_scale
