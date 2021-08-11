# cDMRG

cDMRG enables the application of the Density-Matrix Renormalization Group (DMRG) to continuous quantum
systems by spatial partitioning into segments with continuous basis functions.
It is distinct from cMPS.
The code is written for bosons with contact interactions trapped in a box with a sinusoidal potential.

## Modules

cDMRG consists of two main modules:
- DMRG module using the ITensor library 
- input module written in Mathematica

### DMRG Module

The DMRG module is built from 3 source files and 5 header files:

- `cDMRG.cc` - main block
- `cio.h`, `cio.cc` - input/output
- `readvec.h`, `readvec.cc` - reading vectors
- `cSite.h` - site set
- `cDMRGeps.h` - DMRG sweeps
- `localmpoproj.h` - Hamiltonian MPOs

First install the C++ version of ITensor: https://itensor.org/.
Change the address in the make file `make_cDMRG` to point to the built ITensor library.
Then run the make file from command line. It should produce an executable file named `cDMRG`.

### Input Module

The DMRG module is called by the Mathematica package `cDMRG_input.m`, which generates the necessary
inputs (local basis and operators) for DMRG.
A separate documentation of `cDMRG_input.m` is provided in the notebook `cDMRG_input_documentation.nb`.
The package can be modified by changing the notebook `cDMRG_input.nb` and saving it as a package `.m`.

One calls the program from command line as 

```consol
math -noprompt -run '<<cDMRG_input.m' N gamma Nwell V0 M basisid epsid saveid &
```

where

- `N` - number of particles
- `gamma` - dimensionless interaction strength
- `Nwell` - number of potential minima
- `V0` - potential depth in units of recoil
- `M` - number of segments
- `basisid` - ID for basis parameters
- `epsid` - ID for DMRG parameters
- `saveid` - ID for save parameters.

The results are stored in the directory `Runs`.

## IDs / Stored Inputs

The parameters specifying the local basis, DMRG sweeps, and what results to save are stored in
`basistable.m`, `epstable.m`, and `savetable.m` in the `Parameters` directory,
with unique integer pointers that act as IDs.
One can add more entries using a text editor or the function `storenewparam` in `cDMRG_input_documentation.nb`.

## Sample Outputs

Results from a sample run are stored in the `Runs` directory and illustrated with plots in the notebook `cDMRG_example_output.nb`.

## Basis splitting

Codes for splitting a segment at an arbitrary point, and calculating the bipartite entanglement
are provided in the notebook `cDMRG_basis_splitting.nb`.
An example is given using the output MPS from the DMRG module.
