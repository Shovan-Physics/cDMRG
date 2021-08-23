# cDMRG

cDMRG enables the application of the Density-Matrix Renormalization Group (DMRG) to continuous quantum
systems by spatial partitioning into segments with continuous basis functions.
It is based on [arXiv:2108.05366](https://arxiv.org/abs/2108.05366) and distinct from cMPS.
The code is written for bosons with contact interactions trapped in a box with a sinusoidal potential,
divided into segments of equal width.

## Installation

In order to get this code running on your machine follow these steps.

### Prerequisites

- We've tested this package on Unix only, however running on Win10 should be straightforward
  with the build tools installed.
  We recommend the package manager for Windows [Chocolatey](https://chocolatey.org/),
  then you can simply run `choco install make`.
- [Install and build ITensor 3 (C++ version)](https://itensor.org/).
- [Install Mathematica 12.3](https://www.wolfram.com/mathematica/) or newer.
  We did not test older versions and it's likely they work as well.

### Build

- Clone this repository on you machine.
- Update `DIR` in `Makefile` that it points to the ITensor directory.
- Build binaries:
  ```console
  make build
  ```
- Among other files you should find `cDMRG` executable that will be used by the front-end Mathematica script,
  which are described below.

## Code Structure

cDMRG consists of two main modules:
- DMRG Module using the ITensor library 
- Input Module written in Mathematica

### DMRG Module

The DMRG module consists of the following files:

- `cDMRG.cc` - main block
- `cio.h`, `cio.cc` - input/output
- `readvec.h`, `readvec.cc` - reading vectors
- `cSite.h` - site set
- `cDMRGeps.h` - DMRG sweeps
- `localmpoproj.h` - Hamiltonian MPOs

### Input Module

The DMRG module is called by the Mathematica package `cDMRG_input.m`, which generates the necessary
inputs (local basis and operators) for DMRG.
A separate documentation of `cDMRG_input.m` is provided in the notebook `cDMRG_input_documentation.nb`.
The package can be modified by changing the notebook `cDMRG_input.nb` and saving it as a package `.m`.

### IDs / Stored Inputs

Parameters specifying the local basis, DMRG sweeps, and what results to save are stored in
`basistable.m`, `epstable.m`, and `savetable.m` in the `Parameters` directory,
with unique integer pointers that act as IDs.
One can add more entries using a text editor or the function `storenewparam` in `cDMRG_inputIDgen.nb`.

## How to Use

### Store/find Input IDs

Use the notebook `cDMRG_inputIDgen.nb` to 
(1) specify parameters for the local basis, DMRG parameters, and what results to save.
(2) check if these are already stored as inputs and find their IDs; if not, create new IDs.

### Run the Program

Call the program from the command line as 
```consol
math -noprompt -run '<<cDMRG_input.m' N gamma Nwell V0 M basisid epsid saveid &
```

where

- `N` - `int` number of particles
- `gamma` - `float` dimensionless interaction strength
- `Nwell` - `float` number of potential minima
- `V0` - `float` potential depth in units of recoil
- `M` - `int` number of segments
- `basisid` - `int` ID for basis parameters
- `epsid` - `int` ID for DMRG parameters
- `saveid` - `int` ID for save parameters

For instance a short calculation can be initiated by
```consol
math -noprompt -run '<<cDMRG_input.m' 4 10 4 1 8 1 1 1 &
```

The results are stored in the directory `Runs`.

### Visualize Outputs

Results from a sample run are stored in the `Runs` directory and illustrated with plots in the notebook `cDMRG_example_output.nb`.

### Basis splitting

Codes for splitting a segment at an arbitrary point, and calculating the bipartite entanglement
are provided in the notebook `cDMRG_basis_splitting.nb`.
An example is given using the output MPS stored in `Runs`.
