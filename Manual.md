# Manual


Besides the routine inputs in `CONTROL` file of `ShengBTE`, `FourPhonon` requires a fourth-order force constants and some new name-lists in `CONTROL` file. Check [ShengBTE website](https://bitbucket.org/sousaw/shengbte/src/master/README.md) for definition of other name-lists.

## CPU Parallel environment

FourPhonon uses **MPI+OpenMP hybrid parallelism** for efficient parallel execution on modern HPC systems:

- **MPI**: Distributes work across nodes and processes
- **OpenMP**: Provides thread-level parallelism within each MPI process for memory-intensive operations

### Compilation flags

Configure `Src/arch.make` for CPU compilation:

```makefile
export CPU_COMPILER = mpiifort
CPU_FFLAGS = -qopenmp -traceback -O2 -fpp -DCPU_VERSION
```

Make sure to add `-qopenmp` (Intel) for OpenMP support.

### Example CPU compilation

```bash
cd Src/
make cpu
# or simply
make
```

This creates the `ShengBTE_cpu` executable in the root directory.

### Running with MPI+OpenMP

Set the number of MPI processes and OpenMP threads:

```bash
export OMP_NUM_THREADS=8
export OMP_STACKSIZE=1G
mpirun -np 4 ./ShengBTE_cpu
```

This example uses 4 MPI processes with 8 OpenMP threads each (total 32 cores).

**Important**: `OMP_STACKSIZE=1G` is required on some clusters to allocate sufficient stack memory for OpenMP threads.

**Note**: Earlier versions (1.0-1.1) used MPI-only parallelism. Version 1.2 migrated to OpenMP support for the iterative solver to handle large memory requirements.

## GPU acceleration

FourPhonon supports GPU acceleration using OpenACC directives for NVIDIA GPUs. The GPU version is compiled separately from the CPU version and provides significant speedup for large systems. For more details on GPU implementation, see [arXiv:2510.00518](https://arxiv.org/abs/2510.00518).

### GPU compilation requirements

- **Compiler**: NVIDIA HPC SDK with `mpif90` compiler
- **GPU architecture**: NVIDIA GPU with compute capability 6.0 or higher
- **OpenACC**: Version 2.7 or later support required

### GPU compilation flags

Configure `Src/arch.make` for GPU compilation:

```makefile
export GPU_COMPILER = mpif90
GPU_FFLAGS = -acc -Minfo=accel -Mpreprocess -g -cuda -gpu=ptxinfo,cc80 -mp -O2 -DGPU_VERSION -DGPU_ALL_MODE_PARALLELIZATION
```

### GPU parallelization modes

Two mutually exclusive GPU parallelization strategies are available:

- `GPU_ALL_MODE_PARALLELIZATION` (default): Parallelizes across all phonon modes simultaneously. Recommended for most systems.
- `GPU_MODE_BY_MODE_PARALLELIZATION`: Parallelizes mode-by-mode. May be more memory efficient for very large systems.

### Current GPU restriction

- **RTA calculations only**: The current GPU version only supports Relaxation Time Approximation (RTA) calculations. Iterative BTE solvers (`convergence=.true.` and `four_phonon_iteration=.true.`) are not implemented on GPU
- **Sampling methods not supported**: The GPU version does not support sampling acceleration methods (`num_sample_process_*ph` parameters must be `-1`)
- **Memory requirements**: GPU calculations require sufficient GPU memory. Large systems may need high-memory GPUs
- **Single GPU Support Only**: The GPU version does not support multi-GPU parallelization. However, the CPU part of the calculation can still be accelerated using MPI and OpenMP for parallel processing.



### Example GPU compilation

```bash
cd Src/
make gpu
```

This creates the `ShengBTE_gpu` executable in the root directory.

## 4th-IFCs files: FORCE_CONSTANTS_4TH

This file contains the fourth-order interatomic force constant matrix, and uses a sparse description to save space. To construct 4th-IFCs, one can refer to the `Fourthorder` python scripts. The format of this force constants is a direct extension of third-order force constants, and contains `nb` blocks of such matrix:

- A blank line
- A 1-based sequential index (from 1 to `nb`)
- A line with the Cartesian coordinates of the second unit cell in Å
- A line with the Cartesian coordinates of the third unit cell in Å
- A line with the Cartesian coordinates of the fourth unit cell in Å
- A line with the 1-based indices of the four atoms involved, each from 1 to `natoms`
- 81 lines of force constant matrix in $\frac{\textrm{eV}}{\textrm{Å}^4}$. The indexes at the beginning label the Cartesian axes.

An example block of this file:

```

1
0.0000000000e+00 0.0000000000e+00 0.0000000000e+00
0.0000000000e+00 0.0000000000e+00 0.0000000000e+00
0.0000000000e+00 0.0000000000e+00 0.0000000000e+00
     1      1      1      1
 1  1  1  1       -42.0002584218
 1  1  1  2         0.0000000000
 1  1  1  3         0.0000000000
 1  1  2  1         0.0000000000
 1  1  2  2        40.3875689003
 1  1  2  3         0.0000000000
 1  1  3  1         0.0000000000
 1  1  3  2         0.0000000000
 1  1  3  3        40.3875689003
 1  2  1  1         0.0000000000
 1  2  1  2        40.3875689003
 1  2  1  3         0.0000000000
 1  2  2  1        40.3875689003
 1  2  2  2         0.0000000000
 1  2  2  3         0.0000000000
 1  2  3  1         0.0000000000
 1  2  3  2         0.0000000000
 1  2  3  3         0.0000000000
 1  3  1  1         0.0000000000
 1  3  1  2         0.0000000000
 1  3  1  3        40.3875689003
 1  3  2  1         0.0000000000
 1  3  2  2         0.0000000000
 1  3  2  3         0.0000000000
 1  3  3  1        40.3875689003
 1  3  3  2         0.0000000000
 1  3  3  3         0.0000000000
 2  1  1  1         0.0000000000
 2  1  1  2        40.3875689003
 2  1  1  3         0.0000000000
 2  1  2  1        40.3875689003
 2  1  2  2         0.0000000000
 2  1  2  3         0.0000000000
 2  1  3  1         0.0000000000
 2  1  3  2         0.0000000000
 2  1  3  3         0.0000000000
 2  2  1  1        40.3875689003
 2  2  1  2         0.0000000000
 2  2  1  3         0.0000000000
 2  2  2  1         0.0000000000
 2  2  2  2       -42.0002584218
 2  2  2  3         0.0000000000
 2  2  3  1         0.0000000000
 2  2  3  2         0.0000000000
 2  2  3  3        40.3875689003
 2  3  1  1         0.0000000000
 2  3  1  2         0.0000000000
 2  3  1  3         0.0000000000
 2  3  2  1         0.0000000000
 2  3  2  2         0.0000000000
 2  3  2  3        40.3875689003
 2  3  3  1         0.0000000000
 2  3  3  2        40.3875689003
 2  3  3  3         0.0000000000
 3  1  1  1         0.0000000000
 3  1  1  2         0.0000000000
 3  1  1  3        40.3875689003
 3  1  2  1         0.0000000000
 3  1  2  2         0.0000000000
 3  1  2  3         0.0000000000
 3  1  3  1        40.3875689003
 3  1  3  2         0.0000000000
 3  1  3  3         0.0000000000
 3  2  1  1         0.0000000000
 3  2  1  2         0.0000000000
 3  2  1  3         0.0000000000
 3  2  2  1         0.0000000000
 3  2  2  2         0.0000000000
 3  2  2  3        40.3875689003
 3  2  3  1         0.0000000000
 3  2  3  2        40.3875689003
 3  2  3  3         0.0000000000
 3  3  1  1        40.3875689003
 3  3  1  2         0.0000000000
 3  3  1  3         0.0000000000
 3  3  2  1         0.0000000000
 3  3  2  2        40.3875689003
 3  3  2  3         0.0000000000
 3  3  3  1         0.0000000000
 3  3  3  2         0.0000000000
 3  3  3  3       -42.0002584218
```

## Settings in `CONTROL` file

This `CONTROL` file contains all the user-specified settings and parameters, including crystal structural information, broadening factor, q-mesh, temperature, functions, etc. 

### General usage

To call FourPhonon capabilities, one should add a new `&flags` name-list:

- `four_phonon` (logical, default=.false.): compute four-phonon phase space and four-phonon scattering rates.
- `four_phonon_iteration=.true.`(logical, default=.false.): compute four-phonon scattering and solve BTE exactly. This only works when `four_phonon=.true.`.

Here, we show some practical usage of this flag in combination with other flags.

- `onlyharmonic=.true.` and `four_phonon=.true.`: only compute four-phonon phase space.
- `convergence=.false.` and `four_phonon=.true.`: compute thermal conductivity at RTA level for both three- and four-phonon scatterings; this setting will output contributions from different four-phonon channels.
- `four_phonon=.true.`(`convergence` is default to be .true.): compute thermal conductivity with three-phonon iterative scheme but treat four-phonon scattering at RTA level.
- `four_phonon=.true.` and `four_phonon_iteration=.true.`: compute thermal conductivity with exact solution of BTE for both three- and four-phonon scatterings.

*Note that: all other parameters in `CONTROL` file, like temperature or q-mesh, apply to both three- and four-phonon processes. `nanowires` function is not supported in `FourPhonon` package.

### Sampling method for accelerated RTA solution

For estimating scattering rate from a sample of scattering processes, please also cite [Z. Guo *et al.*, [npj Comput. Mater. 10, 31 (2024).](https://www.nature.com/articles/s41524-024-01215-8)] and use the tags below:

- `num_sample_process_3ph` and `num_sample_process_3ph_phase_space` (int, default=`-1`): the number of sample taken from each mode for estimating three-phonon phase space and scattering rate.  `-1` means not taking any sample and performing rigorous calculation. Note that when `convergence=.true.`, i.e., using iterative scheme for three-phonon scattering calculation, `num_sample_process_3ph` must be `-1`, since the sampling method works on relaxation time approximation.
- `num_sample_process_4ph` and `num_sample_process_4ph_phase_space` (int, default=`-1`): the number of sample taken from each mode for estimating four-phonon phase space and scattering rate.  `-1` means not taking any sample and performing rigorous calculation. Note that when `four_phonon=.false.`, `num_sample_process_4ph` and `num_sample_process_4ph_phase_space` must both be `-1`.
- Example: compute 3ph scattering rates exactly (iterative), and compute 4ph scattering at RTA level with 100000 sampled processes
    
    ```fortran
    &parameters
    	T=300
      scalebroad=1
      num_sample_process_3ph_phase_space = -1
      num_sample_process_3ph = -1
      num_sample_process_4ph_phase_space = 100000
      num_sample_process_4ph = 100000         
    &end
    ```
    

### Input of TDEP IFCs

For harmonic force constants in [TDEP](https://github.com/tdep-developers/tdep) format, use the tags/scripts described below:

- `tdep=.true.` activates the processing of TDEP format **harmonic** IFCs and this input file name should always be `infile.forceconstant`. If the system is a polar material, apart from the normal parameters like Born effective charges and dielectric constants, we also need users to manually input coupling parameter $\lambda$  `Ewald` that is written out near the end of TDEP IFCs like this one:
    
    ```fortran
    0.675614929199219       # Coupling parameter in Ewald summation
    ```
    
- Attention: the input harmonic IFCs are bare IFCs **without** polar correction so that `FourPhonon` can process it and add correction internally. The module has largely replied on TDEP source codes and we cannot guarantee the success. We advise users to check processed frequencies and 3ph scattering rates.
- For non-polar materials, we have also developed a script to convert TDEP format to Phonopy format, `tdep2phonopy.f90`. This way you do not need to activate `tdep` flag. After Fortran compilation, run the executable (this is an example):
    
    ```bash
    ./tdep2phonopy outfile.forceconstant 10 10 1 2
    ```
    
    this converts TDEP IFCs (`outfile.forceconstant`) to Phonopy format with a $10\times10\times1$ supercell constructed from a two-atom unitcell. Check the converted phonon dispersion. We suggest using a large enough supercell size to store all pair interactions in 2nd-IFCs.
    
- Example: compute 3ph scattering rates exactly (iterative),and compute 4ph scattering at RTA level. A TDEP format harmonic IFCs is used and non-analytical correction is applied with Ewald parameter $\lambda=0.675614929199219$. Born effective charges are also needed (not shown here), similar to what `ShengBTE` does for polar systems.
    
    ```fortran
    &parameters
    	T=300
      scalebroad=0.1
    	Ewald=0.675614929199219
    &end
    &flags
    	tdep=.TRUE.
      nonanalytic=.TRUE. ! this tag now needs to be specified explicitly
    	convergence=.TRUE.
    	four_phonon=.TRUE.
    &end
    ```
    

For 3rd/4th-IFCs in TDEP format, we developed a standalone script to convert them to ShengBTE format, `tdep2ShengBTE.f90`. We need `infile.ucposcar` in the same directory. After Fortran compilation, run the executable:

- Example:
    
    ```bash
    ./tdep2ShengBTE 3 2 outfile.forceconstant_thirdorder
    ```
    
    this converts TDEP 3rd-IFCs (`outfile.forceconstant_thirdorder`) of a two-atom unitcell to ShengBTE format.
    

## Output files

Besides the routine output files from previous `ShengBTE` program, `FourPhonon` generates these output files:

- `BTE.Numprocess_4ph`: number of allowed four-phonon scattering processes, for each irreducible q point and phonon band
- `BTE.P4`: phase space available for four-phonon processes, for each irreducible q point and phonon band
- `BTE.P4_total`: total volume in phase space available for four-phonon processes
- `BTE.P4_plusplus*`, `BTE.P4_plusminus*`, `BTE.P4_minusminus*`: similar to BTE.P4 but only includes contributions from ++/+-/- - processes

*Note for four-phonon scatterings, there are three different channels: recombination (++), redistribution (+-) and splitting (- -) processes.

Under temperature-dependent directories: (the unit of scattering rates is $\rm ps^{-1}$ if not specified)

- `BTE.WP4`: weighted phase space available for four-phonon processes ([rad/ps]$^{-5}$, 2nd column) vs angular frequency (rad/ps, 1st column) for each irreducible q point and phonon band
- `BTE.WP4_plusplus*`, `BTE.WP4_plusminus*`, `BTE.WP4_minusminus*`: similar to BTE.WP4 but only includes contributions from ++/+-/- - processes
- `BTE.w_3ph`: three-phonon scattering rates for each irreducible q point and phonon band, this file replaces the original output `BTE.w_anharmonic`. Absorption and emission processes are written out into `BTE.w_3ph_plus` and `BTE.w_3ph_minus`
- `BTE.w_4ph`: four-phonon scattering rates for each irreducible q point and phonon band. Similarly, we provide the contributions from different channels: `BTE.w_4ph_plusplus`, `BTE.w_4ph_plusminus` and `BTE.w_4ph_minusminus`
- `BTE.kappa*`: thermal conductivity results are written out as usual
