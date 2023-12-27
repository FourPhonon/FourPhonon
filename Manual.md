# Manual


Besides the routine inputs in `CONTROL` file of `ShengBTE`, `FourPhonon` requires a fourth-order force constants and some new name-lists in `CONTROL` file. Check [ShengBTE website](https://bitbucket.org/sousaw/shengbte/src/master/README.md) for definition of other name-lists.

## Parallel environment

Version 1.1 and 1.0 were written for MPI parallelism. Starting from Version 1.2 that supports iterative solver for four-phonon scattering, we have migrated to **OpenMP** to handle large memory required for this iterative solver. In the future we may support MPI+OpenMP hybrid parallelism to allow more flexibility. Make sure to add `-qopenmp` in compilation or:

```makefile
export FFLAGS=-qopenmp -traceback -debug -O2 -static_intel
```

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

For estimating scattering rate from a sample of scattering processes, please also cite [Z. Guo *et al.*, [arXiv:2311.12935 (2023)](https://doi.org/10.48550/arXiv.2311.12935)] and use the tags below:

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
- `BTE.w_4ph_normal`: four-phonon scattering rates of normal processes, for each irreducible q point and phonon band. This file has 5 columns and each one represents: angular frequency in rad/ps, recombination channel, redistribution channel, splitting channel, overall scattering rates from normal processes
- `BTE.w_4ph_Umklapp`: four-phonon scattering rates of Umklapp processes, for each irreducible q point and phonon band. The format is the same as `BTE.w_4ph_normal`
- `BTE.kappa*`: thermal conductivity results are written out as usual
