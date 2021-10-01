# `FourPhonon`: An extension module to `ShengBTE` for computing four-phonon scattering rates and thermal conductivity

## Authors and references for `FourPhonon`:

- Zherui Han [zrhan@purdue.edu](mailto:zrhan@purdue.edu)
- Xiaolong Yang [xiaolongyang1990@gmail.com](mailto:xiaolongyang1990@gmail.com)
- Wu Li [wu.li.phys2011@gmail.com](mailto:wu.li.phys2011@gmail.com)
- Tianli Feng [Tianli.Feng2011@gmail.com](mailto:Tianli.Feng2011@gmail.com)
- Xiulin Ruan [ruan@purdue.edu](mailto:ruan@purdue.edu)

References: please refer to our [GitHub homepage](https://github.com/FourPhonon)

1. Feng, T. & Ruan, X. Quantum mechanical prediction of four-phonon scattering rates and reduced thermal conductivity of solids. Phys Rev B **93**, 045202 (2016).
2. Feng, T., Lindsay, L. & Ruan, X. Four-phonon scattering significantly reduces intrinsic thermal conductivity of solids. Phys Rev B **96**, 161201 (2017).
3. Han, Z., Yang, X., Li, W., Feng, T. & Ruan, X. FourPhonon: An extension module to ShengBTE for computing four-phonon scattering rates and thermal conductivity. Comput Phys Commun **270** (2022) 108179, https://doi.org/10.1016/j.cpc.2021.108179.

## The original authors of `ShengBTE` and references:

- Wu Li [wu.li.phys2011@gmail.com](mailto:wu.li.phys2011@gmail.com)
- Jesús Carrete Montaña [jcarrete@gmail.com](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:jcarrete@gmail.com)
- Nebil A. Katcho [nebil.katcho@gmail.com](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:nebil.katcho@gmail.com)
- Natalio Mingo [natalio.mingo@cea.fr](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:natalio.mingo@cea.fr)

References: please refer to the ShengBTE [link](http://www.shengbte.org/how-to-cite).

## How to download and compile `FourPhonon`

`FourPhonon` is built within `ShengBTE` and updates `ShengBTE` to a new version. The codes are hosted at GitHub, and you can download the latest distribution from this repository: [https://github.com/FourPhonon](https://github.com/FourPhonon) (you can also find this link from `ShengBTE` [website](http://www.shengbte.org/home)). The compilation of this new version is the same as the previous `ShengBTE`: after setting proper paths in `arch.make`, one can then run `make` in the `Src` subdirectory. An executable `ShengBTE` will appear in the root directory of this distribution.

## How to call `FourPhonon` function

Besides the routine inputs and `CONTROL` file of `ShengBTE`, `FourPhonon` requires a fourth-order force constants and some new flags in `CONTROL` file.

### 4th-IFCs files: FORCE_CONSTANTS_4TH

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

### Flags in `CONTROL` file

This file contains all the user-specified settings and parameters, including crystal structural information, broadening factor, q-mesh, temperature, etc. To call FourPhonon capabilities, one should add a new `&flags` namelist:

- `four_phonon` (logical, default=.false.): compute four-phonon phase space and four-phonon scattering rates.

Here, we show some practical usage of this flag in combination with other existing flags.

- `onlyharmonic=.true.` and `four_phonon=.true.`: only compute four-phonon phase space
- `convergence=.false.` and `four_phonon=.true.`: compute thermal conductivity at RTA level for both three- and four-phonon scatterings
- `four_phonon=.true.`(`convergence` is default to be .true.): compute thermal conductivity with three-phonon iterative scheme but treat four-phonon scatterings at RTA level

*Note that: all other parameters in `CONTROL` file, like temperature or q-mesh, apply to both three- and four-phonon processes. `nanowires` function is not supported in `FourPhonon` package.

## Output files

Besides the routine output files from previous `ShengBTE` program, `FourPhonon` generates these output files:

- `BTE.Numprocess_4ph`: number of allowed four-phonon scattering processes, for each irreducible q point and phonon band
- `BTE.P4`: phase space available for four-phonon processes, for each irreducible q point and phonon band
- `BTE.P4_total`: total volume in phase space available for four-phonon processes
- `BTE.P4_plusplus*`, `BTE.P4_plusminus*`, `BTE.P4_minusminus*`: similar to BTE.P4 but only includes contributions from ++/+-/- - processes

*Note for four-phonon scatterings, there are three different channels: recombination (++), redistribution (+-) and splitting (- -) processes.

Under temperature-dependent directories: (the unit of scattering rates is ps$^{-1}$ if not specified)

- `BTE.WP4`: weighted phase space available for four-phonon processes ([rad/ps]$^{-5}$, 2nd column) vs angular frequency (rad/ps, 1st column) for each irreducible q point and phonon band
- `BTE.WP4_plusplus*`, `BTE.WP4_plusminus*`, `BTE.WP4_minusminus*`: similar to BTE.WP4 but only includes contributions from ++/+-/- - processes
- `BTE.w_3ph`: three-phonon scattering rates for each irreducible q point and phonon band, this file replaces the original output `BTE.w_anharmonic`. Absorption and emission processes are written out into `BTE.w_3ph_plus` and `BTE.w_3ph_minus`
- `BTE.w_4ph`: four-phonon scattering rates for each irreducible q point and phonon band. Similarly, we provide the contributions from different channels: `BTE.w_4ph_plusplus`, `BTE.w_4ph_plusminus` and `BTE.w_4ph_minusminus`
- `BTE.w_4ph_normal`: four-phonon scattering rates of normal processes, for each irreducible q point and phonon band. This file has 5 columns and each one represents: angular frequency in rad/ps, recombination channel, redistribution channel, splitting channel, overall scattering rates from normal processes
- `BTE.w_4ph_Umklapp`: four-phonon scattering rates of Umklapp processes, for each irreducible q point and phonon band. The format is the same as `BTE.w_4ph_normal`
- `BTE.kappa*`: thermal conductivity results are written out as usual
