# `FourPhonon`: An extension module to `ShengBTE` for computing four-phonon scattering rates and thermal conductivity

`FourPhonon` is a Boltzmann transport equation and phonon scattering solver, which incorporates four-phonon scattering formalism developed by [Ruan group](https://engineering.purdue.edu/NANOENERGY/) (Dr. Tianli Feng and Dr. Xiulin Ruan). This program, based on `ShengBTE` platform, can compute four-phonon scattering rates in crystals, give exact solution of linearized phonon BTE and the resulted thermal conductivity. More details can be found in our [paper](https://doi.org/10.1016/j.cpc.2021.108179) and [manual](https://github.com/FourPhonon/FourPhonon/blob/main/Manual.md). This software is distributed under GPL-3.0 license.

## How to download and compile `FourPhonon`

`FourPhonon` is built within `ShengBTE` and has its standalone development. The codes are hosted at GitHub, and you can download the latest distribution from this repository: https://github.com/FourPhonon or 

```bash
git clone https://github.com/FourPhonon/FourPhonon.git
```

The compilation of `FourPhonon` is the same as the previous `ShengBTE`: after setting proper paths in `arch.make`, one can then run `make` in the `Src` subdirectory. An executable `ShengBTE` will appear in the root directory of this distribution.

## Authors and references for `FourPhonon`:

- Zherui Han [zrhan@purdue.edu](mailto:zrhan@purdue.edu); Xiaolong Yang [xiaolongyang1990@gmail.com](mailto:xiaolongyang1990@gmail.com); Wu Li [wu.li.phys2011@gmail.com](mailto:wu.li.phys2011@gmail.com); Tianli Feng [Tianli.Feng2011@gmail.com](mailto:Tianli.Feng2011@gmail.com); Xiulin Ruan [ruan@purdue.edu](mailto:ruan@purdue.edu); Ziqi Guo [gziqi@purdue.edu](mailto:gziqi@purdue.edu); Guang Lin [guanglin@purdue.edu](mailto:guanglin@purdue.edu); Wenjiang Zhou [wjzhou@stu.pku.edu.cn](mailto:wjzhou@stu.pku.edu.cn); Abdulaziz Alkandari [aalkanda@purdue.edu](mailto:aalkanda@purdue.edu)

**References:** 

1. T. Feng and X. Ruan, [Phys. Rev. B **93**, 045202 (2016).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.045202)
2. T. Feng, L. Lindsay, and X. Ruan, [Phys. Rev. B **96**, 161201 (2017).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.161201)
3. Z. Han *et al.*, [Comput. Phys. Commun. **270**, 108179 (2022).](https://doi.org/10.1016/j.cpc.2021.108179)
4. Z. Guo *et al.*, [npj Comput. Mater. 10, 31 (2024).](https://www.nature.com/articles/s41524-024-01215-8)

## The original authors of `ShengBTE` and references:

- Wu Li [wu.li.phys2011@gmail.com](mailto:wu.li.phys2011@gmail.com); Jesús Carrete Montaña [jcarrete@gmail.com](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:jcarrete@gmail.com); Nebil A. Katcho [nebil.katcho@gmail.com](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:nebil.katcho@gmail.com); Natalio Mingo [natalio.mingo@cea.fr](https://bitbucket.org/sousaw/shengbte/src/d1ecffd7d8c43f1746eae47e8ddf2cf537bdcf01/mailto:natalio.mingo@cea.fr)

**References:** W. Li *et al.*, Comput Phys Commun **185**, 1747 (2014). Please refer to the ShengBTE [website](http://www.shengbte.org/how-to-cite).

## Other acknowledgement

The iterative solver implemented is based/inspired by [M. Omini and A. Sparavigna, Phys B Condens Matter **212**, 101 (1995)] and [G. Fugallo *et al.*, Phys. Rev. B **88**, 045430 (2012)]. The TDEP interface is based on its source code, referring to its original paper [O. Hellman *et al.*, Phys. Rev. B **87**, 104111 (2013)].

We thank the following scholars for their comments and help during the development of this tool:

- Prof. Jesús Carrete; Prof. Te-Huan Liu

Contributions from third-party are welcomed! (Submit new branch request/pull request)

We acknowledge the [NSF CSSI Elements](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2311848&HistoricalAwards=false) program for its support (award # 2311848).
