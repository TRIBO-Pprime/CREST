Table of Contents
-----------------

- [Brief description](#brief-description)
- [Building CREST](#building-crest)
- [Dependencies](#dependencies)
- [Third party components](#third-party-components)
- [Typical use](#typical-use)
- [Full description](#full-description)
- [License](#license)

Brief description
-----------------

**CREST** - **C**omputational **R**esource for **E**roded **S**urface **T**opology - is designed to generate numerical rough surfaces that respect:

* size: $n \times m$ regular grid
* statistical moments: Skewness $Ssk$ and Kurtosis $Sku$
* principal correlation lengths: $\tau_1$ and $\tau_2$
* asperity orientation: $\alpha$ angle
* periodicity: with/without

The programs are written in recent Fortran (2003+)

[top](#table-of-contents)

Building CREST
--------------

Make sure that the dependencies described below are built.

Debug mode:

```bash
make debug
```
Normal mode:

```bash
make
```

Dependencies
------------

CREST needs some components that are available through the following packages:

* *TOOLIB*

  Some general tools like FFT, file handling, minimization, sorting, etc.

* *TPGLIB*

  Some more specific programs like filtering, anisotropy analysis, derivation, etc.

As a consequence the three packages have to be downloaded.

[top](#table-of-contents)

Third party components
----------------------

TOOLIB also uses external codes such as:

+ *FFTW3.3*, [Fastest Fourier Transform in the West](https://www.fftw.org/)

  GNU General Public License

+ *Pikaia_oop*, [Modern Fortran Edition of the Pikaia Genetic Algorithm](http://github.com/jacobwilliams/pikaia)

  BSD like

+ *GNUFOR*, [Gnuplot Fortran](https://people.math.sc.edu/Burkardt/f_src/gnufor/gnufor.html)

  GNU General Public License

+ *Bspline-fortran*, [Multidimensional B-Spline Interpolation of Data on a Regular Grid](https://github.com/jacobwilliams/bspline-fortran)

  BSD like

Typical use
-----------


[top](#table-of-contents)

Full description
----------------



License
-------

[^1]:
    CeCILL-C license
[^2]:
    GNU GPL
[^3]:
    New BSD-3
[^4]:
    HSL software is strictly intended for Personal academic use on the download page
