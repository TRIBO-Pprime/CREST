Table of Contents
-----------------

- [Brief description](#brief-description)
- [Building CREST](#building-crest)
- [Dependencies](#dependencies)
- [Third party components](#third-party-components)
- [Typical use](#typical-use)
- [Full description](#full-description)
- [Documentation](#documentation)
- [License](#license)

Brief description
-----------------

**CREST** - **C**omputational **R**esource for **E**roded **S**urface **T**opology - is designed to generate numerical rough surfaces that respect:

* size: \(n \times m\) regular grid
* statistical moments: Skewness \(Ssk\) and Kurtosis \(Sku\)
* principal correlation lengths: \(\tau_1\) and \(\tau_2\)
* asperity orientation: \(\alpha\) angle
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

* *TOOLIB* ... Some general tools like FFT, file handling, minimization, sorting, etc.

* *TPGLIB* ... Some more specific programs like filtering, anisotropy analysis, derivation, etc.

As a consequence the three packages have to be downloaded.

[top](#table-of-contents)

Third party components
----------------------

TOOLIB also uses external codes such as:

+ *FFTW3.3* ... [Fastest Fourier Transform in the West](https://www.fftw.org/) ... GNU General Public License

+ *Pikaia_oop* ... [Modern Fortran Edition of the Pikaia Genetic Algorithm](http://github.com/jacobwilliams/pikaia) ... BSD like

+ *GNUFOR* ... [Gnuplot Fortran](https://people.math.sc.edu/Burkardt/f_src/gnufor/gnufor.html) ... GNU General Public License

+ *Bspline-fortran* ... [Multidimensional B-Spline Interpolation of Data on a Regular Grid](https://github.com/jacobwilliams/bspline-fortran) ... BSD like

`.sur` surface files can be visualized and analyzed with [Gwyddion software](http://gwyddion.net/download.php), a modular program for SPM (scanning probe microscopy) data visualization and analysis.

Typical use
-----------

The program `main` reads a script file `my_script.md` where the following parameters are defined:

+ image size (pix)
+ surface size (m)
+ periodicity
+ correlation lengths
+ roughness orientation
+ number of available threads
+ statistical parameters
+ output surface name
+ ...

as well as the algorithm for rough surface generation.

Run:

```bash
./main cfg/my_script.md
```


[top](#table-of-contents)

Full description
----------------

## Test 1

| size *m* | size *n* | width (m) | height (m) | periodic? | \(\tau_1\) | \( \tau_2 \) | Ssk | Sku |
|:--------:|:--------:|:---------:|:----------:|:---------:|:----------:|:------------:|:---:|:---:|
| 1024     | 512      | 200.e-6   | 100.e-6    | **True**  | 30.e-6     | 10.e-6       | -3. | 15. |

Run:

```bash
./main cfg/test01.md
```

Result:
[test01_img](./media/test01.jpg)

## Test 2

| size *m* | size *n* | width (m) | height (m) | periodic? | \(\tau_1\) | \( \tau_2 \) | Ssk | Sku |
|:--------:|:--------:|:---------:|:----------:|:---------:|:----------:|:------------:|:---:|:---:|
| 1024     | 512      | 200.e-6   | 100.e-6    | **False** | 30.e-6     | 10.e-6       | -3. | 15. |

Run:

```bash
./main cfg/test02.md
```

Result:
[test02_img](./media/test02.jpg)

[top](#table-of-contents)

Documentation
-------------
The documentation is automatically generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford), an automatic documentation generator for modern Fortran programs.

[top](#table-of-contents)

License
-------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warrenty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

You should have received a [copy](https://github.com/TRIBO-Pprime/CREST/LICENSE) of the GNU General Public License along with this program. If not, see the [GNU website](https://www.gnu.org/licenses/gpl.html).
