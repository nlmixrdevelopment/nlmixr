
# nlmixr: an R package for population PKPD modeling
***
![nlmixr](logo.png)

##### Authors: Matthew Fidler, Yuan Xiong, Rik Schoemaker, Justin Wilkins, Teun Post, Richard Hooijmaijers, Wenping Wang

[![R build status](https://github.com/nlmixrdevelopment/nlmixr/workflows/R-CMD-check/badge.svg)](https://github.com/nlmixrdevelopment/nlmixr/actions)

***

`nlmixr` is an R package for fitting general dynamic models,
pharmacokinetic (PK) models and pharmacokinetic-pharmacodynamic (PKPD)
models in particular, with either individual data or population
data. The nlme and SAEM estimation routines can be accessed using a
universal user interface (UUI), that provides universal model and
parameter definition syntax and results in a fit object that can be
used as input into the `Xpose` package. Running nlmixr using the UUI
is described in the vignette:
ihttps://github.com/nlmixrdevelopment/nlmixr/blob/master/vignettes/running_nlmixr.Rmd

Under the hood `nlmixr` has five main modules:  

1. `dynmodel()` and its mcmc cousin `dynmodel.mcmc()` for nonlinear
   dynamic models of individual data;
2. `nlme_lin_cmpt()`for one to three linear compartment models of
   population data with first order absorption, or i.v. bolus, or
   i.v. infusion using the nlme algorithm;
3. `nlme_ode()` for general dynamic models defined by ordinary
   differential equations (ODEs) of population data using the nlme
   algorithm;
4. `saem_fit` for general dynamic models defined by ordinary differential equations (ODEs) of population data by the Stochastic Approximation Expectation-Maximization (SAEM) algorithm;  
5. `gnlmm` for generalized non-linear mixed-models (possibly defined
   by ordinary differential equations) of population data by the
   adaptive Gaussian quadrature algorithm.

A few utilities to facilitate population model building are also included in `nlmixr`.

For a brief vignette describing the modules, please see:
https://github.com/nlmixrdevelopment/nlmixr/blob/master/inst/nlmixr-intro.pdf

The examples in the vignette can be run using VignetteDemo.R and the associated data files available at:
https://github.com/nlmixrdevelopment/nlmixr/tree/master/vignettes

For PKPD modeling (with ODE and dosing history) with
[Stan](http://mc-stan.org/), check out Yuan's package PMXStan:
https://github.com/yxiong1/pmxstan.

# Using a Docker Image for running nlmixr

One of the easiest way to setup nlmixr is to docker image.  For more details see:

https://github.com/RichardHooijmaijers/nlmixr.docker


# Windows installer
For those not interested in customized installation on Windows, we
 **recommend** you download a Windows installer for your platform from
 the following link:
 https://github.com/nlmixrdevelopment/nlmixr/releases/
 
 
 **STOP AT THIS POINT IF YOU USED THE INSTALLER**
 
 Sometimes trying to install yourself will corrupt the installation and nlmixr will not run.

# Installation on Windows
To replicate the environment that was used in Windows for `nlmixr`
development, you will need administrator rights, and you should
perform the following steps:

1. Install R 3.5.0 (or later) from the R website.
   - For best results, we suggest you use `C:\R\R-3.5.0`, but
     you can also use the default location (`C:\Program Files\R\R-3.5.0`) as well, if really needed.
   - For 64-bit Windows, it is best practice to include *only* the
     64-bit version.  If you include 32-bit files, some packages may
     not run correctly.  Additionally, both the 32- and 64-bit
     binaries have to be compiled for every package. Similarly, if on
     32-bit Windows, install only the 32-bit version of R (and Python,
     and Rtools).
	 
2. Install the appropriate version of Rtools for Windows, currently version 3.5, from [here](https://cran.r-project.org/bin/windows/Rtools/).
   - This is an absolute requirement, since it includes C++ and related compilers not usually available under Windows.
   - For best results, use the default location of `c:\Rtools`
     - `RxODE`, a required component of `nlmixr`, checks and sets up the path based on the following:
	    a. `Rtools` is in the path (fastest and recommended option)
		b. `Rtools` was installed with information saved to the Windows registry, and `RxODE` can 
		   find the installation.
		c. `Rtools` is on a hard drive installed in either `Rtools` or `RBuildTools`
     - If you are on 64-bit windows, please *do not install* the R
       3.3.x 32-bit toolchain.  These files can interfere with some
       packages that compile binaries, with unpredictable
       consequences.  Similarly, only install 32-bit Rtools on 32-bit
       versions of Windows.
   - Make sure the compilers have been added to the Windows `PATH`
     environment variable, or `RxODE` and `nlmixr` will not work (this
     should be done automatically during installation).
3. Install a version of Python for Windows.
   - This is used for its symbolic algebra package [SymPy](http://sympy.org/).
   - If you already have a python installed, it is much easier to
     piggy-back on that system than install two different python
     versions.
   - If you are new to python, a very robust Python distribution that
     includes [SymPy](http://sympy.org/) and many packages that may be
     useful to the data scientist and/or pharmacometrician is
     [Anaconda](https://www.anaconda.com/download/). Although very
     straightforward and easy to install, it is quite a large download
     and contains much more than you will need to run `nlmixr`. When
     installing, use the Python 3.6 version. During the installation,
     Anaconda provides the option of adding itself to the `PATH`
     environment variable, but advises against it; please do this
     anyway (despite the red warning).
   - Another option is to use [official Python](http://python.org),
     although you will need to install [SymPy](http://sympy.org/)
     separately if you go this route, which is sometimes not
     straightforward under Windows 10 owing to folder permissions (see
     [here](https://stackoverflow.com/questions/31172719/pip-install-access-denied-on-windows)
     for a few workarounds). Nonetheless, see
     [here](http://simpy.readthedocs.io/en/latest/simpy_intro/installation.html)
     for instructions for installation from source or using
     `pip`. 
   - Regardless of the option you choose, please use like with like
     (64-bit Python for 64-bit Windows, for example).
   - Once again, make sure Python has been added to the Windows `PATH`
     environment variable, or `RxODE` and `nlmixr` will not work, no
     matter what Anaconda might say.
3. Install `devtools`.
   - This package is required to install packages from Github, amongst other things.
   - This can be done from a clean R session by `install.packages("devtools")`.
4. Load `devtools` using `library(devtools)`.
5. Install `RxODE`.
   - Currently the new version of `RxODE` is in the process of being
     uploaded to CRAN.  `nlmixr` needs this newer version of `RxODE` to
     function correctly. To install this version, use the command:
     `install_github("nlmixrdevelopment/RxODE")`.
   - Once installed, type `RxODE::rxWinPythonSetup()` to install the
     required package `SnakeCharmR`and to make sure Python and SymPy
     are working properly.  Since `SnakeCharmR` is orphaned, it may
     not work without installing directly from github with the
     following command:
     `devtools::install_github("nlmixrdevelopment/SnakeCharmR")`
   - Restart your R session.
   - As a quick test, you can make sure that R and Python can
     communicate by typing the command `library(SnakeCharmR)`.
   - To validate or test the installation of `RxODE` completely, you
     can type the following `library(RxODE); rxTest();` and it will
     run all of the unit tests in RxODE to make sure it is running
     correctly on your system. (Note that the `testthat` package is
     required for this, and it will take a long time.).  This also
     requires installing `RxODE` with the tests.  Devtools changed so
     that tests are no longer installed by default.  The correct
     command to install the test suite is:
	 
```{r}
devtools::install_github("nlmixrdevelopment/RxODE", INSTALL_opts="--install-tests")
```
	 
6. Install `nlmixr`.
   - Load `devtools` again using `library(devtools)`
   - Since `nlmixr` shares some stan headers, nlmixr now needs `C++14`
     support.  To ease the setup you can run the command
     `RxODE::rxC14()` which creates the following build file in
     `~/.R/Makevars.win`.  If there are problems you may have to edit
     the file manually to match your configuration.
	 
```{make}
# Makvars.win
CXX14=$(BINPREF)g++ $(M_ARCH)
CXX14STD=-std=c++1y
CXX14FLAGS=-O2 -Wall
```

   - Finally, install `nlmixr` by running
     `devtools::install_github("nlmixrdevelopment/nlmixr")`

# Installation on Linux
Instructions for Ubuntu-alike distributions are given here
(specifically, [Ubuntu 16.04 Xenial
Xerus](http://releases.ubuntu.com/16.04/)), but all current Linux
distributions are supported, in principle.

1. Install R 3.4.1 (or later) from an appropriate repository (Ubuntu Xenial shown below, based on instructions provided [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2)).
   - You will need administrator privileges (i.e. access to `sudo`). Provide your admin password when asked. 
   - Add the official CRAN repository for Ubuntu: 
     `sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9`
   - Add the official CRAN repository for Ubuntu: 
     `sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'`. If you aren't using Ubuntu Xenial, change `xenial` to match your distribution's codename.
   - Now refresh the package list using `sudo apt-get update`. 
   - We can now install base R and required development libraries, and their dependencies: `sudo apt-get install r-base r-base-dev libssl-dev` 
2. Install Python dependencies.
   - Enter this: `sudo apt-get install python-sympy python-pip python-setuptools python3-pip python-dev python3-dev`.
3. Install `devtools` and dependencies.
   - This package is required to install packages from Github, amongst other things.
   - Some Linux distributions don't include build tools out of the box. To be safe, check this: `sudo apt-get install build-essential`
   - Install `devtools` from a clean R session by entering `install.packages("devtools")`.
4. In R, load `devtools` using `library(devtools)`.
5. Install `RxODE`.
   - Currently the new version of `RxODE` is in the process of being
     uploaded to CRAN.  `nlmixr` needs this newer version of `RxODE` to
     function correctly. To install this version, use the command:
     `install_github("nlmixrdevelopment/RxODE")`.
   - Install `SnakeCharmR` using `install_github("nlmixrdevelopment/SnakeCharmR")`.
   - Restart your R session.
   - As a quick test, you can make sure that R and Python can
     communicate by typing the command `library(SnakeCharmR)`.
   - To validate or test the installation of `RxODE` completely, you
     can type the following `library(RxODE); rxTest();` and it will
     run all of the unit tests in RxODE to make sure it is running
     correctly on your system. (Note that the `testthat` package is
     required for this, and it will take a long time.)
	 - Note this requires a special install option to make sure the
       tests are installed:  `remotes::install_github("nlmixrdevelopment/RxODE", INSTALL_opts="--install-tests")`
6. Install `nlmixr`.
   - This can be done by `install_github("nlmixrdevelopment/nlmixr")`
 
# Installation on macOS

Instructions for macOS 10.13 Mojave are provided here. They should be
broadly extensible to all recent releases of macOS, however.

1. Install R 3.6.1 from the R website.
   - Download and install `R-3.6.1.pkg` (or later) from [CRAN](http://cran.r-project.org/bin/macosx/).
2. Install Python dependencies.
   - Install `pip` from the macOS terminal prompt: `sudo easy_install pip`.
     - If `easy_install` doesn't work, you can download by `curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py` and then run `sudo python get-pip.py`
   - Install `sympy` using `pip`: `sudo -H pip install sympy`.
3. Install build tools.
   - Install the clang 8.0+ and gfortran 6.1+ packages from [CRAN development
     tools website](http://cran.r-project.org/bin/macosx/tools/)
   - Install Mac OSX command like tools by typing `xcode-select
     --install` and then choose `install`.  This allows you to compile
     the package without downloading the entire (large) Xcode
     package. Make sure to perform this step even if the Xcode package
     is already installed.
   - If you are using macOS >= 10.14 you have to also type the
     following in the terminal to install the compiler headers:

```sh
sudo installer -pkg \
/Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg \
-target /
```
   - Change `~/.Rprofile` to have the following line to adjust the
     paths in R to use the CRAN compilers (if not present, create the
     .Rprofile file):
   
```R
# make sure to adjust the version number of clang based on the downloaded version in step 3 (e.g. clang8/bin)
Sys.setenv("PATH"=paste0("/usr/local/clang8/bin:/usr/local/gfortran/bin:", Sys.getenv("PATH")))
```

**Make Sure you have an extra line at the end of this file; Otherwise R will not pick up the correct paths**.

You also need to create a directory `~/.R` and then put a file `Makevars` in it;  The context of the file should contain:

```sh
# ~/.R/Makevars
CXX14FLAGS=-O2 -march=native -mtune=native -arch x86_64 -ftemplate-depth-256 -fPIC
```

4. Install `devtools` and dependencies.
   - This package is required to install packages from Github, amongst other things.
   - Install `devtools` from a clean R session by entering `install.packages("devtools")`.
5. In R, load `devtools` using `library(devtools)`.
6. Install `RxODE`.
   - Currently the new version of `RxODE` is in the process of being
     uploaded to CRAN.  `nlmixr` needs this newer version of `RxODE` to
     function correctly. To install this version, use the command:
     `install_github("nlmixrdevelopment/RxODE")`.
   - Install `reticulate` using `install.packages("reticulate")`.
   - Restart your R session.
   - To validate or test the installation of `RxODE` completely, you
     can type the following `library(RxODE); rxTest();` and it will
     run all of the unit tests in RxODE to make sure it is running
     correctly on your system. (Note that the `testthat` package is
     required for this, and it will take a long time.)
	 - Note this requires a special install option to make sure the
       tests are installed:  `remotes::install_github("nlmixrdevelopment/RxODE", INSTALL_opts="--install-tests")`
7. Install `nlmixr`.
   - This can be done by `install_github("nlmixrdevelopment/nlmixr")`
   
# Notes on Python configuration

If you have multiple versions of python installed on your system, you
may run into some more issues.  To be more careful in your install,
you can install SnakeCharmR as follows:

```R
Sys.setenv(SNAKECHARMR_PYTHON_VERSION="/path/to/python/with/sympy/installed")
devtools::install_github("nlmixrdevelopment/SnakeCharmR")
```

You may also have to adjust the following python variables to match
the version with sympy installed:

Environmental Variable |  Correct value
--------------------------|---------------------------
`PYTHON_EXE`  | Path where the python with sympy is installed
`PYTHONHOME` | Path where the python with sympy is installed
`PYTHON_INCLUDE` | Path where the python libraries are installed; In windows this is `PYTHONHOME\include`
`PYTHON_LIB` | Path where python libraries are installed; In windows this is `PYTHONHOME\libs`
`PYTHONPATH` | Path where python searches.  In windows this is a path-style variable including `PYTHONHOME\DLLs`, `PYTHONHOME\Lib` and `PYTHONHOME\Lib\site-packages`.  
`PYTHONSTARTUP` | In windows, this value is unset if present

I believe you could also unset some of these variables and python can
figure them out.  Please be warned if these variables are not setup
correctly a call to python will abort python.  The python abort does
not respect R's way of doing things and will also abort/crash R.

To test without running a model you may wish to try

```R
library(SnakeCharmR)
```
If it doesn't crash R, then python is likely setup correctly.

After that point you can then type:
```R
library(RxODE)
rxSymPyVersion()
```
It will show the version of sympy that you are using if SnakeCharmR is setup correctly.

# Testing the install

Once nlmixr is installed, you can test if everything is working well (as well as running shinyMixR)
by [running this install checker](https://raw.githubusercontent.com/nlmixrdevelopment/nlmixr/master/build/test_install.R) by Richard Hooijmaijers & Matt Fidler.

Once downloaded you would source the code;  If successful you should see something similar to the following:
```
> source("test_install.R")
Correct R version: Yes, R version 3.5.0 (2018-04-23)
RxODE installed: Yes
Python installed: Yes, Python 3.6.2
sympy installed: Yes
devtools package installed: Yes
Rtools installed: Yes
RxODE package installed: Yes
SnakeCharmR package installed: Yes
nlmixr package installed: Yes
xpose.nlmixr package installed: Yes
shinyMixR package installed: Yes
Loading required package: nlme
nlmixr run under nlme: Yes
nlmixr run under saem: Yes
---- Installation test finished! ----
```


