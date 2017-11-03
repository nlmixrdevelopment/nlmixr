---
output: html_document
---
![alt tag](https://github.com/nlmixrdevelopment/nlmixr/blob/master/logo.png)

# nlmixr: an R package for population PKPD modeling
***  

##### Authors: Matthew Fidler, Yuan Xiong, Rik Schoemaker, Justin Wilkins, Mirjam Trame, Wenping Wang

***
`nlmixr` is an R package for fitting general dynamic models, pharmacokinetic (PK) models and pharmacokinetic-pharmacodynamic (PKPD) models in particular, with either individual data or population data. `nlmixr` has five main modules:  

1. `dynmodel()` and its mcmc cousin `dynmodel.mcmc()` for nonlinear dynamic models of individual data; 
2. `nlme_lin_cmpt()`for one to three linear compartment models of population data with first order absorption, or i.v. bolus, or i.v. infusion using the nlme algorithm; 
3. `nlme_ode()` for general dynamic models defined by ordinary differential equations (ODEs) of population data using the nlme algorithm; 
4. `saem_fit` for general dynamic models defined by ordinary differential equations (ODEs) of population data by the Stochastic Approximation Expectation-Maximization (SAEM) algorithm;  
5. `gnlmm` for generalized non-linear mixed-models (possibly defined by ordinary differential equations) of population data by the adaptive Gaussian quadrature algorithm.

A few utilities to facilitate population model building are also included in `nlmixr`.

For a brief vignette, please see:
https://github.com/nlmixrdevelopment/nlmixr/blob/master/inst/nlmixr-intro.pdf

The examples in the vignette can be run using VignetteDemo.R and the associated data files available at:
https://github.com/nlmixrdevelopment/nlmixr/tree/master/vignettes

For PKPD modeling (with ODE and dosing history) with [Stan](http://mc-stan.org/), check out Yuan's package PMXStan: https://github.com/yxiong1/pmxstan. 

# Installation on Windows
To replicate the environment that was used in Windows for `nlmixr` development, you will need administrator rights, and you should perform the following steps:

1. Install R 3.4.2 (or better) from the R website.
   - For best results, we suggest you use `C:\R\R-3.4.2`, but
     you can also use the default location (`C:\Program Files\R\R-3.4.2`) as well, if really needed.
   - For 64-bit Windows, it is best practice to include *only* the 64-bit version. 
     If you include 32-bit files, some packages may not
     run correctly.  Additionally, both the 32- and
     64-bit binaries have to be compiled for every package. Similarly, if on 32-bit Windows, install only the 32-bit version of R (and Python, and Rtools).

2. Install the appropriate version of Rtools for Windows, currently version 3.4, from [here](https://cran.r-project.org/bin/windows/Rtools/).
   - This is an absolute requirement, since it includes C++ and related compilers not usually available under Windows.
   - For best results, use the default location of `c:\Rtools`
     - `RxODE`, a required component of `nlmixr`, checks and sets up the path based on the following:
	    a. `Rtools` is in the path (fastest and recommended option)
		b. `Rtools` was installed with information saved to the Windows registry, and `RxODE` can 
		   find the installation.
		c. `Rtools` is on a hard drive installed in either `Rtools` or `RBuildTools`
     - If you are on 64-bit windows, please *do not install* the R
       3.3.x 32-bit toolchain.  These files can interfere with some
       packages that compile binaries, with unpredictable consequences.  Similarly, only install 32-bit
       Rtools on 32-bit versions of Windows.
   - Make sure the compilers have been added to the Windows `PATH` environment variable, or `RxODE` and `nlmixr` will not work (this should be done automatically during installation).
3. Install a version of Python for Windows.
   - This is used for its symbolic algebra package [SymPy](http://sympy.org/).
   - A very robust Python distribution that includes [SymPy](http://sympy.org/) and
     many packages that may be useful to the data scientist and/or
     pharmacometrician
     is [Anaconda](https://www.anaconda.com/download/). Although very straightforward and easy to install, it is quite a large download and contains much more than you will need to run `nlmixr`. When installing, use the Python 3.6 version. During the installation, Anaconda provides the option of adding itself to the `PATH` environment variable, but advises against it; please do this anyway (despite the red warning).
   - Another option is to use [official Python](http://python.org), although you will need to install [SymPy](http://sympy.org/) separately if you go this route, which is sometimes not straightforward under Windows 10 owing to folder permissions. Nonetheless, see [here](http://simpy.readthedocs.io/en/latest/simpy_intro/installation.html) for instructions for installation from source or using `pip`. Note that if you approach us for support, we are going to recommend that you use [Anaconda](https://www.anaconda.com/download/).
   - Regardless of the option you choose, please use like with like (64-bit Python for 64-bit Windows, for example).
   - Note that using the official Python may result in some issues with write permissions on Windows 10 - see [here](https://stackoverflow.com/questions/31172719/pip-install-access-denied-on-windows) for a few workarounds.
   - Once again, make sure Python has been added to the Windows `PATH` environment variable, or `RxODE` and `nlmixr` will not work, no matter what Anaconda might say.   
3. Install `devtools`.
   - This package is required to install packages from Github, amongst other things.
   - This can be done from a clean R session by `install.packages("devtools")`.
4. Load `devtools` using `library(devtools)`
5. Install `RxODE`.
   - Currently the new version of `RxODE` is in the process of being
     uploaded to CRAN.  `nlmixr` needs this newer version of `RxODE` to
     function correctly. To install this version, use the command:
     `install_github("nlmixrdevelopment/RxODE")`.
   - Once installed, type `RxODE::rxWinPythonSetup()` to install the required package `SnakeCharmR`and to make sure Python and SymPy are working properly.
   - Restart your R session.
   - As a quick test, you can make sure that R and Python can
     communicate by typing the command `library(SnakeCharmR)`.
   - To validate or test the installation of `RxODE` completely, you
     can type the following `library(RxODE); rxTest();` and it will
     run all of the unit tests in RxODE to make sure it is running
     correctly on your system.
6. Install `nlmixr`.
   - This can be done by `install_github("nlmixrdevelopment/nlmixr")`

# Installation on Linux
Installation on Linux is much more straightforward, since many prerequisites (such as compilers and Python) are already available. Details TBA.

# Installation on macOS
Installation on macOS is much more straightforward, since many prerequisites (such as compilers and Python) are already available. Details TBA.

