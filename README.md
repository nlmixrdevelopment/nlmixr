![alt tag](https://github.com/nlmixrdevelopment/nlmixr/blob/master/logo.png)

# nlmixr: an R package for population PKPD modeling
***  

#####Authors: Matthew Fidler, Yuan Xiong, Rik Schoemaker, Justin Wilkins, Mirjam Trame, Wenping Wang

***
`nlmixr` is an R package for fitting general dynamic models, pharmacokinetic (PK) models and pharmacokinetic-pharmacodynamic (PKPD) models in particular, with either individual data or population data. `nlmixr` has five main modules:  1) `dynmodel()` and its mcmc cousin `dynmodel.mcmc()` for nonlinear dynamic models of individual data; 2) `nlme_lin_cmpt()`for one to three linear compartment models of population data with first order absorption, or i.v. bolus, or i.v. infusion using the nlme algorithm; 3) `nlme_ode()` for general dynamic models defined by ordinary differential equations (ODEs) of population data using the nlme algorithm; 4) `saem_fit` for general dynamic models defined by ordinary differential equations (ODEs) of population data by the Stochastic Approximation Expectation-Maximization (SAEM) algorithm;  5) `gnlmm` for generalized non-linear mixed-models (possibly defined by ordinary differential equations) of population data by the adaptive Gaussian quadrature algorithm.

A few utilities to facilitate population model building are also included in `nlmixr`.

For a brief vignette, please see:
https://github.com/nlmixrdevelopment/nlmixr/blob/master/inst/nlmixr-intro.pdf

The examples in the vignette can be run using VignetteDemo.R and the associated data files available at:
https://github.com/nlmixrdevelopment/nlmixr/tree/master/vignettes

For PKPD modeling (with ODE and dosing history) with Stan, check out Yuan's package PMXStan: https://github.com/yxiong1/pmxstan

# Installation in Windows
To replicate the environment that was used in windows for nlmixr development, you will need administrator rights, and you should perform the following steps:

1. Install R 3.4.2 from the R website
   - Install R to a user writable location; I use `c:\R\R-3.4.2`, but
     you can also use the default location `C:\Program
     Files\R\R-3.4.2`, as long as it is user writable.  If you have
     admin access you can make it user writable by the following procedure:
	 - In Windows Explorer, right click the directory and Select
       "Properties/Security", 
	 - then click the "Edit" button with the shield next to it, 
	 - then click "Users", click the check box under "Full control", click "Apply", and "OK"
	   twice.
   - For 64 bit windows, it is best practice to *uncheck* the 32 bit
     installation files.  If you include them, some packages may not
     run correctly.  Additionally, you have to compile both the 32 and
     64 bit binaries for every package.

2. Install Rtools for windows version 3.4
   - This allows for fast solving of ODEs and faster estimation
   - For best results, use the default location of `c:\Rtools`
     - `RxODE`, a required component of `nlmixr` checks and sets up the path based on the following:
	    a. `Rtools` is in the path (fastest and recommended option)
		b. `Rtools` was installed with information saved to the windows registry, and `RxODE` can 
		   find the installation.
		c. `Rtools` is on a hard drive installed in either `Rtools` or `RBuildTools`
     - If you are on 64 bit windows, please *do not install* the R
       3.3.x 32 bit toolchain.  These files can interfere with some
       packages that compile binaries.  Similarly, only install 32 bit
       on 32 bi windows
3. Install python for windows;
   - This is used for its symbolic algebra package [SymPy](http://sympy.org/).
   - A very robust python distribution that includes [SymPy](http://sympy.org/) and
     many packages that may be useful to the data scientist and/or
     pharmacometrician
     is [anacodna](https://www.anaconda.com/download/).
   - Another option is to use python from
     the [official pyhton](http:://python.org) website.
   - Regardless of the option you choose, please use 64 bit python for 64 bit windows.
   - Also, like R, make sure that the users have full control of this
     directory.  If you have admin access you can adjust this as
     follows:
      - In the install location, using  Windows Explorer, right click the directory
      - Select "Properties/Security", 
      - then click the "Edit" button with the shield next to it, 
      - then click "Users", 
      - click the check box under "Full control", 
	  - click "Apply", and "OK" twice.
3. Install devtools
   - This package is required to install packages off of the github website.
   - This can be done from a clean R session by `install.packages("devools")`
4. Load devtools by `library(devtools)`
6. Install RxODE
   - Currently the new version of RxODE is in the process of being
     sent to CRAN.  nlmixr needs this newer version of RxODE to
     function correctly. To install this version, use the command:
     `install_github("nlmixrdevelopment/RxODE")`
   - Once installed, type `RxODE::rxWinPythonSetup()`
   - Restart your R session
   - As a quick test, you can make sure that R and python can
     communicate by typing the command `library(SnakeCharmR)`.
   - To validate or test the installation of `RxODE` completely, you
     can type the following `library(RxODE); rxTest();` and it will
     run all of the unit tests in RxODE to make sure it is running
     correctly on your system.
7. Install nlmixr
   - This can be done by `install_github("nlmixrdevelopment/nlmixr", ref="dparser-saem")`
