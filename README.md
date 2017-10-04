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

1. Install R 3.4.1 from the R website
   - Install R to a user writable location; I use `c:\R\R-3.4.1`, but you can also use the default location `C:\Program Files\R\R-3.4.1`. In either case check to make sure you have rights to write in this directory: In Windows Explorer, right click the directory and Select “Properties/Security”, then click the “Edit” button with the shield next to it, then click “Users”, click the check box under “Full control”, click “Apply”, and “OK” twice.
   - For 64 bit windows, make sure to *uncheck* the 32 bit installation files.  They have been known to interfere with R in the past.
2. Install Rtools for windows version 3.4
   - This allows for fast solving of ODEs and faster estimation
   - For best results, use the default location of `c:\Rtools`
   - Please *do not install* the R 3.3.x 32 bit toolchain.  These files can interfere with nlmixr
3. Install python for windows 
   - This is used for its symbolic algebra package [SymPy](http://sympy.org/).
   - In 64-bit windows, the best Python to install is Python 3.6.2 using
     `python-3.6.2-amd64.exe`
     https://www.python.org/downloads/release/python-362/ scroll down to select the ‘Windows x86-64 executable installer’ and run the installer. 
   - Ensure ‘Install launcher for all users’ and ‘Add Python 3.6 to PATH’ are checked. 
   - Then choose ‘Customize installation’: make sure all options are checkmarked on the next page, click ‘Next’ and then checkmark ‘Install for all users’. This changes the Python install location to `C:\Program Files\Python36`. Then click ‘Install’.
   - Next, add the PYTHONPATH environment variable, as follows: In the Windows search box, type “environment”, and choose “Edit the system environment variables” when it appears as a suggested option. In the System Properties box that appears, click on “Environment Variables…” In the “System variables” box (lower half of the window), click “New” and create PYTHONPATH, with value `C:\Program Files\Python36; C:\Program Files\Python36\Lib` (assuming the default installation path has been used), and click “OK”. 
   - Please check that the environmental variable `PYTHONHOME` is not set and `PYTHONPATH` is set to the correct location.
   - The default permissions of the Python directory under Program Files must be adjusted to ensure smooth installation of packages. In Windows Explorer, right click the directory `C:\Program Files\Python36` and Select “Properties/Security”, then click the “Edit” button with the shield next to it, then click “Users”, click the check box under “Full control”, click “Apply”, and “OK” twice.  
3. Install devtools
   - This package is required to install packages off of the github website.
   - This can be done from a clean R session by `install.packages("devools")`
4. Load devtools by `library(devtools)`
5. Install dparser
   - This package is available on CRAN. 
     - For vanilla R, you can use the command: `install.packages("dparser")`
     - For MRAN, you should use devtools to install the latest version
       of dparser since the package was not available when R 3.4.1 was
       released: `install_github("nlmixrdevelopment/dparser-R")`
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
