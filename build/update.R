.libPaths(.libPaths()[regexpr("nlmixr",.libPaths())!=-1])

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win") #
if (!file.exists(M)) file.create(M)

## In CRAN and installer should use -O2 and drop the -march
if (Sys.getenv("binOpt") == "false"){
    cat("\n",
        "\nCXX14FLAGS=-O2",
        "CXX14 = $(BINPREF)g++ -std=c++1y",
        "CXX11FLAGS=-O2",
        file = M, sep = "\n", append = TRUE)
} else {
    cat("\n",
        "\nCXX14FLAGS=-O3 -march=native",
        "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
        "CXX11FLAGS=-O3 -march=native",
        file = M, sep = "\n", append = TRUE)
}
if (Sys.getenv("useCRAN") == "true"){
    install.packages("RxODE", type="source")
    install.packages("nlmixr", type="source")
} else {
    if (Sys.getenv("rxodeRef") == ""){
        devtools::install_github("nlmixrdevelopment/RxODE", dependencies=FALSE)
    } else {
        devtools::install_github("nlmixrdevelopment/RxODE", ref=Sys.getenv("rxodeRef"), dependencies=FALSE)
    }
    if (Sys.getenv("nlmixrRef") == ""){
        devtools::install_github("nlmixrdevelopment/nlmixr", dependencies=FALSE)
    } else {
        devtools::install_github("nlmixrdevelopment/nlmixr", ref=Sys.getenv("nlmixrRef"), dependencies=FALSE)
    }
}
## Install CRAN from source
## Install


