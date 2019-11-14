message("Optimize:", Sys.getenv("binOpt"))
message("CRAN:", Sys.getenv("useCRAN"))
message("RxODE:", Sys.getenv("rxodeRef"))
message("nlmixr:", Sys.getenv("nlmixrRef"))

.rtoolsBase <- file.path(getwd(), "R", "Rtools")
.extraBin <- normalizePath(file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")))
Sys.setenv(BINPREF=gsub("\\\\", "/", paste0(.extraBin, "/")))
.extraBin <- paste0(normalizePath(file.path(.rtoolsBase, "bin")), ";", .extraBin)
Sys.setenv(PATH=paste0(.extraBin, ";", Sys.getenv("PATH")))

#invisible(readline(prompt="Press [enter] to continue"))
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
        "CXX11FLAGS=-O2\n",
        file = M, sep = "\n")
} else {
    cat("\n",
        "\nCXX14FLAGS=-O3 -march=native",
        "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
        "CXX11FLAGS=-O3 -march=native\n",
        file = M, sep = "\n")
}
message(M)
message("================================================================================")
message(paste0(readLines(M), collapse="\n"))

if (Sys.getenv("useCRAN") == "true"){
    if (Sys.getenv("rxodeRef") == ""){
        install.packages("RxODE", type="source", repos="https://cloud.r-project.org")
    } else {
        devtools::install_version("RxODE", version=Sys.getenv("rxodeRef"), repos="https://cloud.r-project.org")
    }
    if (Sys.getenv("nlmixrRef") == ""){
        install.packages("nlmixr", type="source", repos="https://cloud.r-project.org")
    } else {
        devtools::install_version("nlmixr",version=Sys.getenv("nlmixrRef"), repos="https://cloud.r-project.org")
   }
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


