## This is only for nlmixr
for (f in c("src/RcppExports.cpp")) {
  l <- readLines(f)
  w <- which(regexpr("^[#]include <RcppArmadillo.h>", l) != -1)
  if (length(w) == 1) {
    l <- l[-w]
    message("Excluding RcppArmadillo from", f)
    writeLines(l, f)
  }
}

if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  writeLines(gsub("@ISYSTEM@", "I",
                  gsub("@CXX14STD@", "CXX14STD = -std=c++1y",
                       suppressWarnings(readLines("src/Makevars.in")))),
             "src/Makevars.win")
} else {
  writeLines(gsub("@ISYSTEM@", "isystem",
                  gsub("@CXX14STD@", "CXX14STD = -std=gnu++14",
                       suppressWarnings(readLines("src/Makevars.in")))),
             "src/Makevars")
}
