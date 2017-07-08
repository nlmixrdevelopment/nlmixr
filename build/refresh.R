## message("Update init.c for R based registration of routines.")
## unlink(devtools::package_file("src/init.c"))
## tools::package_native_routine_registration_skeleton(devtools::package_file(""),
##                                                     devtools::package_file("src/init.c"))

## lines <- readLines(devtools::package_file("src/init.c"))
## ## These are the c functions not created in nlmixr
## ## FIXME: why?
## reg <- rex::rex(or("RxODE_mod_m1_ode_solver"))

## lines <- lines[regexpr(reg, lines) == -1]
## writeLines(lines, devtools::package_file("src/init.c"));
