# Installation test function
test_install <- function(){
  # Test 1: Correct R version
  if(sessionInfo()$R.version$major=="3" & as.numeric(sessionInfo()$R.version$minor)>=4.1){
    cat(paste0("Correct R version: Yes, ",sessionInfo()$R.version$version.string,"\n"))
  }else{
    cat("Correct R version: No, go to https://cloud.r-project.org/ and install latest version\n")
    return("Installation not complete!")
  }
  tmp <- try(library(RxODE), silent=TRUE);
  if (inherits(tmp, "try-error"))  {
      cat("RxODE installed: No; check https://github.com/nlmixrdevelopment/nlmixr for instructions to do so\n")
      return("Installation not complete!")
  } else {
      cat(paste0("RxODE installed: Yes\n"))
  }
  # Test 2: Python installation
  if(Sys.info()['sysname']%in%c('Darwin','Linux')){
    try(system("python -V > pv.txt 2>&1",intern=TRUE),silent=TRUE)
  }else if(Sys.info()['sysname']=='Windows'){
    try(shell("python -V > pv.txt 2>&1",intern=TRUE),silent=TRUE)
  }
  res1 <- readLines("pv.txt")
    if(any(grepl("not found|recognized",res1))){
    cat("Python installed: No, check https://github.com/nlmixrdevelopment/nlmixr for instructions to do so\n")
    return("Installation not complete!")
  }else{
    cat(paste0("Python installed: Yes, ",res1,"\n"))
  }

  # Test 3: sympy installation
  if(Sys.info()['sysname']%in%c('Darwin','Linux')){
    try(system("python -c \"import sympy\" > pv.txt 2>&1",intern=TRUE),silent=TRUE)
  }else if(Sys.info()['sysname']=='Windows'){
    try(shell("python -c \"import sympy\" > pv.txt 2>&1",intern=TRUE),silent=TRUE)
  }
  res2 <- readLines("pv.txt")
  if(length(res2)==0){
    cat("sympy installed: Yes\n")
  }else{
    cat("sympy installed: No, check https://github.com/nlmixrdevelopment/nlmixr for instructions to do so\n")
    return("Installation not complete!")
  }

  # Test 4: devtools installed
  if("devtools" %in% rownames(installed.packages())){
    cat("devtools package installed: Yes\n")
  }else{
    cat("devtools package installed: No, use install.packages('devtools') to do so\n")
    return("Installation not complete!")
  }

  # Test 5: Rtools installed
  if(devtools::find_rtools()){
    cat("Rtools installed: Yes\n")
  }else{
    cat("Rtools installed: No, check https://github.com/nlmixrdevelopment/nlmixr for instructions to do so\n")
    return("Installation not complete!")
  }

  # Test 6: RxODE installed
  if("RxODE" %in% rownames(installed.packages())){
    cat("RxODE package installed: Yes\n")
  }else{
    cat("RxODE package installed: No, use install_github('nlmixrdevelopment/RxODE') to do so\n")
    return("Installation not complete!")
  }

  # Test 7: SnakeCharmR installed
  if("SnakeCharmR" %in% rownames(installed.packages())){
    cat("SnakeCharmR package installed: Yes\n")
  }else{
    cat("SnakeCharmR package installed: No, use install_github('nlmixrdevelopment/SnakeCharmR') to do so\n")
    return("Installation not complete!")
  }

  # Test 8: nlmixr installed
  if("nlmixr" %in% rownames(installed.packages())){
    cat("nlmixr package installed: Yes\n")
  }else{
    cat("nlmixr package installed: No, use install_github('nlmixrdevelopment/nlmixr') to do so\n")
    return("Installation not complete!")
  }

  # Test 9: xpose.nlmixr installed
  if("xpose.nlmixr" %in% rownames(installed.packages())){
    cat("xpose.nlmixr package installed: Yes\n")
  }else{
    cat("xpose.nlmixr package installed: No, use install_github('nlmixrdevelopment/xpose.nlmixr') to do so\n")
    return("Installation not complete!")
  }

  # Test 10: shinyMixR installed
  if("shinyMixR" %in% rownames(installed.packages())){
    cat("shinyMixR package installed: Yes\n")
  }else{
    cat("shinyMixR package installed: No, use install_github('richardhooijmaijers/shinyMixR') to do so\n")
    return("Installation not complete!")
  }

  # Test 11: nlme test for theophylline
  library(nlmixr)
  testmod <- function() {
    ini({
      tka <- .5
      tcl <- -3.2
      tv <- -1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      add.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.err)
    })
  }
  if(Sys.info()['sysname']%in%c('Darwin','Linux')){
    {sink("/dev/null"); fit1 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="nlme"))); sink(); }
    {sink("/dev/null"); fit2 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="saem"))); sink(); }
  }else{
    {sink("NUL"); fit1 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="nlme"))); sink(); }
    {sink("NUL"); fit2 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="saem"))); sink(); }
  }
  if(all(class(fit1)%in%c("nlmixr.ui.nlme","focei.fit","data.frame"))){
    cat("nlmixr run under nlme: Yes\n")
  }else{
    cat("nlmixr run under nlme: No, contact nlmixr team\n")
    return("Installation not complete!")
  }
  if(all(class(fit2)%in%c("nlmixr.ui.saem","focei.fit","data.frame"))){
    cat("nlmixr run under saem: Yes\n")
  }else{
    cat("nlmixr run under saem: No, contact nlmixr team\n")
    return("Installation not complete!")
  }
  cat("---- Installation test finished! ----\n")
}

try({test_install()})

cat("Begin Session Info:\n")
print(sessionInfo())
cat("\n")
