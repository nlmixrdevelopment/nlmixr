##' nlmixr fits population PK and PKPD non-linear mixed effects models.
##'
##' nlmixr is an R package for fitting population pharmacokinetic (PK)
##' and pharmacokinetic-pharmacodynamic (PKPD) models.
##' importFrom(Rcpp, evalCpp)
##' @importFrom brew brew
##' @importFrom lattice xyplot
##' @importFrom lattice trellis.par.get
##' @importFrom nlme nlme
##' @importFrom nlme groupedData
##' @importFrom nlme getData
##' @importFrom nlme pdDiag
##' @importFrom RxODE RxODE
##' @importFrom graphics abline lines matplot plot points title
##' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start
##' @importFrom utils assignInMyNamespace getFromNamespace head stack
##' @importFrom parallel mclapply
##' @importFrom lbfgs lbfgs
##' @importFrom methods is
##' @importFrom Rcpp evalCpp
##' @useDynLib nlmixr, .registration=TRUE
"_PACKAGE"

rex::register_shortcuts("nlmixr");
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values"));




nlmixr.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|"




nlmixr.logo.full <- "                                  %%%%%%(                                                    #%%%%%#,\n                                  %%%%%%(                                                   #%%%%%%%(,\n                                  %%%%%%(                                                   (%%%%%%%/,\n                                  %%%%%%(                                                    /%%%%%/,\n                                  %%%%%%( \n    (%%%%%%   /%%%%%%%%,          %%%%%%(     ,%%%%%%,  ,#%%%%%%%(      .*#%%%%%%%(.         %%%%%%%   (%%%%%%%,       . #%%%%%%#   #%%%%%  *%%%%%/        \n    (%%%%%%*%%%%%%%%%%%%%#        %%%%%%(     ,%%%%%%,(%%%%%%%%%%%%%,  #%%%%%%%%%%%%%*       %%%%%%%     %%%%%%%(      ,%%%%%%%(    #%%%%%,%%%%%%#         \n    (%%%%%%%%%%%%%%%%%%%%%%       %%%%%%(     ,%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/      %%%%%%%,     *%%%%%%%.   (%%%%%%%      #%%%%%%%%%%%%          \n    (%%%%%%%%/     .%%%%%%%(      %%%%%%(     ,%%%%%%%%#      %%%%%%%%%%#  .  .%%%%%%%%.     %%%%%%%       .%%%%%%%#.%%%%%%%,       #%%%%%%%*              \n    (%%%%%%%.       .%%%%%%%      %%%%%%(     ,%%%%%%%/        %%%%%%%%/     .  #%%%%%%,     %%%%%%%         (%%%%%%%%%%%%%.        #%%%%%%/               \n    (%%%%%%,         #%%%%%%      %%%%%%(     ,%%%%%%(         *%%%%%%#        .*%%%%%%*     %%%%%%%          ,%%%%%%%%%%/          #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%*         .%%%%%%*         .%%%%%%*     %%%%%%%            #%%%%%%%*           #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%           #%%%%%%%%%,          #%%%%%#                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%          %%%%%%%%%%%%#         #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%        (%%%%%%%*%%%%%%%.       #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%       %%%%%%%,   %%%%%%%#      #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%     /%%%%%%%.     /%%%%%%%,    #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%   .%%%%%%%/         %%%%%%%(   #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%  ,%%%%%%%.           (%%%%%%%  #%%%%%%                \n";

##' Messages the nlmixr logo...
##'
##' @param str String to print
##' @author Matthew L. Fidler
nlmixrLogo <- function(str="", version=sessionInfo()$otherPkgs$nlmixr$Version){
    message(sprintf(nlmixr.logo, str, version));
}
