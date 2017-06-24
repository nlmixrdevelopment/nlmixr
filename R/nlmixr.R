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
##' @useDynLib nlmixr
"_PACKAGE"

rex::register_shortcuts("nlmixr");


nlmixr.logo <- "                                  %%%%%%(                                                    #%%%%%#,\n                                  %%%%%%(                                                   #%%%%%%%(,\n                                  %%%%%%(                                                   (%%%%%%%/,\n                                  %%%%%%(                                                    /%%%%%/,\n                                  %%%%%%( \n    (%%%%%%   /%%%%%%%%,          %%%%%%(     ,%%%%%%,  ,#%%%%%%%(      .*#%%%%%%%(.         %%%%%%%   (%%%%%%%,       . #%%%%%%#   #%%%%%  *%%%%%/        \n    (%%%%%%*%%%%%%%%%%%%%#        %%%%%%(     ,%%%%%%,(%%%%%%%%%%%%%,  #%%%%%%%%%%%%%*       %%%%%%%     %%%%%%%(      ,%%%%%%%(    #%%%%%,%%%%%%#         \n    (%%%%%%%%%%%%%%%%%%%%%%       %%%%%%(     ,%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/      %%%%%%%,     *%%%%%%%.   (%%%%%%%      #%%%%%%%%%%%%          \n    (%%%%%%%%/     .%%%%%%%(      %%%%%%(     ,%%%%%%%%#      %%%%%%%%%%#  .  .%%%%%%%%.     %%%%%%%       .%%%%%%%#.%%%%%%%,       #%%%%%%%*              \n    (%%%%%%%.       .%%%%%%%      %%%%%%(     ,%%%%%%%/        %%%%%%%%/     .  #%%%%%%,     %%%%%%%         (%%%%%%%%%%%%%.        #%%%%%%/               \n    (%%%%%%,         #%%%%%%      %%%%%%(     ,%%%%%%(         *%%%%%%#        .*%%%%%%*     %%%%%%%          ,%%%%%%%%%%/          #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%*         .%%%%%%*         .%%%%%%*     %%%%%%%            #%%%%%%%*           #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%           #%%%%%%%%%,          #%%%%%#                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%          %%%%%%%%%%%%#         #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%        (%%%%%%%*%%%%%%%.       #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%       %%%%%%%,   %%%%%%%#      #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%     /%%%%%%%.     /%%%%%%%,    #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%   .%%%%%%%/         %%%%%%%(   #%%%%%%                \n    (%%%%%%          (%%%%%%      %%%%%%(     ,%%%%%%,         .%%%%%%*         .%%%%%%,     %%%%%%%  ,%%%%%%%.           (%%%%%%%  #%%%%%%                \n";

##' Messages the nlmixr logo...
##'
##' @author Matthew L. Fidler
nlmixrLogo <- function(){
    message(nlmixr.logo);
}
