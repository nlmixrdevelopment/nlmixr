
preconditionFit <- function(fit){
  pre <- preCondInv(fit$R)
  P <- symengine::Matrix(pre)
  d0 <- dimnames(fit$R)[[1]]
  d <- paste0("nlmixrPre_", dimnames(fit$R)[[1]])
  d2 <- symengine::Matrix(d)
  modExtra <- paste(paste0(d0, "=", sapply(as.vector(P%*%d2),as.character)), collapse="\n")
  preInv <- solve(pre)
  return(list(modExtra, setNames(as.vector(preInv %*% matrix(fit$theta[d0])), d)))
}
