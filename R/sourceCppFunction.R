## same named function copied from Rcpp
## Rcheck hack to avoid complains
sourceCppFunction <- function(func, isVoid, dll, symbol) {
  args <- names(formals(func))
  body <- quote(CALL_PLACEHOLDER(EXTERNALNAME, ARG))[c(
    1:2,
    rep(3, length(args))
  )]
  for (i in seq(along.with = args)) body[[i + 2]] <- as.symbol(args[i])
  body[[1L]] <- .Call
  body[[2L]] <- getNativeSymbolInfo(symbol, dll)$address
  if (isVoid) {
    body <- call("invisible", body)
  }
  body(func) <- body
  func
}
