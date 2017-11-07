summary.apc <-
function(object, ...){
  cat("Error control: ", object$ErrorType, " at ", object$Level, "\n")
  print(object$Results, quote = FALSE, row.names = FALSE)
  cat("---", "\n")
  cat("Minimum APC:", object$apc, "\n")
  cat("Penalties:", object$Penalties, "\n")
  cat("Active Effects:", object$ActEffs, "\n")
  cat("Non-active Effects:", object$NonActEffs, "\n")
}
