plot.apc <-
function(x, elabs = TRUE,  ...){
 plot(c(-.5, x$m+.5), c(min(x$Results[,4],na.rm=TRUE), max(x$Results[,4],na.rm=TRUE)), type = "n", xlab = "Model size", ylab = "APC*")
 if (isTRUE(elabs)) text(0:x$m, x$Results[1:(x$m+1), 4], labels = names(x$Ests)[1:(x$m+1)]) else points(0:x$m, x$Results[1:(x$m+1), 4], pch = 19)
}
