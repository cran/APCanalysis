apc <-
function(y, x, maxsize, level = 0.05, method = 1, data = NULL, effnames = NULL, reps = 50000, dp = 4){

  if(is.element("formula", class(y))){
    formula <- y
    rm(y)
    options(contrasts=c('contr.sum','contr.poly'))
    mod <- model.frame(formula, data = data)
    x <- model.matrix(formula, data = data)[ , -1]
    y <- model.extract(mod, "response")
  }

  m <- as.integer(maxsize)
  method <- as.integer(method)

  nr <- nrow(x)
  k <- ncol(x)
 
  xtx<-t(x)%*%x
  if(!isTRUE(all.equal(xtx[lower.tri(xtx)],(rep(0,((k*(k-1))/2)))))) stop("model matrix is not orthogonal")

  effnames <- colnames(x)
  if (is.null(effnames)) effnames <- paste("C", 1:nr, sep = "")

  
  if (method == 1) cs <- IERpenalties(n=nr, k = k, m = m, ier = level, reps = reps)
  if (method == 2) cs <- EERpenalties(n=nr, k = k, m = m, eer = level, reps = reps)
  if (method == 3) cs <- FDRpenalties(n=nr, k = k, m = m, fdr = level, reps = reps)
  
# cts <- cs[(m+1):1]
  bhat <- lm(y ~ x)$coefficients
  rss0 <- sum((lm(y ~ 1)$residuals) ^ 2)
  rssE <- (bhat[-1] ^ 2) * diag(xtx)
  ord.rssE <- order(rssE, decreasing = TRUE)
  
  Ests <- bhat[c(1, 1 + ord.rssE)]
  Col <- as.integer(c(0, ord.rssE))
  APC <- c(log(c(rss0,rss0-cumsum(rssE[ord.rssE])[1:m])) + cs,rep(NA,k-m))
  K <-  as.integer(which.min(APC))
  Active <- rep(" ", m + 1)
  Active[1:K] <- "*  "
  apc <- APC[K]

  Results <- data.frame(Model = c("intercept",paste("+", effnames[ord.rssE[1:m]])), 
                        Size = 0:m, RSS = round(c(rss0, rss0 - cumsum(rssE[ord.rssE])[1:m]), dp),
                        APC = round(APC[1:(m + 1)], dp), Active = Active)

  effnames.ord <- effnames[ord.rssE]
  if(isTRUE(identical(K,as.integer(1)))) active <- "none" else  active <- paste(sort(effnames.ord[1:(K-1)]), collapse=", ")
  if(isTRUE(identical(K,as.integer(k+1)))) nonactive <- "none" else  nonactive <- paste(sort(effnames.ord[K:length(effnames)]), collapse=", ")
  
out <- list(Results = Results, Penalties = cs, Level = level, ErrorType = c("IER","EER","FDR")[method], k = k, maxsize = maxsize, apc = apc, Ests = Ests, ActEffs = active, NonActEffs = nonactive)

class(out) = "apc"
return(out)

}
