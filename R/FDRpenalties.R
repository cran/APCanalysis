FDRpenalties <-
function(n, k = n - 1, m = min(n - 2, k), fdr = .1, reps = 50000, rnd = 3){
# Checks on conditions that must be satisfied
if(!isTRUE(all.equal(n%%4,0))) stop("n must be a multiple of 4")
if(!isTRUE(all.equal(k%%1,0))) stop("k must be an integer")
if(!isTRUE(all.equal(m%%1,0))) stop("m must be an integer")
if(k > n - 1) stop("k cannot be greater than n-1")
if(k < 1) stop("k cannot be less than 1")
if(m > k) stop("m cannot be greater than k")
if(m < 1) stop("m cannot be less than 1")
if(m > (n-2)) stop("m cannot be greater than n-2")
if(m > k) stop("m cannot be greater than k")
if(fdr <= 0) stop("fdr must be greater than 0")
if(fdr >= 1) stop("fdr must be less than 1")

if (fdr > 1/m) warning("fdr > 1/m which results in some penalties being equal")
# Algorythm callculates differences between penalties 
# starting with diff(m-1) = pen(m) - pen(m-1) 

cs <- NULL
startj <- m - 1

# If fdr >= 1 then diff(j)=0 for sufficiently large j 

if((fdr * m) >= 1){
   num <- floor(fdr * m)
   cs <- rep(0, num)
   startj <- m - 1 - num
}


# Loop that estimates diff(j) for j = starj, startj-1, ... 1, 0.
# Estimate of diff(j) is based on assuming j large active effects.

for(j in startj:0){

  # Create matrix of squared random N(0,1) observations. 
  # Number of columns is n-1-j which equals inactive columns (k-j) plus unused columns (n-1-k).

  sqres <- matrix(rnorm(reps * (n - 1 - j)) ^ 2, reps, n - 1 - j)

  # If there is more than one inactive column (k!=j+1) then sort entries for inactive columns.
  # Inactive  columns are the last (k-j) columns.

  if((n - k) != (n - 1 - j)) sqres[, (n - k):(n - 1 - j)] <- t(apply(sqres[ ,(n - k):(n - 1 - j)], 1, sort))

  # Find RSS for models containing just the j active effects, the j-effect model + 1, ... the j-effect model + m-j.

  lRSS <- log(apply(sqres, 1, cumsum)[(n - m - 1):(n - 1 - j), ])
  d1<-dim(lRSS)[1]

  # If d1==2 then m = j+1. In this case at most one variable is being added. 
  # The differences in log(RSS) are found and the relevant quantile taken to estimate diff(j). 

  if(d1 == 2){
    out <- lRSS[2, ]- lRSS[1, ]
    cs <- as.numeric(quantile(out, 1 - (fdr * m))) 
  }

  # If m> j+1 then the maximum number of additional variables is >=2. 
  # The models that add >=1 variable are compared and the one that will minimize 
  # APC* identified. For this model the difference in log(RSS) between this 
  # this model and the j-variable model plus its current penalty is recorded in out
  # and the number of additional variables in wts. The number of allowable errors are 
  # calculated (toterrs) and the value of diff(j) that allows this to be achieved 
  # is estimated (newc) and the current list of penalties is updated.

  if(d1 > 2){
    out <- lRSS[d1, ] - apply((lRSS[-d1, ] + c(cs, 0)), 2, min)
    wts <- d1 - apply((lRSS[-d1, ] + c(cs, 0)), 2, order)[1, ]
    wts <- wts / (wts + j)
    ord <- order(out, decreasing = TRUE)
    oout <- out[ord]
    owts <- cumsum(wts[ord])
    toterrs <- reps * fdr
    newc <- min(oout[owts <= toterrs])
    cs <- c(cs + newc, newc) 
  }
}


# The set of estimated penalties is returned.

cs <- round(c(0, cs[length(cs):1]), rnd)
attributes(cs) <- NULL
return(cs)
}
