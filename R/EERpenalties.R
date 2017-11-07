EERpenalties <-
function(n, k = n - 1, m = min(n - 2, k), eer = .20, reps = 50000, rnd = 3){
  if(!isTRUE(all.equal(n %% 4, 0))) stop("n must be a multiple of 4")
  if(!isTRUE(all.equal(k %% 1, 0))) stop("k must be an integer")
  if(!isTRUE(all.equal(m %% 1, 0))) stop("m must be an integer")
  if(k > n - 1) stop("k cannot be greater than n-1")
  if(k < 1) stop("k cannot be less than 1")
  if(m > k) stop("m cannot be greater than k")
  if(m < 1) stop("m cannot be less than 1")
  if(m > (n - 2)) stop("m cannot be greater than n-2")
  if(eer <= 0) stop("EER must be greater than 0")
  if(eer >= 1) stop("EER must be less than 1")

# Algorythm callculates differences between penalties 
# starting with diff(m-1) = pen(m) - pen(m-1) 

  cs <- NULL
  startj <- m - 1

# The value of diff(m-1) can (under a certain condition) be calculated analytically.

  if(qf((1 - eer / (k + 1 - m)), 1, (n - 1 - m)) > n - 1 - m){
    cs <- log(qf((1 - eer / (k + 1 - m)), 1, (n - 1 - m)) / (n - 1 - m) + 1)
    startj <- m - 2
  }

# Stop if m=1 and diff(0) was calculated above.

  if(startj < 0){
    cs <- round(c(0, cs), rnd)
    return(cs)
  }

# Loop that estimates diff(j) for j = starj, startj-1, ... 1, 0.
# Estimate of diff(j) is based on assuming j large active effects.

  for(j in startj:0){

  # Create matrix of squared random N(0,1) observations. 
  # Number of columns is n-1-j which equals inactive columns (k-j) plus unused columns (n-1-k).

    sqres<-matrix(rnorm(reps*(n-1-j))^2,reps,n-1-j)

  # If there is more than one inactive column (k!=j+1) then sort entries for inactive columns.
  # Inactive  columns are the last (k-j) columns.

    if((n - k) != (n - 1 - j)) sqres[ , (n - k):(n - 1 - j)] <- t(apply(sqres[ , (n - k):(n - 1 - j)], 1, sort))

  # Find RSS for models containing just the j active effects, the j-effect model + 1, ... the j-effect model + m-j.

    lRSS <- log(apply(sqres, 1, cumsum)[(n - m - 1):(n - 1 - j), ])
    d1 <- dim(lRSS)[1]

  # If d1==2 then m = j+1. In this case at most one variable is being added. 
  # The differences in log(RSS) are found and the relevant quantile taken to estimate diff(j). 

    if(d1 == 2){
      out <- lRSS[2, ] - lRSS[1, ]
      cs <- as.numeric(quantile(out, 1 - eer))
    }

  # If m> j+1 then the maximum number of additional variables is >=2. 
  # The models that add >=1 variable are compared and the one that will minimize 
  # APC* identified. For this model the difference in log(RSS) between this 
  # this model and the j-variable model plus its current penalty is recorded in out.
  # The value of diff(j) that allows the specified EER to be achieved 
  # is estimated (newc) and the current list of penalties is updated.

    if(d1 > 2){
      out <- lRSS[d1, ] - apply((lRSS[-d1, ] + c(cs, 0)), 2, min)
      newc <- quantile(out, 1 - eer)
      cs <- c(cs + newc, newc)
    }
  }

cs <- round(c(0, cs[length(cs):1]), rnd)
attributes(cs) <- NULL
return(cs)
}
