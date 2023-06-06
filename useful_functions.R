#### adstock_geometric

adstock_geometric <- function(x, theta) {
  stopifnot(length(theta) == 1)  # check whether len(theta) == 1 
  if (length(x) > 1) {
    x_decayed <- c(x[1], rep(0, length(x) - 1))
    for (xi in 2:length(x_decayed)) {
      x_decayed[xi] <- x[xi] + theta * x_decayed[xi - 1]
    }
    thetaVecCum <- theta
    for (t in 2:length(x)) {
      thetaVecCum[t] <- thetaVecCum[t - 1] * theta
    } # plot(thetaVecCum)
  } else {
    x_decayed <- x
    thetaVecCum <- theta
  }
  inflation_total <- sum(x_decayed) / sum(x)
  return(list(x = x, x_decayed = x_decayed, thetaVecCum = thetaVecCum, inflation_total = inflation_total))
}


##### adstock weibull

adstock_weibull <- function(x, shape, scale, windlen = length(x), type = "cdf") {
  stopifnot(length(shape) == 1)
  stopifnot(length(scale) == 1)
  if (length(x) > 1) {
    check_opts(tolower(type), c("cdf", "pdf"))
    x_bin <- 1:windlen
    scaleTrans <- round(quantile(1:windlen, scale), 0)
    if (shape == 0 | scale == 0) {
      x_decayed <- x
      thetaVecCum <- thetaVec <- rep(0, windlen)
      x_imme <- NULL
    } else {
      if ("cdf" %in% tolower(type)) {
        thetaVec <- c(1, 1 - pweibull(head(x_bin, -1), shape = shape, scale = scaleTrans)) # plot(thetaVec)
        thetaVecCum <- cumprod(thetaVec) # plot(thetaVecCum)
      } else if ("pdf" %in% tolower(type)) {
        thetaVecCum <- .normalize(dweibull(x_bin, shape = shape, scale = scaleTrans)) # plot(thetaVecCum)
      }
      x_decayed <- mapply(function(x_val, x_pos) {
        x.vec <- c(rep(0, x_pos - 1), rep(x_val, windlen - x_pos + 1))
        thetaVecCumLag <- lag(thetaVecCum, x_pos - 1, default = 0)
        x.prod <- x.vec * thetaVecCumLag
        return(x.prod)
      }, x_val = x, x_pos = seq_along(x))
      x_imme <- diag(x_decayed)
      x_decayed <- rowSums(x_decayed)[seq_along(x)]
    }
  } else {
    x_decayed <- x_imme <- x
    thetaVecCum <- 1
  }
  inflation_total <- sum(x_decayed) / sum(x)
  return(list(
    x = x,
    x_decayed = x_decayed,
    thetaVecCum = thetaVecCum,
    inflation_total = inflation_total,
    x_imme = x_imme)
  )
}

.normalize <- function(x) {
  if (diff(range(x)) == 0) {
    return(c(1, rep(0, length(x) - 1)))
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}

####  creating a vector with sample function
# output include original, decayed, decay rate, inflation total
test = runif(n=100, min=1, max=20)

re = adstock_geometric(test, 0.5)

shape = 0.5
scale = 0.5
re_cdf = adstock_weibull(test, shape = shape, scale = scale, type = "cdf")
re_pdf = adstock_weibull(test, shape = shape, scale = scale, type = "pdf")



  