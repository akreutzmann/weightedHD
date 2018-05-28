whdquantile <- function (x, probs = seq(0, 1, 0.25), se = FALSE, na.rm = FALSE, 
                                names = TRUE, weights = FALSE, smp_weight = NULL) 
{
  if (na.rm) {
    na <- is.na(x)
    if (any(na)) 
      x <- x[!na]
  }
  order <- order(x)
  x <- x[order]
  
  if (is.null(smp_weight)) { # Unweighted HD
    n <- length(x)
    m <- n + 1
    ps <- probs[probs > 0 & probs < 1]
    qs <- 1 - ps
    a <- outer((0:n)/n, ps, function(x, p, m) pbeta(x, p * m, 
                                                    (1 - p) * m), m = m)
  } else { # Weighted HD
    smp_weight <- smp_weight[order]
    n <- sum(smp_weight)
    if (n < 2) 
      return(rep(NA, length(probs)))
    m <- n + 1
    ps <- probs[probs > 0 & probs < 1]
    qs <- 1 - ps
    # Usage of cumulative sum of weights divided by sum of weights
    a <- outer(c(0,cumsum(smp_weight))/sum(smp_weight), ps, 
               function(x, p, m) pbeta(x, p * m, (1 - p) * m), m = m) 
  }
  w <- a[-1, ] - a[-(length(x) + 1), ]
  r <- drop(x %*% w)
  rp <- range(probs)
  pp <- ps
  if (rp[1] == 0) {
    r <- c(x[1], r)
    pp <- c(0, pp)
  }
  if (rp[2] == 1) {
    r <- c(r, x[length(x)])
    pp <- c(pp, 1)
  }
  r <- r[match(pp, probs)]
  if (names) 
    names(r) <- format(probs)
  if (weights) 
    attr(r, "weights") <- structure(w, dimnames = list(NULL, 
                                                       format(ps)))
  if (!se) 
    return(r)
}
