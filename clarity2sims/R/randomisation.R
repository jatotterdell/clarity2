# %%%%% Functions used for treatment assignment %%%%%

#' @title permuted_block_rand
#' @description Permuted block randomisation
#' @param target_alloc The target allocation ratio
#' @param sample_size The number of allocations to generate
#' @param blocksize The size of each block
#' @return A list
permuted_block_rand <- function(
  target_alloc,
  sample_size,
  blocksize
) {
  # Check that blocksize is integer ratio of weight sums
  # stopifnot()
  arms <- length(target_alloc)
  # Normalise weights
  w <- target_alloc / sum(target_alloc)
  R <- sum(target_alloc)
  L <- R / blocksize
  # Random numbers
  rand_num <- stats::runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation history
  n <- matrix(0, sample_size + 1, arms)
  # Constant history
  k <- numeric(sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)
  for(i in 1:sample_size) {
    k[i] <- floor((i - 1) / blocksize)
    p[i, ] <- (w * blocksize + w * blocksize * k[i] - n[i, ]) / (blocksize + blocksize * k[i] - (i - 1))
    trt[i] <- findInterval(rand_num[i], c(0, cumsum(p[i, ])))
    n[i + 1, ] <- n[i, ]
    n[i + 1, trt[i]] <- n[i, trt[i]] + 1
    d[i] <- sqrt(sum((n[i + 1, ] - i*w)^2))
  }
  return(list(
    imbalance = d,
    rand_num = rand_num,
    trt = trt,
    sample_size = n[-1, ],
    selection_prob = p))
}
