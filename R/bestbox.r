#' bestbox function.
#'
#' @description Function for box with lowest probability for the target shift,
#' among boxes with probability target or higher, for shift = 0. Details are
#' described i  https://doi.org/10.1371/journal.pone.0233920 .
#' @param pt0 joint probabilities inthe symmetric case.
#' @param pts joint probabilities for the target shift.
#' @param target minimal box probability.
#' @param n1 sequence length.
#' @param mult multiplier for joint probabilities.
#' @param prec mpft precision.
#' @return c1 and l1, minimal value of Cand maximal value of L for the box
bestbox <- function(pt0,
                    pts,
                    target = 0.925,
                    n1     = 100,
                    mult   = 2,
                    prec   = 120) {
  nill    <- Rmpfr::mpfr(0, prec)
  one     <- Rmpfr::mpfr(1, prec)
  two     <- Rmpfr::mpfr(2, prec)
  multm   <- Rmpfr::mpfr(mult, prec)
  targetm <- Rmpfr::mpfr(target, prec)
  targt   <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n    <- pt0[[n1]]
  ptsn    <- pts[[n1]]
  bpt0    <- boxprobt(pt0n)  # box probabilities for no shift
  bpttarg <- boxprobt(ptsn)  # box probabilities for target shift
  boxprt  <-
    two * (multm ^ (n1 - 1)) # initialize to impossible high value
  for (cc in 0:(n1 - 1))
    for (ll in 1:n1) {
      if (pt0n[cc + 1, ll] > nill &
          bpt0[cc + 1, ll] >= targt &
          bpttarg[cc + 1, ll] < boxprt) {
        c1 <- cc
        l1 <- ll
        boxprt <- bpttarg[cc + 1, ll]
      }
    }
  return(c(c1, l1))
} # end function bestbox
