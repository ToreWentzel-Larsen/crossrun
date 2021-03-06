#' Joint probabilities, n=14, symmetric case
#'
#' Joint probabilities of the number C of crossings
#' (0, ... 13) and the longest run L (1, ..., 14) in a
#' series of n=14 independent Bernoulli observations for the
#' symmetric case (success probability 0.5). The probabilities
#' are stored in the "times" representations, multiplied by
#' \eqn{2^{14-1}=8192} and are integers in the symmetric
#' case. Only the joint distributions for n=14, 60, 100
#' and success probabilities 0.5 and 0.6 are included in
#' the package to avoid excessive storage, but many more
#' cases may be generated by the function crossrunsymm.
#'
#' @format matrix, 14 rows and 14 columns
#' @source generated by the function crossrunsymm and
#' transformed from an Rmpfr array to a matrix
"joint14symm"

#'
#' Joint probabilities, n=14, success probability 0.6
#'
#' The joint probabilities of the number C og crossings
#' (0, ... 13) and the longest run L (1, ..., 14) in a
#' series of n=14 independent Bernoulli observations for
#' success probability 0.6. The probabilities are stored
#' in the "times" representations, multiplied by
#' \eqn{2^{14-1}=8192}. Only the joint distributions for
#' n=14, 60, 100 and success probabilities 0.5 and 0.6 are
#' included in the package to avoid excessive storage, but
#' many more cases may be generated by the function crossrunbin.
#'
#' @format matrix, 14 rows and 14 columns
#' @source generated by the function crossrunbin and
#' transformed from an Rmpfr array to a matrix
"joint14.6"

#' Joint probabilities, n=14, around the empirical median
#'
#' Joint probabilities of the number C of crossings
#' (1, ... 13) and the longest run L (1, ..., 17) in a
#' series of n=60 Bernoulli observations around its
#' empirical median. The probabilities  are stored in
#' the "times" representations, multiplied by
#' (60 by 30)/2, the number of constellations starting
#' above the median, and are integers. About the empirical
#' median there is at least one crossing, and the longest
#' run cannot exceed 14/2=7. Only the joint distributions
#' for n=14, 60 are included in the package to avoid excessive
#' storage, but many more cases may be generated by the function
#' 'crossrunem. Since these computations are demanding in terms
#' of storage and computation time, they are at present not
#' performed for n much above 60.
#'
#' @format matrix, 13 rows and 7 columns
#' @source generated by the function crossrunsymm and
#' transformed from an Rmpfr array to a matrix
"joint14em"

#' Joint probabilities, n=60, symmetric case
#'
#' The joint probabilities of the number C og crossings
#' (0, ... 59) and the longest run L (1, ..., 60) in a
#' series of n=60 independent Bernoulli observations for the
#' symmetric case (success probability 0.5). The probabilities
#' are stored in the "times" representations, multiplied by
#' \eqn{2^{60-1}} and are integers in the symmetric
#' case. Only the joint distributions for n=15, 60, 100
#' and success probabilities 0.5 and 0.6 are included in
#' the package to avoid excessive storage, but many more
#' cases may be generated by the function crossrunsymm.
#'
#' @format matrix, 60 rows and 60 columns
#' @source generated by the function crossrunsymm and
#' transformed from an Rmpfr array to a matrix
"joint60symm"

#' Joint probabilities, 60, success probability 0.6
#'
#' The joint probabilities of the number C og crossings
#' (0, ... 59) and the longest run L (1, ..., 60) in a
#' series of n=60 independent Bernoulli observations for
#' success probability 0.6. The probabilities are stored
#' in the "times" representations, multiplied by
#' \eqn{2^{60-1}}. Only the joint distributions for
#' n=15, 60, 100 and success probabilities 0.5 and 0.6 are
#' included in the package to avoid excessive storage, but
#' many more cases are generated in the script crossrun1.R.
#'
#' @format matrix, 60 rows and 60 columns
#' @source generated by the function crossrunbin and
#' transformed from an Rmpfr array to a matrix
"joint60.6"

#' Joint probabilities, n=100, symmetric case
#'
#' The joint probabilities of the number C og crossings
#' (0, ... 99) and the longest run L (1, ..., 100) in a
#' series of n=100 independent Bernoulli observations for the
#' symmetric case (success probability 0.5). The probabilities
#' are stored in the "times" representations, multiplied by
#' \eqn{2^{100-1}} and are integers in the symmetric
#' case. Only the joint distributions for n=15, 60, 100
#' and success probabilities 0.5 and 0.6 are included in
#' the package to avoid excessive storage, but many more
#' cases may be generated by the function crossrunsymm.
#'
#' @format matrix, 100 rows and 100 columns
#' @source generated by the function crossrunsymm and
#' transformed from an Rmpfr array to a matrix
"joint100symm"

#' Joint probabilities, n=100, success probability 0.6
#'
#' The joint probabilities of the number C og crossings
#' (0, ... 99) and the longest run L (1, ..., 100) in a
#' series of n=100 independent Bernoulli observations for
#' success probability 0.6. The probabilities are stored
#' in the "times" representations, multiplied by
#' \eqn{2^{100-1}}. Only the joint distributions for
#' n=15, 60, 100 and success probabilities 0.5 and 0.6 are
#' included in the package to avoid excessive storage, but
#' many more cases may be generated by the function crossrunbin.
#'
#' @format matrix, 100 rows and 100 columns
#' @source generated by the function crossrunbin and
#' transformed from an Rmpfr array to a matrix
"joint100.6"


#' Joint probabilities, n=60, around the empirical median
#'
#' Joint probabilities of the number C of crossings
#' (1, ... 59) and the longest run L (1, ..., 30) in a
#' series of n=14  Bernoulli observations around its
#' empirical median. The probabilities  are stored in
#' the "times" representations, multiplied by
#' (14 by 7)/2=1716, the number of constellations starting
#' above the median, and are integers. About the empirical
#' median there is at least one crossing, and the longest
#' runcannot exceed 60/2=30. Only the joint distributions
#' for n=14, 60 are included in the package to avoid excessive
#' storage, but many more cases may be generated by the function
#' 'crossrunem. Since these computations are demanding in terms
#' of storage and computation time, they are at present not
#' performed for n much above 60.
#' '#'
#' @format matrix, 59 rows and 30 columns
#' @source generated by the function crossrunem and
#' transformed from an Rmpfr array to a matrix
"joint60em"
