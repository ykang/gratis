#' Specify parameters for a Mixture Autoregressive model
#'
#' This function allows the parameters of a mixture of k Gaussian AR processes
#' to be specified. The output is used in \code{\link{simulate.mar}()}.
#' The model is of the form
#' \deqn{y_t = \phi_{0,i} + \phi_{1,i}y_{t-1} + \cdots + \phi_{p,i}y_{t-p} + \sigma_{i,t}\epsilon_t}
#' with probability \eqn{\alpha_i}, where \eqn{\epsilon_t} is a N(0,1) variate.
#' If any argument is \code{NULL}, the corresponding parameters are randomly selected. 
#' The AR parameters may be non-stationary. When randomly selected, they are chosen
#' to be the pi weights of an ARIMA(p,d,0)(P,D,0) process where p is in \{0,1,2,3\},
#' d is in \{0,1,2\}, P is in \{0,1,2\} and D is in \{0,1\}. The model orders and the
#' parameters are uniformly sampled. The sigmas are uniformly sampled on (1,5)
#' and the weights are uniformly sampled on (0,1). The number of components is
#' uniformly sampled on \{1,2,3,4,5\}.
#' @param ar A (p+1) x k numeric matrix containing the AR parameters
#' (\eqn{\phi_{0,i}, \phi_{1,i},\dots,\phi_{p,i}}), \eqn{i=1,\dots,k} for each component.
#' @param sigmas A numeric vector of length k or a list of k GARCH specifications.
#' If it is a vector, it is assumed \eqn{\sigma_{i,t} = \sigma_i} and
#' \code{sigmas} = \eqn{\sigma_1,\dots,\sigma_k}.
#' If it is a list, each element should be the output from \code{fGarch::\link[fGarch]{garchSpec}()}.
#' @param weights A numeric vector of length k containing the probability of
#' each of the component processes, \eqn{\alpha_1,\dots,\alpha_k}.
#' @param frequency The number of seasonal periods (e.g., 12 for monthly data).
#' @return A `mar` object containing a list of \code{k}, \code{p}, \code{ar},
#' \code{sigmas} and \code{weights}.
#' @author Rob J Hyndman
#' @seealso \code{\link{generate_mar}}
#' @examples
#' n <- 100
#' # Seasonal MAR model with randomly selected parameters
#' mar_model(frequency = 4)
#' # MAR model with constant variances
#' # containing an AR(1) component and an AR(2) component
#' ar <- cbind(c(0, 0.8, 0), c(0, 0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model1 <- mar_model(ar = ar, sigmas = c(1, 2), weights = weights)
#'
#' # MAR model with heteroskedastic errors
#' sigmas.spec <- list(
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#' )
#' model2 <- mar_model(ar = ar, sigmas = sigmas.spec, weights = weights)
#' @export
mar_model <- function(ar = NULL, sigmas = NULL, weights = NULL, frequency = 1) {
  # How many components?
  if (is.null(ar) & is.null(sigmas) & is.null(weights)) {
    # choose a random number of components between 1 and 5
    k <- sample(seq(5), 1, replace = TRUE)
  } else if (!is.null(weights)) {
    k <- length(weights)
  } else if (!is.null(sigmas)) {
    k <- length(sigmas)
  } else {
    k <- NCOL(ar)
  }

  # weights
  if (is.null(weights)) {
    # Choose a random set of weights that sum to 1
    weights <- runif(k)
    weights <- weights / sum(weights)
  } else {
    if (length(weights) != k) {
      stop("Dimension of weights does not match other components")
    }
  }
  # sigmas
  if (is.null(sigmas)) {
    # Choose k variances in [1,5]
    sigmas <- runif(k, 1, 5)
  } else {
    if (length(sigmas) != k) {
      stop("Dimension of sigmas does not match other components")
    }
  }
  # AR parameters
  seasonal <- frequency > 1
  if (!is.null(ar)) {
    if (NCOL(ar) != k) {
      stop("Dimension of ar does not match other components")
    }
  } else {
    ar <- matrix(0, ncol = k, nrow = 4 * frequency + 4)
    p <- sample(c(0, 1, 2, 3), k, replace = TRUE)
    if (seasonal) {
      d <- sample(c(0, 1), k, replace = TRUE)
      D <- sample(c(0, 1), k, replace = TRUE)
      P <- sample(c(0, 1, 2), k, replace = TRUE)
    } else {
      D <- P <- rep(0,k)
      d <- sample(c(0, 1, 2), k, replace = TRUE)
    }
    for (i in seq(k)) {
      phi0 <- phi <- Phi <- 0
      if (d[i] + D[i] <= 1) 
        phi0 <- runif(1, -3, 3)
      if (p[i] > 0) 
        phi <- stationary_ar(p[i])
      if (seasonal & P[i] > 0) 
        Phi <- stationary_ar(P[i])
      pi_weights <- c(
        phi0,
        pi_coefficients(ar = phi, sar = Phi, d = d[i], D = D[i], m = frequency)
      )
      ar[seq_along(pi_weights), i] <- pi_weights
    }
  }
  # Trim redundant zeros
  nonzero <- which(rowSums(abs(ar)) > 0)
  ar <- ar[seq(max(nonzero)),,drop=FALSE]
  
  # Return mar object
  structure(list(ar = ar, sigmas = sigmas, weights = weights),
    class = "mar"
  )
}

# Function to return random coefficients from a stationary AR(p) process
stationary_ar <- function(p) {
  p <- as.integer(p)
  if(p < 1)
    stop("p must be a positive integer")
  else if(p == 1L)
    phi <- runif(1, -1, 1)
  else if(p == 2L) {
    phi2 <- runif(1, -1, 1)
    phi1 <- runif(1, phi2-1, 1-phi2)
    phi <- c(phi1, phi2)
  } else {
    unit_root <- TRUE
    while (unit_root) {
      phi <- runif(p,c(-3,-3,-1),c(3,1,1))
      roots <- polyroot(c(1, -phi))
      unit_root <- min(abs(roots)) < 1
    }
  }
  return(phi)
}
