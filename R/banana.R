#' banana
#'
#' Let \eqn{(u, v)} be distributed bivariate normal with mean 0 and correlation \eqn{r}
#' and let \eqn{(x, y)} be distributed banana with \deqn{x = u a} and \deqn{y = v/a + b (u^2 + a^2).}
#'
#' @param a Input parameter \eqn{a} defined above
#' @param b Input parameter \eqn{b} defined above
#' @param r Input parameter \eqn{r} defined above
#'
#' @return list with functions \code{target}, \code{grad}, and \code{sample}
#' @export
#'
#' @details This function takes \eqn{a}, \eqn{b}, and \eqn{r} as inputs and returns a
#' list containing a function to evaluate the PDF of this distribution,
#' a function to evaluate the gradient of this distribution, and a
#' function to draw random samples from this distribution.
#' @examples
#' ban <- banana(a = 1.25, b = .5, r = .95)
#' xy <- ban$sample(1000)
#' plot(xy)
#' z <- ban$target(xy[,1], xy[,2])
#' # Show the density
#' library(ggplot2)
#' xyz <- as.data.frame(cbind(xy, z))
#' ggplot(xyz, aes(x, y, alpha = z)) + geom_point(colour = 'black')
#' # Simulate Hamiltonian behavior starting at 0 with step size .06
#' f <- function(g){ .06 * ban$grad(g[1], g[2])}
#' path <- matrix(0, ncol = 2, nrow = 40)
#' for (i in 1 + seq_len(39)) {
#'   path[i,] <- path[i-1,] + f(path[i-1,])
#' }
#' ggplot(NULL, aes(x, y)) +
#'   geom_point(aes(alpha = z), data = xyz, colour = 'black', show.legend = FALSE) +
#'   geom_path(aes(x = path[,1], y = path[,2], alpha = 1), data = NULL, colour ='red', show.legend = FALSE) +
#'   coord_cartesian(xlim = c(-1, .1), ylim = c(-.1, 1))
banana <-
  function(a, b, r) {
    PDF <-
      function(x, y) {
      exp((
        a ^ 8 * b ^ 2 - 2 * a ^ 6 * b * y - 2 * a ^ 2 * x * y * (b * x + r) +  x ^
          2 * (1 + b ^ 2 * x ^ 2 + 2 * b * x * r) + a ^ 4 * (2 * b ^ 2 * x ^ 2 + y ^
                                                               2 + 2 * b * x * r)
      ) / (2 * a ^ 2 * (-1 + r ^ 2))) / (2 * pi * sqrt(1 - r ^ 2))
      }
    U <-
      function(q) {
        x <- q[1]
        y <- q[2]

      -(((a^8*b^2 + x^2*(1 + b*x*(2*r + b*x)) - 2*a^6*b*y - 2*a^2*x*(r + b*x)*y +
    a^4*(2*b*x*(r + b*x) + y^2))/(a^2*(-1 + r^2)) - 2*log(pi) - log(4 - 4*r^2))/2)
      }
    dx <- function(x, y) {
      -(exp((a^8*b^2 + x^2*(1 + 2*b*r*x + b^2*x^2) - 2*a^6*b*y - 2*a^2*x*(r + b*x)*y +
      a^4*(2*b*r*x + 2*b^2*x^2 + y^2))/(2*a^2*(-1 + r^2)))*
   (a^4*b*(r + 2*b*x) + x*(1 + 3*b*r*x + 2*b^2*x^2) - a^2*(r + 2*b*x)*y))/
 (2*a^2*pi*(1 - r^2)^(3/2))
    }
    dy <- function(x, y) {
     (exp((a^8*b^2 + x^2*(1 + 2*b*r*x + b^2*x^2) - 2*a^6*b*y - 2*a^2*x*(r + b*x)*y +
     a^4*(2*b*r*x + 2*b^2*x^2 + y^2))/(2*a^2*(-1 + r^2)))*
  (a^4*b + x*(r + b*x) - a^2*y))/(2*pi*(1 - r^2)^(3/2))
    }
    d_log_target_x <- function(x, y) {
      (a^4*b*(r + 2*b*x) + x*(1 + 3*b*r*x + 2*b^2*x^2) - a^2*(r + 2*b*x)*y)/(a^2*(-1 + r^2))
    }
    d_log_target_y <- function(x, y) {
      (-(a^4*b) - x*(r + b*x) + a^2*y)/(-1 + r^2)
    }
    grad <-
      function(x, y) {
        c(d_target_x(x, y), d_target_y(x, y))
      }
    grad_U <-
      function(q) {
        x <- q[1]
        y <- q[2]
        -c(d_log_target_x(x, y), d_log_target_y(x, y))
      }
    sample <- function(n) {
      U <- chol(matrix(c(1, r, r, 1), ncol = 2))
      uv <- matrix(rnorm(2 * n, 0, 1), ncol = 2) %*% U
      u <- uv[,1]
      v <- uv[,2]
      x <- u * a
      y <- v/a + b * (u^2 + a^2)
      as.matrix(cbind(x, y))
    }
    list(PDF = PDF,
         grad = grad,
         U = U,
         grad_U = grad_U,
         sample = sample)
  }


