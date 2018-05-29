

 single_HMC_sample <- function(U, grad_U, epsilon, L, current_q, start_p = NULL)
{
  q <- current_q
  if (is.null(start_p))
    p <- rnorm(length(q), 0, 1)  # independent standard normal variates
  else
    p <- start_p

  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # store results
  intermediate_q <- vector('list', L + 1)
  intermediate_p <- vector('list', L + 1)
  intermediate_q[[1]] <- q
  intermediate_p[[1]] <- p
  # Alternate full steps for position and momentum
  for (i in 1:L) {
    # Make a full step for the position
    q <- q + epsilon * p
    intermediate_q[[1 + i]] <- q
    intermediate_p[[1 + i]] <- p
    # Make a full step for the momentum, except at end of trajectory
    if (i != L)
      p <- p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p <- p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p ^ 2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p ^ 2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (log(runif(1)) < current_U - proposed_U + current_K - proposed_K) {
    accept <- TRUE #return (q)  # accept
  } else {
    #return (current_q)  # reject
    q <- current_q
    accept <- FALSE
  }
  list(new_q = q,
       intermediate_q = intermediate_q,
       intermediate_p = intermediate_p,
       accept = accept)
}

HMC_sampler <- function(U, grad_U, epsilon, L, T, start_q, start_p, force_reject = FALSE) {
  current_q <- start_q

  accept <- rep(0, T)
  intermediate_q <- vector(mode = 'list', length = T)
  saved_q <- intermediate_q

  for (t in seq_len(T)) {
    if (t == 1)
      output <- single_HMC_sample(U, grad_U, epsilon, L, current_q, start_p)
    else
      output <- single_HMC_sample(U, grad_U, epsilon, L, current_q)
    if (!force_reject)
      new_q <- output$new_q
    else
      new_q = current_q
    intermediate_q[[t]] <- output$intermediate_q
    accept[t] <- output$accept
    current_q <- new_q
    saved_q[[t]] <- current_q
  }
  intermediate <- purrr::map(intermediate_q, ~dplyr::bind_rows(purrr::map(., ~as.data.frame(t(.)))))
  out <- dplyr::bind_rows(purrr::map(seq_along(intermediate), ~dplyr::mutate(intermediate[[.]], draw = ., step = row_number() - 1)))
  out <- tibble::as_tibble(dplyr::transmute(out, x = V1, y = V2, draw, step, accept_path = accept[draw]))
  # out <- dplyr::bind_rows(tibble::data_frame(x = start_q[1], y = start_q[2], draw = 0, step = L, accept_path = 1), out)
  dplyr::mutate(out, accept_path = accept_path == 1,
                accept_point = accept_path & step == L,
                reject_point = !accept_path & step == L,
                draw = as.integer(draw),
                step = as.integer(step),
                path_start = step == 0,
                path_end = step == L)
}



#' animate_HMC
#'
#' @param b An object representing a banana-shaped distribution
#' @param epsilon HMC step size
#' @param n_steps Number of HMC leapfrog integration steps
#' @param n_samples Number of MCMC samples
#' @param start_q Starting position (parameter value)
#' @param start_p Starting momentum (random noise)
#' @param force_reject I can't remember what this does
#'
#' @return A list with a ggplot object and a data frame
#' @export
#'
animate_HMC <- function(b, epsilon, n_steps, n_samples, start_q, start_p, force_reject = FALSE) {
  hmc <- HMC_sampler(b$U, b$grad_U, epsilon, n_steps, n_samples, start_q, start_p, force_reject) %>%
    dplyr::bind_rows(., dplyr::mutate(dplyr::filter(., draw == 1 & step == 0),
                                      path_end = TRUE, accept_point = TRUE,
                                      accept_path = TRUE, step = n_steps,
                                      draw = 0)) %>%
    dplyr::arrange(draw, step)


  hmc2 <-
    purrr::map(seq_len(n_steps), ~dplyr::mutate(dplyr::filter(hmc, step <= .), superframe = draw, subframe = . + 1)) %>%
    dplyr::bind_rows()

  hmc3 <-
    purrr::map(seq_len(n_samples), ~dplyr::mutate(dplyr::filter(hmc, draw < . & path_end), superframe = ., step = 0)) %>%
    dplyr::bind_rows() %>%
    tidyr::crossing(subframe = seq_len(n_steps + 1))


  hmc <- dplyr::bind_rows(hmc2, hmc3) %>%
    dplyr::arrange(superframe, subframe) %>%
    dplyr::mutate(frame = superframe + subframe / (10^(1 + ceiling(log10(n_steps))))) %>%
    dplyr::arrange(frame, draw, step) %>%
    dplyr::group_by(frame) %>%
    dplyr::mutate(cl = max(step) - step) %>%
    dplyr::ungroup()



  g <- ggplot2::ggplot(hmc, ggplot2::aes(x = x, y = y, frame = frame)) +
    ggplot2::geom_path(ggplot2::aes(group = factor(draw):factor(frame), colour = cl), show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(shape = accept_path),
                        data = hmc %>% dplyr::filter(accept_point | reject_point),
                        show.legend = FALSE) +
    ggplot2::labs(x = expression(theta[1]), y = expression(theta[2])) +
    viridis::scale_colour_viridis(alpha = 1, begin = 0, end = 1, direction = 1, option = 'D', discrete = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::geom_blank(ggplot2::aes(shape = accept_path),
                        data = dplyr::mutate(hmc[1:2,], accept_path = c(FALSE, TRUE))) +
    ggplot2::scale_shape_manual(values = c(4, 21))
  list(gg = g, data = hmc)
}


#' get_par_df
#'
#' @param x a stanfit object
#'
#' @return a data.frame
#' @export
#'
get_par_df <- function(x) {
  map(seq_len(x@sim$chains), ~mutate(as.data.frame(x@sim$samples[[.x]]), Chain = .x)) %>%
    map(~filter(.x, row_number() > x@sim$warmup / x@sim$thin)) %>%
    bind_rows()
}



