## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ----cache = FALSE-------------------------------------------------------
library(gentleHMC)
library(tidyverse)
library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ------------------------------------------------------------------------
b <- banana(a = 1.25, b = .5, r = .95)
int_length <- 1.9
n_steps <- 14
n_samples <- 100
dpi <- 300

set.seed(100)
out <- animate_HMC(b, int_length/n_steps, n_steps, n_samples,
                               c(0, 4),
                               c(1, 1))

# gganimate::gganimate(out$gg, filename = here::here('vignettes/part_1_files/ani.mp4'),
#                      saver = 'mp4', title_frame = FALSE, interval = 1/n_steps,
#                      ani.height = 6/16*9, ani.width = 6,
#                      ani.dev = function(...) grDevices::png(..., res = dpi, units = "in", 
#                                                             type = 'cairo', bg = 'white'),
#                      other.opts = '-c:v libx264 -profile:v high -pix_fmt yuv420p -crf 1')

## ----message = FALSE-----------------------------------------------------
library(rstan)
sm <- stan_model(model_code = '
data {}
transformed data {
  real a;
  real b;
  real r;

  a = 1.25;
  b = .5;
  r = .95;
}
parameters {
  real x;
  real y;
}
model {
  target += (((a^8*b^2 + x^2*(1 + b*x*(2*r + b*x)) - 2*a^6*b*y - 2*a^2*x*(r + b*x)*y +
    a^4*(2*b*x*(r + b*x) + y^2))/(a^2*(-1 + r^2)) - 2*log(pi()) - log(4 - 4*r^2))/2);
}
           ')
fit <- sampling(sm, control = list(adapt_delta = .9))
ssp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
purrr::map(seq_along(ssp), ~dplyr::mutate(as.data.frame(ssp[[.]]), chain = .)) %>%
  dplyr::bind_rows() %>%
  tibble::as_tibble() %>%
  dplyr::summarise_at(dplyr::vars(stepsize__, n_leapfrog__), dplyr::funs(mean))
pd <- fit@sim$samples
pddf <- 
  purrr::map(seq_along(pd), ~dplyr::mutate(as.data.frame(pd[[.]]), Chain = .)) %>%
  dplyr::bind_rows()
pddf %>%
  dplyr::select(-lp__) %>%
  bayesplot::mcmc_scatter(alpha = .1) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = expression(theta[1]), y = expression(theta[2]))
pddf %>%
  bayesplot::mcmc_acf(pars = c('x', 'y'))
pddf %>%
  bayesplot::mcmc_combo(pars = c('x', 'y'))
pddf %>%
  bayesplot::mcmc_hex(pars = c('x', 'y'))

## ------------------------------------------------------------------------
n_steps <- 150
int_length <- 2*pi * (n_steps - 1)/n_steps
dpi <- 300

set.seed(100)

U <- function(q) {(1/2)*(q[1]^2 + q[2]^2) + log(2) + log(pi)}
grad_U <- function(q)q

track <- function(p = c(1, .5), q = c(.5, .5)) {
  out <- gentleHMC:::single_HMC_sample(U, grad_U, int_length/n_steps, n_steps,
                                       q,
                                       p)
  td <-
    bind_cols(out %>%
                pluck('intermediate_q') %>%
                map(~as.data.frame(t(.))) %>%
                bind_rows() %>%
                rename(x.q = V1, y.q = V2),
              out %>%
                pluck('intermediate_p') %>%
                map(~as.data.frame(t(.))) %>%
                bind_rows() %>%
                rename(x.p = V1, y.p = V2)) %>%
    mutate(z = sqrt(x.p^2 + y.p^2),
           step = row_number()) %>%
    as_tibble %>%
    select(-y.p, -x.p)
  td
}
td <- track()
bd <-   
  crossing(x.q = seq(-1, 1, by = .01)*2, y.q = seq(-1, 1, by = .01)*2) %>%
  mutate(z = exp(-.5 * (x.q^2 + y.q^2))) %>%
  crossing(step = seq_len(n_steps + 1))

g <- 
  ggplot(NULL, aes(x = x.q, y = y.q, frame = step)) +
  geom_tile(aes(fill = z), bd) +
  geom_point(aes(size = -z, alpha = -z), td, colour = 'black', fill = 'black') +
  scale_fill_gradient(low = 'white', high = gray(.25)) +
  scale_alpha_continuous(range = c(.5, 1)) +
  theme_minimal() +
  coord_cartesian(xlim = c(-1, 1)*2, ylim = c(-1, 1) * 9/16 * 2) +
  guides(alpha = 'none', fill = 'none', colour = 'none', size = 'none') +
  labs(x = expression(theta[1]), y = expression(theta[2])) 

# gganimate::gganimate(g, filename = here::here('vignettes/part_1_files/one_particle.mp4'), cumulative = TRUE,
#                      saver = 'mp4', title_frame = FALSE, interval = 1/16,
#                      ani.height = 6/16*9, ani.width = 6,
#                      ani.dev = function(...) grDevices::png(..., res = dpi, units = "in", 
#                                                             type = 'cairo', bg = 'white'),
#                      other.opts = '-c:v libx264 -profile:v high -pix_fmt yuv420p -crf 1')

## ------------------------------------------------------------------------
f <- function(t, q, p) {
  int_length <- pi/3
  n_steps <- 20
  n_samples <- 1
  dpi <- 300
  
  b <- list(U = function(q) {(1/2)*(q[1]^2 + q[2]^2) + log(2) + log(pi)}, grad_U = function(q)q)
  
  out <- gentleHMC:::animate_HMC(b, int_length/n_steps, n_steps, n_samples,
                                 q,
                                 p)
  gg <- out$gg + coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
  # gganimate::gganimate(gg, filename = here::here(str_glue('hmc_t_{t}.mp4')),
  #                      saver = 'mp4', title_frame = FALSE, interval = 1/n_steps,
  #                      ani.height = 6/16*9, ani.width = 6,
  #                      ani.dev = function(...) grDevices::png(..., res = dpi, units = "in",
  #                                                             type = 'cairo', bg = 'white'),
  #                      other.opts = '-c:v libx264 -profile:v high -pix_fmt yuv420p -crf 1')
}
f(1, c(0, 0), c(.5, -.5))
f(2, c(0.433, -0.433), c(-0.00997, 0.616))
f(3, c(0.208, 0.317), c(.5, .15))
f(4, c(0.537, 0.288), c(-.5, .15))

## ------------------------------------------------------------------------
int_length <- .99 * 2 * pi
n_steps <- 20
n_samples <- 10
dpi <- 300

b <- list(U = function(q) {(1/2)*(q[1]^2 + q[2]^2) + log(2) + log(pi)}, grad_U = function(q)q)

out <- animate_HMC(b, int_length/n_steps, n_steps, n_samples,
                               c(.5, .5),
                               c(-.5, .25))
gg <- out$gg + coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))
# gganimate::gganimate(gg, filename = here::here('vignettes/part_1_files/badL.mp4'),
#                      saver = 'mp4', title_frame = FALSE, interval = 1/n_steps,
#                      ani.height = 6/16*9, ani.width = 6,
#                      ani.dev = function(...) grDevices::png(..., res = dpi, units = "in", 
#                                                             type = 'cairo', bg = 'white'),
#                      other.opts = '-c:v libx264 -profile:v high -pix_fmt yuv420p -crf 1')


## ------------------------------------------------------------------------
b <- banana(a = 1.25, b = .5, r = .99)
int_length <- 1.9
n_steps <- 14
n_samples <- 50
dpi <- 300

set.seed(100)
out <- gentleHMC:::animate_HMC(b, .5 * int_length/n_steps, n_steps, n_samples,
                               start_q = c(1, 1),
                               start_p = NULL)

hmc <- out$data %>% filter(subframe > 1)

set.seed(100)
out <- gentleHMC:::animate_HMC(b, 1.0225 * int_length/n_steps, 1, n_steps * n_samples,
                               start_q = c(1, 1),
                               start_p = NULL) 
rw <- 
  out$data %>%
  filter(subframe > 1)

pd <-
  bind_rows(hmc %>% mutate(type = 'Hamiltonian'), 
            rw %>% mutate(type = 'Metropolis')) %>%
  group_by(type) %>%
  mutate(frame = dense_rank(frame)) %>%
  ungroup

gg <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = y, frame = frame)) +
  ggplot2::geom_point(ggplot2::aes(shape = accept_path),
                      data = pd %>% dplyr::filter(accept_point | reject_point),
                      show.legend = FALSE,
                      fill = '#eb3454',
                      colour = gray(.5)) +
  ggplot2::geom_path(ggplot2::aes(group = factor(draw):factor(frame), colour = cl), show.legend = FALSE) +
  ggplot2::labs(x = expression(theta[1]), y = expression(theta[2])) +
  viridis::scale_colour_viridis(alpha = 1, begin = 0, end = 1, direction = 1, option = 'D', discrete = FALSE) +
  ggplot2::facet_wrap(~type) + 
  ggplot2::theme_minimal() +
  ggplot2::geom_blank(ggplot2::aes(shape = accept_path),
                      data = dplyr::mutate(pd[1:2,], accept_path = c(FALSE, TRUE))) +
  ggplot2::scale_shape_manual(values = c(4, 21)) + 
  ggplot2::theme(strip.text = element_text(colour = '#eb3454', face = 'bold', size = rel(1.2)))


# gganimate::gganimate(gg, filename = here::here('vignettes/part_1_files/hmc.rw.mp4'),
#                      saver = 'mp4', title_frame = FALSE, interval = 2/n_steps,
#                      ani.height = 6/16*9, ani.width = 6,
#                      ani.dev = function(...) grDevices::png(..., res = dpi, units = "in",
#                                                             type = 'cairo', bg = 'white'),
#                      other.opts = '-c:v libx264 -profile:v high -pix_fmt yuv420p -crf 1')

## ----message = FALSE-----------------------------------------------------
sm <- stan_model(model_code = '
functions {
  real pc_exponential_rate(real U, real a) {return - log(a)/U;}
}
data {
  int<lower = 0> nr;
  int<lower = 0> nc;
  matrix[nr, nc] y;
}
transformed data {
  real sigma_rate; 
  sigma_rate = pc_exponential_rate(.5, .1);
}
parameters {
  vector[nr] theta_;
  vector<lower = 0>[nr] sigma_;
}
transformed parameters {
  vector[nr] theta;
  vector[nr] sigma;
  sigma = sigma_ / sigma_rate;
  theta = theta_ + 10;
}
model {
  theta_ ~ normal(0, 1);
  sigma_ ~ exponential(1);
  for (i in 1:nr) {
    y[i,] ~ normal(theta[i], sigma[i] / fabs(theta[i]));
  }
}
')
nc <- 30
nr <- 1
sigma <- matrix(rexp(nr) / (- log(.1)/.5), nrow = nr, ncol = nc)
theta <- matrix(10 + rnorm(nr), nrow = nr, ncol = nc)
y <- matrix(rnorm(nr * nc), nrow = nr) * sigma / abs(theta) + theta

fit <- sampling(sm, data = list(y = y, nc = nc, nr = nr), control = list(adapt_delta = .9),
                chains = 1)
fit2 <- sampling(sm, data = list(y = y, nc = nc, nr = nr), 
                 algorithm = 'HMC',
                 control = with(list(e = 1/1000), 
                                list(adapt_engaged = TRUE, 
                                     stepsize = e,
                                     int_time = 1e-64,
                                     metric = 'unit_e')),
                 iter = 20000, thin = 10, chains = 1)


ggplot_the_par_df <- function(fit) {
  get_par_df(fit) %>%
    mutate(Iteration = row_number()) %>%
    select(Iteration, `1` = sigma.1., `2` = theta.1.) %>%
    gather(par, val, -Iteration) %>%
    ggplot(aes(x = Iteration, y = val)) + 
    geom_line(colour = '#eb3454') +
    facet_wrap(~par, scale = 'free_y', ncol = 1, labeller = label_bquote(cols = theta[.(par)])) +
    labs(y = NULL)
}
ggplot_the_par_df(fit)
# ggsave(here::here('vignettes/part_1_files/trace_good.png'), width = 6, height = 3.5)
ggplot_the_par_df(fit2)
# ggsave(here::here('vignettes/part_1_files/trace_bad.png'), width = 6, height = 3.5)

