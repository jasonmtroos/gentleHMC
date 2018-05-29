## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ----echo = FALSE, warning = FALSE, message = FALSE, include = FALSE, cache = FALSE----
library(rstan)
library(tidyverse)
library(brms)
library(rstanarm)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ----eval = FALSE--------------------------------------------------------
#  sm <- stan_model(model_code = ...)
#  fit <- sampling(sm)

## ----echo = FALSE, warning = FALSE, message = FALSE, include = FALSE-----
rstan_options(auto_write = TRUE)
readr::write_file(path = here::here('vignettes/part_1_files/sm1.stan'), x = 'parameters {
  real theta;
}
model {
  theta ~ normal(0, 1);
}\n\n')
sm <- stan_model(file = here::here('vignettes/part_1_files/sm1.stan'))
fit <- sampling(sm)

## ----fig.height = 4------------------------------------------------------
stan_trace(fit)

## ----echo = FALSE, include = FALSE---------------------------------------
sm2 <- '
data {
  int<lower = 0> n;
  int<lower = 0> p;
  vector[n] y;
  matrix[n, p] X;
}
parameters {
  real alpha;
  vector[p] beta;
  real<lower = 0> sigma;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ exponential(1);

  y ~ normal(alpha + X * beta, sigma);
}
'
readr::write_file(x = paste0(sm2, '\n\n'), path = here::here('vignettes/part_1_files/sm2.stan'))

## ----echo = FALSE, output = 'asis', comment=''---------------------------
cat(readr::read_file(here::here('vignettes/part_1_files/sm2.stan')))

## ----echo = FALSE, warning = FALSE, message = FALSE, include = FALSE-----
sm2 <- stan_model(here::here('vignettes/part_1_files/sm2.stan'))
set.seed(2)
n <- 100
p <- 4
alpha <- 1
beta <- rnorm(p, 0, 1)
sigma <- rexp(1)
X <- matrix(rnorm(n * p), ncol = p)
y <- as.vector(alpha + X %*% beta + rnorm(n, 0, sigma))
fit <- sampling(sm2, data = list(n = n, p = p, X = X, y = y))

## ----eval = FALSE--------------------------------------------------------
#  library(rstan)
#  sm <- stan_model(file = 'my_model.stan')
#  
#  X <- ...
#  y <- ...
#  d <- list(n = nrow(X), p = ncol(X),
#            X = X, y = y)
#  
#  fit <- sampling(sm, data = d)

## ----fig.height = 3, fig.width = 6---------------------------------------
stan_trace(fit)

## ----fig.height = 3, fig.width = 6---------------------------------------
stan_ac(fit)

## ----fig.height = 3, fig.width = 6, warning = FALSE, message = FALSE-----
stan_plot(fit)

## ----echo = FALSE, include = FALSE---------------------------------------
d <- data_frame(y) %>%
  bind_cols(as_data_frame(X)) %>%
  rename_at(vars(starts_with('V')), funs(str_replace(., 'V', 'X')))

## ----eval = FALSE--------------------------------------------------------
#  library(rstanarm)
#  fit <- stan_glm(y ~ 1 + X1 + X2 + X3 + X4, data = d,
#             prior = normal(0, 1),
#             prior_intercept = normal(0, 1),
#             prior_aux = exponential(1))

## ----echo = FALSE, include = FALSE---------------------------------------
fit <- stan_glm(y ~ 1 + X1 + X2 + X3 + X4, data = d,
           prior = normal(0, 1),
           prior_intercept = normal(0, 1),
           prior_aux = exponential(1))

## ----eval = FALSE--------------------------------------------------------
#  summary(fit)

## ----fig.height = 3, fig.width = 6, warning = FALSE, message = FALSE-----
fit %>% as.array() %>% bayesplot::mcmc_dens_overlay()

## ----fig.height = 3, fig.width = 6, warning = FALSE, message = FALSE-----
pp_check(fit)

## ----eval = FALSE--------------------------------------------------------
#  loo(fit)

