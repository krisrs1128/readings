#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Experiment using Variational Autoencoder
##
## Taken from https://tensorflow.rstudio.com/keras/articles/examples/variational_autoencoder.html
##
## author: sankaran.kris@gmail.com
## date: 12/20/2017

###############################################################################
## Setup
###############################################################################

library("keras")
K <- keras::backend()

params <- list(
  "batch_size" = 100L,
  "original_dim" = 784L,
  "latent_dim" = 2L,
  "intermediate_dim" = 256L,
  "epochs" = 50L,
  "epsilon_std" = 1.0
)

###############################################################################
## Model definition
###############################################################################
x <- layer_input(shape = c(params$original_dim))
h <- layer_dense(x, params$intermediate_dim, activation = "relu")
z_mean <- layer_dense(h, params$latent_dim)
z_log_var <- layer_dense(h, params$latent_dim)

## function to return variational lower bound
sampling <- function(arg) {
  z_mean <- arg[, 1:2]
  z_log_var <- arg[, 3:4]

  eps <- k_random_normal(
    shape = c(k_shape(z_mean)[[1]]),
    mean = 0,
    stddev = params$epsilon_std
  )

  z_mean + k_exp(z_log_var / 2) * eps
}

z <- layer_concatenate(
  list(z_mean, z_log_var)
) %>%
  layer_lambda(sampling)

decoder_h <- layer_dense(
  units = params$intermediate_dim,
  activation = "relu"
)
decoder_mean <- layer_dense(
  units = params$original_dim,
  activation = "sigmoid"
)

h_decoded <- decoder_h(z)
x_decoded_mean <- decoder_mean(h_decoded)

###############################################################################
## Write full encoder
###############################################################################
vae <- keras_model(x, x_decoded_mean)
encoder <- keras_model(x, z_mean) # inputs to latent space
decoder_input <- layer_input(shape = params$latent_dim)
h_decoded_2 <- decoder_h(decoder_input)
x_decoded_mean_2 <- decoder_mean(h_decoded_2)
generator <- keras_model(decoder_input, x_decoded_mean_2)

vae_loss <- function(x, x_decoded_mean) {
  cross_ent <- (params$original_dim / 1.0) * loss_binary_crossentropy(x, x_decoded_mean)
  kl_loss <- - 0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
  cross_ent + kl_loss
}

compile(vae, "rmsprop", vae_loss)

###############################################################################
## Apply this model to data
###############################################################################
mnist <- dataset_mnist()
x_train <- mnist$train$x/255
x_test <- mnist$test$x/255
x_train <- x_train %>% apply(1, as.numeric) %>% t()
x_test <- x_test %>% apply(1, as.numeric) %>% t()

fit(
  vae,
  x_train,
  x_train,
  epochs = params$epochs,
  batch_size = params$batch_size,
  validation_data = list(x_test, x_test)
)
