
###############################################################################
# Functions for doing variational inference for LDA with counts
#
# Based on Algorithm 27.1 in Murphy's Machine Learning
###############################################################################

## ---- requirements ----
library("reshape2")
library("data.table")
library("plyr")
library("dplyr")

## ---- simulation functions ----
rdirichlet <- function(n, alpha) {
  gammas <- sapply(alpha, function(x) rgamma(n, x))
  gammas <- as.matrix(gammas, nrow = n)
  gammas / rowSums(gammas)
}

#' @importFrom data.table rbindlist setcolorder
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @examples
#' lda_data <- generate_data(100)
#' lda_data$theta
#' lda_data$beta
#' head(lda_data$Nkv)
#' head(lda_data$Nv)
generate_data <- function(n_documents, document_lengths = NULL, alpha = NULL,
                          gamma = NULL) {
  if (is.null(document_lengths)) {
    document_lengths <- rep(10000, n_documents)
  }
  if (is.null(alpha)) {
    alpha <- rep(1, 5) # topics prior
  }
  if (is.null(gamma)) {
    gamma <- rep(1, 200) # vocab prior
  }

  # Document and higher level parameters
  K <- length(alpha)
  V <- length(gamma)
  theta <- rdirichlet(n_documents, alpha)
  beta <- rdirichlet(K, gamma)

  # generate word counts for each topic
  Nkv <- vector(length = n_documents, mode = "list")
  for (i in seq_len(n_documents)) {
    Nkv[[i]] <- generate_word_counts(document_lengths[i], theta[i, ], beta)
    Nkv[[i]][, "document" := i]
  }
  Nkv <- Nkv %>%
    rbindlist() %>%
    setcolorder(c("document", "topic", "word", "value"))

  # observed data doesn't see topic labels
  Nv <- Nkv %>%
    group_by(document, word) %>%
    summarise(value = sum(value))

  list(
    "Nkv" = Nkv,
    "Nv" = Nv,
    "theta" = theta,
    "beta" = beta
  )
}

#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table
#' @importFrom dplyr filter
generate_word_counts <- function(document_length, theta, beta) {
  probs <- (theta %*% matrix(1, nrow = 1, ncol = ncol(beta))) * beta
  Nkv <- matrix(
    rmultinom(1, document_length, as.numeric(probs)),
    nrow = nrow(beta),
    )
  melt(Nkv, varnames = c("topic", "word")) %>%
    as.data.table() %>%
    filter(value != 0) %>%
    as.data.table()
}
