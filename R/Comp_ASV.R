#' I am not sure what this does
#'
#' @param sigma2_hat ???
#' @param BI ???
#' @param R ???
#' @param Gamma_hat estimated matrix
#' @param beta_hat matrix of the estimated values of \est{\beta}
#' @param Omega_hat estimated matrix
#' @param I_hat estimated matrix
#' @param I_ind_hat estimated matrix
#' @return ???

Comp_ASV <- function(sigma2_hat, BI, R, Gamma_hat, beta_hat, Omega_hat, I_hat, I_ind_hat) {
  D_tau_bar <- t(BI) %*% diag(Gamma_hat[I_hat]) %*% BI
  V1 <- as.numeric(sigma2_hat + t(beta_hat) %*% D_tau_bar %*% beta_hat)
  Q <- t(solve(crossprod(R), t(R)))
  V2 <- Omega_hat + t(Q) %*% diag(Gamma_hat) %*% Q

  K_hat <- length(I_ind_hat)

  D <- c()
  ms <- c()
  for (a in 1:K_hat) {
    group_a <- unlist(I_ind_hat[[a]])
    m_a <- length(group_a)
    ms[a] <- m_a

    D1 <- (2 * m_a - 1) * (D_tau_bar[a,a] ^ 2) * (m_a ^ 2) - sum(Gamma_hat[group_a] ** 2)
    D[a] <- D1 * (beta_hat[a] ^ 2) / m_a / ((m_a - 1) ^ 2)
  }

  V3 <- t(Q[I_hat,]) %*% diag(rep(D, ms)) %*% Q[I_hat,]
  return(diag(V1 * V2 + V3))
}
