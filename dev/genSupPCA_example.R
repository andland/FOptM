sup_dim_red_log_like <- function(x, y, theta_x, theta_y, family_x, family_y, alpha) {
  -0.5 * genSupPCA::exp_fam_deviance(y, theta_y, family_y) - alpha * 0.5 * genSupPCA::exp_fam_deviance(x, theta_x, family_x)
}

sup_dim_red_deriv <- function(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, beta, U) {
  Ex = genSupPCA:::exp_fam_mean(theta_x, family_x)
  Ey = genSupPCA:::exp_fam_mean(theta_y, family_y)

  term1 = beta[-1] %*% crossprod(y - Ey, eta_centered)
  term2a = crossprod(x - Ex, eta_centered)
  term2 = t(U) %*% (term2a + t(term2a))
  term1 + alpha * term2
}


stiefel_objfun <- function(U, x, y, family_x, family_y, mu, eta_centered, beta, alpha) {
  theta_x = outer(rep(1, nrow(x)), mu) + eta_centered %*% tcrossprod(U)
  theta_y = cbind(1, eta_centered %*% U) %*% W$beta

  f <- sup_dim_red_log_like(x, y, theta_x, theta_y, family_x, family_y, alpha)
  gradient <- sup_dim_red_directional_deriv(W$x, W$y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, beta, U)

  list(f = value,
       G = gradient)
}
