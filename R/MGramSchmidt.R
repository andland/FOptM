#' Gram Schmidt Orthogonormalization of a matrix
#'
#' @param V matrix to be orthogonormalized
#'
#' @return orthogonormalized matrix
#' @export
#'
#' @examples
#' M = matrix(rnorm(20), 10, 2)
#' MGramSchmidt(M)
MGramSchmidt <- function(V) {
  n = nrow(V)
  k = ncol(V)

  for (dj in 1:k) {
    # TODO: not sure seq_len is right
    for (di in seq_len(dj - 1)) {
      V[, dj] <- V[, dj] - proj(V[, di], V[, dj])
    }
    V[, dj] <- V[, dj] / norm(V[, dj, drop = FALSE], "F")
  }
  return(V)
}

#project v onto u
proj <- function(u, v) {
  v <- (crossprod(v, u) / crossprod(u, u)) * u
  return(v)
}
