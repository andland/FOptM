fun <- function(X,  A) {
  G = -(A %*% X)
  f = 0.5 * sum(t(G) %*% X)
  list(f = f, G = G)
}

n = 1000; k = 6;
A = matrix(rnorm(n^2), n, n); A = t(A) %*% A;
opts = list()
opts$record = 1; #
opts$mxitr  = 1000;
opts$xtol = 1e-5;
opts$gtol = 1e-5;
opts$ftol = 1e-8;

X0 = matrix(rnorm(n * k), n, k)
X0 = qr.Q(qr(X0))

res = OptStiefelGBB(fun, X0, opts = opts, A = A)
ee = eigen(A)
sum(ee$values[1:6])
-2 * res$out$fval

fun(res$X, A)$f * -2
fun(ee$vectors[, 1:6], A)$f * -2
