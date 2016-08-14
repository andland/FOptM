#' Gradient based optimization on the Stiefel manifold
#'
#' @param fun function to be minimized. It must take \code{X} as the first argument and
#'   return a list with two named elements:
#'   \code{f}, which is the value of the objective function and \code{G}, which is
#'   the derivative
#' @param X the initial orthonormal
#' @param size if initial value for \code{X} is not given, the number of rows and columns
#'   can be given instead in a two-element vector. Number of columns must be less than the
#'   number of rows
#' @param opts a list of options for the optimization
#' @param ... additional argmuents to pass to \code{fun}
#'
#' @return list with \code{X}, which is the solution found, and \code{out},
#'   which includes supplementary material
#' @references Zaiwen Wen and Wotao Yin. A Feasible method for Optimization with
#'   Orthogonality Constraints, Optimization Online, 11/2010. Also as Rice
#'   CAAM Tech Report TR10-26.
#' @export
#'
#' @examples
#' fun <- function(X,  A) {
#'   G = -(A %*% X)
#'   f = 0.5 * sum(t(G) %*% X)
#'   list(f = f, G = G)
#' }
#'
#' n = 1000; k = 6
#' A = matrix(rnorm(n^2), n, n)
#' A = t(A) %*% A
#'
#' res = OptStiefelGBB(fun, size = c(n, k), A = A)
OptStiefelGBB <- function(fun, X, size, opts, ...) {
  #-------------------------------------------------------------------------
  # curvilinear search algorithm for optimization on Stiefel manifold
  #
  #   min F(X), S.t., X'*X = I_k, where X \in R^{n,k}
  #
  #   H = [G, X]*[X -G]'
  #   U = 0.5*tau*[G, X];    V = [X -G]
  #   X(tau) = X - 2*U * inv( I + V'*U ) * V'*X
  #
  #   -------------------------------------
  #   U = -[G,X];  V = [X -G];  VU = V'*U;
  #   X(tau) = X - tau*U * inv( I + 0.5*tau*VU ) * V'*X
  #
  #
  # Input:
  #           X --- n by k matrix such that X'*X = I
  #         fun --- objective function and its gradient:
  #                 [f, G] = fun(X,  data1, data2)
  #                 f, G are the objective function value and gradient, repectively
  #                 data1, data2 are addtional data, and can be more
  #                 Calling syntax:
  #                   [X, out]= OptStiefelGBB(X0, @fun, opts, data1, data2);
  #
  #        opts --- option structure with fields:
  #                 record = 0, no print out
  #                 mxitr       max number of iterations
  #                 xtol        stop control for ||X_k - X_{k-1}||
  #                 gtol        stop control for the projected gradient
  #                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
  #                             usually, max{xtol, gtol} > ftol
  #
  # Output:
  #           X --- solution
  #         Out --- output information
  #
  # -------------------------------------
  # For example, consider the eigenvalue problem F(X) = -0.5*Tr(X'*A*X);
  #
  # function demo
  #
  # function [f, G] = fun(X,  A)
  #   G = -(A*X);
  #   f = 0.5*sum(dot(G,X,1));
  # end
  #
  # n = 1000; k = 6;
  # A = randn(n); A = A'*A;
  # opts.record = 0; #
  # opts.mxitr  = 1000;
  # opts.xtol = 1e-5;
  # opts.gtol = 1e-5;
  # opts.ftol = 1e-8;
  #
  # X0 = randn(n,k);    X0 = orth(X0);
  # tic; [X, out]= OptStiefelGBB(X0, @fun, opts, A); tsolve = toc;
  # out.fval = -2*out.fval; # convert the function value to the sum of eigenvalues
  # fprintf('\nOptM: obj: #7.6e, itr: #d, nfe: #d, cpu: #f, norm(XT*X-I): #3.2e \n', ...
  #             out.fval, out.itr, out.nfe, tsolve, norm(X'*X - eye(k), 'fro') );
  #
  # end
  # -------------------------------------
  #
  # Reference:
  #  Z. Wen and W. Yin
  #  A feasible method for optimization with orthogonality constraints
  #
  # Author: Zaiwen Wen, Wotao Yin
  #   Version 1.0 .... 2010/10
  #-------------------------------------------------------------------------
  varargin <- list(...) # TODO: not sure if this is needed
  out = list()

  ## Size information
  if (missing(X)) {
    n = size[1]
    k = size[2]
    X = matrix(rnorm(n * k), n, k)
    X = qr.Q(qr(X))
  } else {
    n = nrow(X)
    k = ncol(X)
  }

  stopifnot(n > k)

  if (missing(opts)) {
    opts = list()
  }

  if (!is.null(opts[["xtol"]])){
    if (opts$xtol < 0 || opts$xtol > 1){
      opts$xtol <- 1e-6
    }
  } else {
    opts$xtol <- 1e-6
  }

  if (!is.null(opts[["gtol"]])){
    if (opts$gtol < 0 || opts$gtol > 1){
      opts$gtol <- 1e-6
    }
  } else {
    opts$gtol <- 1e-6
  }

  if (!is.null(opts[["ftol"]])){
    if (opts$ftol < 0 || opts$ftol > 1){
      opts$ftol <- 1e-12
    }
  } else {
    opts$ftol <- 1e-12
  }

  # parameters for control the linear approximation in line search
  if (!is.null(opts[["rho"]])){
    if (opts$rho < 0 || opts$rho > 1){
      opts$rho <- 1e-4
    }
  } else {
    opts$rho <- 1e-4
  }

  # factor for decreasing the step size in the backtracking line search
  if (!is.null(opts[["eta"]])){
    if (opts$eta < 0 || opts$eta > 1){
      opts$eta <- 0.1
    }
  } else {
    opts$eta <- 0.2
  }

  # parameters for updating C by HongChao, Zhang
  if (!is.null(opts[["gamma"]])){
    if (opts$gamma < 0 || opts$gamma > 1){
      opts$gamma <- 0.85
    }
  } else {
    opts$gamma <- 0.85
  }

  if (!is.null(opts[["tau"]])){
    if (opts$tau < 0 || opts$tau > 1e3){
      opts$tau <- 1e-3
    }
  } else {
    opts$tau <- 1e-3
  }

  # parameters for the  nonmontone line search by Raydan
  if (is.null(opts[["STPEPS"]])){
    opts$STPEPS <- 1e-10
  }

  if (!is.null(opts[["nt"]])){
    if (opts$nt < 0 || opts$nt > 100){
      opts$nt <- 5
    }
  } else {
    opts$nt <- 5
  }

  if (!is.null(opts[["projG"]])){
    if (!(opts$projG %in% c(1, 2))) {
      opts$projG <- 1
    }
  } else {
    opts$projG <- 1
  }

  if (!is.null(opts[["iscomplex"]])){
    if (!(opts$iscomplex %in% c(0, 1))) {
      opts$iscomplex <- 0
    }
  } else {
    opts$iscomplex <- 0
  }

  if (!is.null(opts[["mxitr"]])){
    if (opts$mxitr < 0 | opts$mxitr > 2^20){
      opts$mxitr <- 1000
    }
  } else {
    opts$mxitr <- 1000
  }

  if (is.null(opts[["record"]])){
    opts$record <- 0
  }


  #-------------------------------------------------------------------------------
  # copy parameters
  xtol <- opts$xtol
  gtol <- opts$gtol
  ftol <- opts$ftol
  rho  <- opts$rho
  STPEPS <- opts$STPEPS
  eta   <- opts$eta
  gamma <- opts$gamma
  iscomplex <- opts$iscomplex
  record <- opts$record

  nt <- opts$nt;
  crit <- matrix(NA, opts$mxitr, 3) # TODO: not sure

  invH <- TRUE
  if (k < n / 2) {
    invH <- FALSE
    eye2k <- diag(nrow = 2 * k)
  } else {
    eyen <- diag(nrow = n)
  }
  eyek <- diag(nrow = k)


  ## Initial function value and gradient
  # prepare for iterations
  fun_out = fun(X, ...)
  f = fun_out$f
  G = fun_out$G
  out$nfe <- 1
  GX <- t(G) %*% X

  if (invH) {
    GXT <- G %*% t(X)
    H <- 0.5 * (GXT - t(GXT))
    RX <- H %*% X
  } else {
    if (opts$projG == 1) {
      U <-  cbind(G, X)
      V <- cbind(X, -G)
      VU <- t(V) %*% U
    } else if (opts$projG == 2){
      GB <- G - 0.5 * X %*% (t(X) %*% G)
      U <- cbind(GB, X)
      V <- cbind(X, -GB)
      VU <- t(V) %*% U
    }
    #U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
    #VX = VU(:,k+1:end); #VX = V'*X;
    VX <- t(V) %*% X
  }
  dtX <- G - X %*% GX
  nrmG <- norm(dtX, "F")

  Q <- 1; Cval <- f;  tau <- opts$tau

  ## Print iteration header if debug == 1
  if (opts$record == 1){
    fid <- 1
    # TODO: fix the printing out
    cat(fid, '----------- Gradient Method with Line search ----------- \n')
    cat(fid, '%4s %8s %8s %10s %10s\n', 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff')
    #fprintf(fid, '#4d \t #3.2e \t #3.2e \t #5d \t #5d	\t #6d	\n', 0, 0, f, 0, 0, 0);
  }

  ## main iteration
  for (itr in 1:opts$mxitr) {
    XP <- X; FP <- f; GP <- G; dtXP <- dtX
    # scale step size

    nls <- 1; deriv <- rho * nrmG^2 #deriv
    while (TRUE) {
      # calculate G, f,
      if (invH) {
        X <- solve(eyen + tau*H, XP - tau*RX)
      } else {
        # print(tau)
        # print(cbind(eye2k + (0.5 * tau) * VU, NA, VX))
        aa <- solve(eye2k + (0.5 * tau) * VU, VX)
        X <- XP - U %*% (tau*aa)
      }
      #if norm(X'*X - eyek,'fro') > 1e-6; stop('X^T*X~=I'); end
      if (is.complex(X) & !iscomplex) {
        stop('X is complex')
      }

      fun_out = fun(X, ...)
      f = fun_out$f
      G = fun_out$G
      out$nfe <- out$nfe + 1

      if ((f <= Cval - tau * deriv) || nls >= 5){
        break
      }
      tau <- eta * tau
      nls <- nls + 1
    }

    GX <- t(G) %*% X
    if (invH){
      GXT <- G %*% t(X)
      H <- 0.5 * (GXT - t(GXT))
      RX <- H %*% X
    } else {
      if (opts$projG == 1){
        U <-  cbind(G, X)
        V <- cbind(X, -G)
        VU <- t(V) %*% U
      } else if (opts$projG == 2){
        GB <- G - 0.5 * X %*% (t(X) %*% G)
        U <- cbind(GB, X)
        V <- cbind(X, -GB)
        VU <- t(V) %*% U
      }
      #U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
      #VX = VU(:,k+1:end); # VX = V'*X;
      VX <- t(V) %*% X
    }
    dtX <- G - X %*% GX
    nrmG  <- norm(dtX, 'F')

    S <- X - XP; XDiff <- norm(S, 'F') / sqrt(n)
    tau <- opts$tau; FDiff <- abs(FP - f) / (abs(FP) + 1)

    if (iscomplex) {
      #Y = dtX - dtXP; SY = (sum(sum(real(conj(S).*Y))));
      Y <- dtX - dtXP; SY <- abs(sum(Conj(S) * Y))
      if ((itr %% 2) == 0) {
        tau <- sum(Conj(S) * S) / SY
      } else {
        tau <- SY / sum(Conj(Y) * Y)
      }
    } else {
      #Y = G - GP;     SY = abs(sum(sum(S.*Y)));
      Y <- dtX - dtXP; SY <- abs(sum(S * Y))
      #alpha = sum(sum(S.*S))/SY;
      #alpha = SY/sum(sum(Y.*Y));
      #alpha = max([sum(sum(S.*S))/SY, SY/sum(sum(Y.*Y))]);
      if ((itr %% 2) == 0) {
        tau <- sum(S * S) / SY
      } else {
        tau  <- SY / sum(Y^2)
      }

      # #Y = G - GP;
      # Y = dtX - dtXP;
      # YX = Y'*X;     SX = S'*X;
      # SY =  abs(sum(sum(S.*Y)) - 0.5*sum(sum(YX.*SX)) );
      # if mod(itr,2)==0;
      #     tau = SY/(sum(sum(S.*S))- 0.5*sum(sum(SX.*SX)));
      # else
      #     tau = (sum(sum(Y.*Y)) -0.5*sum(sum(YX.*YX)))/SY;
      # end

    }
    tau <- max(min(tau, 1e20), 1e-20)

    if (record >= 1){
      # TODO: make printing out better
      cat(itr, tau, f, nrmG, XDiff, FDiff, nls, "\n")
      #fprintf('#4d  #3.2e  #4.3e  #3.2e  #3.2e (#3.2e, #3.2e)\n', ...
      #    itr, tau, f, nrmG, XDiff, alpha1, alpha2);
    }

    # TODO: I'm not sure what's going on here
    crit[itr, ] <- c(nrmG, XDiff, FDiff)
    mcrit <- colMeans(crit[(itr - min(nt, itr) + 1):itr, , drop = FALSE])
    #if (XDiff < xtol && nrmG < gtol ) || FDiff < ftol
    #if (XDiff < xtol || nrmG < gtol ) || FDiff < ftol
    #if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol
    #if ( XDiff < xtol || FDiff < ftol ) || nrmG < gtol
    if ( (XDiff < xtol && FDiff < ftol ) || nrmG < gtol || all(mcrit[2:3] < 10 * c(xtol, ftol))){
      if (itr <= 2) {
        ftol <- 0.1*ftol
        xtol <- 0.1*xtol
        gtol <- 0.1*gtol
      } else {
        out$msg <- 'converge'
        break
      }
    }

    Qp <- Q; Q <- gamma * Qp + 1; Cval <- (gamma * Qp * Cval + f) / Q
  }

  if (itr >= opts$mxitr){
    out$msg <- 'exceed max iteration'
  }

  out$feasi <- norm(t(X) %*% X - eyek, 'F')
  if (out$feasi > 1e-13) {
    X <- MGramSchmidt(X)
    fun_out = fun(X, ...)
    f = fun_out$f
    G = fun_out$G
    out$nfe <- out$nfe + 1
    out$feasi <- norm(t(X) %*% X - eyek, 'F')
  }

  out$nrmG <- nrmG
  out$fval <- f
  out$itr <- itr
  out$tau = tau

  return(list(X = X, out = out))

}

