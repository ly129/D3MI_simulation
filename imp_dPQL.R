source("DPQL.R")

mice.impute.2l.cts.dPQL <- function(y, ry, x, type,
                                    wy = NULL, intercept = TRUE, ...) {
  
  if (is.null(wy)) wy <- !ry
  if (intercept) {
    x <- cbind(1, as.matrix(x))
    type <- c(2, type)
    names(type)[1] <- colnames(x)[1] <- "(Intercept)"
  }
  
  clust <- names(type[type == -2])
  rande <- names(type[type == 2])
  fixe <- names(type[type > 0])
  
  lev <- unique(x[, clust])
  
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]
  
  fixedpart <- paste(
    "yobs ~ ", paste(fixe[-1L], collapse = "+")
  )
  
  randompart <- paste(
    "~ ", clust
  )
  
  suppressWarnings(fit <- try(
    dpql.fit(fixed = fixedpart,
             random = randompart,
             data = cbind(yobs, xobs),
             family = "gaussian",
             niter = 1),
    silent = TRUE
  ))
  if (!is.null(attr(fit, "class"))) {
    if (length(attr(fit, "class")) == 1) {
      if (attr(fit, "class") == "try-error") {
        warning("dPQL does not run. Simplify imputation model")
        return(y[wy])
      }
    }
  }
  
  # draw sigma*
  # sigmahat <- sigma(fit)
  # df <- nrow(fit@frame) - length(fit@beta)
  sigmahat <- fit$sigma
  df <- fit$N - length(fit$fixef)
  sigma2star <- df * sigmahat^2 / rchisq(1, df)
  
  # draw beta*
  # beta <- lme4::fixef(fit)
  # RX <- lme4::getME(fit, "RX")
  beta <- fit$fixef
  
  ##### sigmahat^2 * chol2inv(lme4::getME(fit, "RX")) is equal to vcov(fit) for lmer
  
  # cov-matrix, i.e., vcov(fit)
  # covmat <- sigma2star * chol2inv(RX)
  covmat <- fit$fix.vcov * sigma2star / sigmahat^2
  rv <- t(chol(covmat))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  
  # draw psi*
  # applying the standard Wishart prior
  # rancoef <- as.matrix(lme4::ranef(fit)[[1]])
  rancoef <- as.matrix(fit$ranef)
  lambda <- t(rancoef) %*% rancoef
  df.psi <- nrow(rancoef)
  temp.psi.star <- stats::rWishart(1, df.psi, diag(nrow(lambda)))[, , 1]
  temp <- MASS::ginv(lambda)
  ev <- eigen(temp)
  if (sum(ev$values > 0) == length(ev$values)) {
    deco <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
    psi.star <- MASS::ginv(deco %*% temp.psi.star %*% t(deco))
  } else {
    try(temp.svd <- svd(lambda))
    if (!inherits(temp.svd, "try-error")) {
      deco <- temp.svd$u %*% diag(sqrt(temp.svd$d), nrow = length(temp.svd$d))
      psi.star <- MASS::ginv(deco %*% temp.psi.star %*% t(deco))
    } else {
      psi.star <- temp
      warning("psi fixed to estimate")
    }
  }
  
  # Calculate myi, vyi and drawing bi per cluster
  for (jj in lev) {
    if (jj %in% unique(xobs[, clust])) {
      Xi <- Xobs[xobs[, clust] == jj, ]
      Zi <- as.matrix(Zobs[xobs[, clust] == jj, ])
      yi <- yobs[xobs[, clust] == jj]
      sigma2 <- diag(sigma2star, nrow = nrow(Zi))
      Mi <- psi.star %*% t(Zi) %*% MASS::ginv(Zi %*% psi.star %*% t(Zi) + sigma2)
      myi <- Mi %*% (yi - Xi %*% beta.star)
      vyi <- psi.star - Mi %*% Zi %*% psi.star
    } else {
      myi <- matrix(0, nrow = nrow(psi.star), ncol = 1)
      vyi <- psi.star
    }
    
    vyi <- vyi - upper.tri(vyi) * vyi + t(lower.tri(vyi) * vyi)
    # generating bi.star using eigenvalues
    deco1 <- eigen(vyi)
    if (sum(deco1$values > 0) == length(deco1$values)) {
      A <- deco1$vectors %*% sqrt(diag(deco1$values, nrow = length(deco1$values)))
      bi.star <- myi + A %*% rnorm(length(myi))
    } else {
      # generating bi.star using svd
      try(deco1 <- svd(vyi))
      if (!inherits(deco1, "try-error")) {
        A <- deco1$u %*% sqrt(diag(deco1$d, nrow = length(deco1$d)))
        bi.star <- myi + A %*% rnorm(length(myi))
      } else {
        bi.star <- myi
        warning("b_", jj, " fixed to estimate")
      }
    }
    
    # imputation
    y[wy & x[, clust] == jj] <- as.vector(
      as.matrix(X[wy & x[, clust] == jj, , drop = FALSE]) %*% beta.star +
        as.matrix(Z[wy & x[, clust] == jj, , drop = FALSE]) %*% as.matrix(bi.star) +
        rnorm(sum(wy & x[, clust] == jj)) * sqrt(sigma2star)
    )
  }
  y[wy]
}

mice.impute.2l.bin.dPQL <- function(y, ry, x, type,
                                    wy = NULL, intercept = TRUE, ...) {
  # mice::install.on.demand("lme4", ...)
  # mice::install.on.demand("MASS", ...)
  
  if (is.null(wy)) wy <- !ry
  if (intercept) {
    x <- cbind(1, as.matrix(x))
    type <- c(2, type)
    names(type)[1] <- colnames(x)[1] <- "(Intercept)"
  }
  
  clust <- names(type[type == -2])
  rande <- names(type[type == 2])
  fixe <- names(type[type > 0])
  
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  
  fixedpart <- paste(
    "yobs ~ ", paste(fixe[-1L], collapse = "+")
  )
  
  randompart <- paste(
    "~ ", clust
  )
  
  suppressWarnings(fit <- try(
    dpql.fit(fixed = fixedpart,
             random = randompart,
             data = cbind(yobs, xobs),
             family = "binomial",
             niter = 1),
    silent = TRUE
  ))
  if (!is.null(attr(fit, "class"))) {
    if (length(attr(fit, "class")) == 1) {
      if (attr(fit, "class") == "try-error") {
        warning("dPQL does not run. Simplify imputation model")
        return(y[wy])
      }
    }
  }
  
  # draw beta*
  # beta <- lme4::fixef(fit)
  # rv <- t(chol(vcov(fit)))
  beta <- fit$fixef
  rv <- t(chol(fit$fix.vcov))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  
  # calculate psi*
  # psi.hat <- matrix(lme4::VarCorr(fit)[[1L]],
  #                   nrow = dim(lme4::VarCorr(fit)[[1L]])[1L]
  # )
  ##### random intercept only, scalar
  psi.hat <- fit$ran.vcov
  
  s <- nrow(psi.hat) * psi.hat
  # rancoef <- as.matrix(lme4::ranef(fit)[[1L]])
  rancoef <- as.matrix(fit$ranef)
  lambda <- t(rancoef) %*% rancoef
  temp <- lambda + s
  if (attr(suppressWarnings(chol(temp, pivot = TRUE)), "rank") != nrow(temp)) {
    warning("The cov matrix is not full rank")
  }
  temp <- MASS::ginv(temp)
  ev <- eigen(temp)
  if (mode(ev$values) == "complex") {
    ev$values <- suppressWarnings(as.numeric(ev$values))
    ev$vectors <- suppressWarnings(matrix(as.numeric(ev$vectors),
                                          nrow = length(ev$values)
    ))
    warning("The cov matrix is complex")
  }
  if (sum(ev$values < 0) > 0) {
    ev$values[ev$values < 0] <- 0
    temp <- ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
  }
  deco <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
  temp.psi.star <- stats::rWishart(
    1,
    nrow(rancoef) + nrow(psi.hat),
    diag(nrow(psi.hat))
  )[, , 1L]
  psi.star <- MASS::ginv(deco %*% temp.psi.star %*% t(deco))
  
  # psi.star positive definite?
  if (!isSymmetric(psi.star)) psi.star <- (psi.star + t(psi.star)) / 2
  valprop <- eigen(psi.star)
  if (sum(valprop$values < 0) > 0) {
    valprop$values[valprop$values < 0] <- 0
    psi.star <- valprop$vectors %*% diag(valprop$values) %*% t(valprop$vectors)
  }
  
  # find clusters for which we need imputes
  clmis <- x[wy, clust]
  
  # the main imputation task
  for (i in clmis) {
    bi.star <- t(MASS::mvrnorm(
      n = 1L, mu = rep(0, nrow(psi.star)),
      Sigma = psi.star
    ))
    idx <- wy & (x[, clust] == i)
    logit <- X[idx, , drop = FALSE] %*% beta.star +
      Z[idx, , drop = FALSE] %*% matrix(bi.star, ncol = 1)
    vec <- rbinom(nrow(logit), 1, as.vector(1 / (1 + exp(-logit))))
    if (is.factor(y)) {
      vec <- factor(vec, c(0, 1), levels(y))
    }
    y[idx] <- vec
  }
  y[wy]
}


mice.df <- function(m, lambda, dfcom, method) {
  if (is.null(dfcom)) {
    dfcom <- 999999
    warning("Large sample assumed.")
  }
  lambda[lambda < 1e-04] <- 1e-04
  dfold <- (m - 1)/lambda^2
  dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
  df <- dfold * dfobs/(dfold + dfobs)
  if (method != "smallsample") 
    df <- dfold  ## Rubin 1987, 3.1.6, Van Buuren 2012, 2.30, added 31/10/2012
  return(df)
}

# Pooling functions
# Pool the parameters of substantive model
poolres.binom <- function(fit,...){
  m <- length(fit$analyses)
  dimf <- length(fit$analyses[[1]]$fixef)
  dimr <- 1 # length(as.data.frame(VarCorr(fit$analyses[[1]]))[,"vcov"])
  qfix <- matrix(NA, m, dimf)
  ufix <- matrix(0,  dimf, dimf)
  qran <- matrix(NA,  m, dimr)  
  for (i in 1:m){
    qfix[i,] <- fit$analyses[[i]]$fixef
    ufix <- ufix + as.matrix(fit$analyses[[i]]$fix.vcov)
    qran[i,] <- fit$analyses[[i]]$ran.vcov# t(chol(fit$analyses[[i]]$modelStruct$reStruct[[1]]))^2 # as.data.frame(VarCorr(fit$analyses[[i]]))[,"vcov"][1:dimr]
  }
  Qbarfix <- apply(qfix, 2, mean)
  Ubarfix <- ufix/m
  Bfix <- cov(qfix)
  Tfix <- Ubarfix + Bfix*(1 + 1/m)
  riv <- (1 + 1/m)*diag(Bfix/Ubarfix)
  lambda <- (1 + 1/m)*diag(Bfix/Tfix)
  df <- mice.df(m, lambda, fit$analyses[[1]]$N-dimf-dimr-1, method = "smallsample")
  fmi <- (riv + 2/(df + 3))/(riv + 1)
  resfix <- t(rbind(Qbarfix, sqrt(diag(Tfix)), Qbarfix/sqrt(diag(Tfix)), dnorm(abs(Qbarfix/sqrt(diag(Tfix)))),
                    df, riv, fmi))
  rownames(resfix) <- attr(fit$analyses[[1]]$coefficients$fixed, "names")
  colnames(resfix) <- c("Estimate", "Std. Error", "t value", "p-value", "df", "riv", "fmi")
  ranef <- apply(qran, 2, mean)
  list (fixef = resfix, ranef = ranef) # ranef has length 2: random intercept sd and residual sd
}

# Pool the parameters of models for x1 and x2
poolres.cts <- function(fit,...){
  m <- length(fit$analyses)
  dimf <- length(fit$analyses[[1]]$fixef)
  dimr <- 1 # length(as.data.frame(VarCorr(fit$analyses[[1]]))[,"vcov"])
  qfix <- matrix(NA, m, dimf)
  ufix <- matrix(0,  dimf, dimf)
  qran <- matrix(NA,  m, dimr)  
  for (i in 1:m){
    qfix[i,] <- fit$analyses[[i]]$fixef
    ufix <- ufix + as.matrix(fit$analyses[[i]]$fix.vcov)
    qran[i,] <- fit$analyses[[i]]$ran.vcov# t(chol(fit$analyses[[i]]$modelStruct$reStruct[[1]]))^2 # as.data.frame(VarCorr(fit$analyses[[i]]))[,"vcov"][1:dimr]
  }
  Qbarfix <- apply(qfix, 2, mean)
  Ubarfix <- ufix/m
  Bfix <- cov(qfix)
  Tfix <- Ubarfix + Bfix*(1 + 1/m)
  riv <- (1 + 1/m)*diag(Bfix/Ubarfix)
  lambda <- (1 + 1/m)*diag(Bfix/Tfix)
  df <- mice.df(m, lambda, fit$analyses[[1]]$N-dimf-dimr-1, method = "smallsample")
  fmi <- (riv + 2/(df + 3))/(riv + 1)
  resfix <- t(rbind(Qbarfix, sqrt(diag(Tfix)), Qbarfix/sqrt(diag(Tfix)), dnorm(abs(Qbarfix/sqrt(diag(Tfix)))),
                    df, riv, fmi))
  rownames(resfix) <- attr(fit$analyses[[1]]$coefficients$fixed, "names")
  colnames(resfix) <- c("Estimate", "Std. Error", "t value", "p-value", "df", "riv", "fmi")
  ranef <- apply(qran, 2, mean)
  list (fixef = resfix, ranef = ranef) # ranef has length 2:  random intercept sd and residual sd
}


# multiple imputation using random effects (mlmi)
mire.dpql <- function(dat,
                      fm,
                      method,
                      binom.out = FALSE,
                      MAX = 10, M = 5, ...){
  imp0 <- mice(dat, print = F, maxit=0)
  pred <- imp0$pred  
  # pred[pred == 1] <- 2 # otherwise, all remaining variables are also included as random effects in imputation models
  # 1: fixed effect, 2: random effect; 0: self; -2: cluster
  pred[pred[,"st"] != 0,"st"] <- -2
  imp <- mice(dat, print = F, pred=pred, method = method, maxit = MAX, m = M)
  # mice::densityplot(imp)
  
  fit <- list(analyses = NULL)
  fit$analyses <- vector(mode = "list", length = M)
  if (binom.out) {
    suppressWarnings(for (i in 1:M) {
      dd <- imp$data
      dd[is.na(dd$x1), which(names(imp$nmis) == "x1")] <- imp$imp$x1[, i]
      dd[is.na(dd$x2), which(names(imp$nmis) == "x2")] <- imp$imp$x2[, i]
      fit$analyses[[i]] <- dpql.fit(fixed = formula(fm),
                                    random = formula("~st"),
                                    data = dd,
                                    family = binomial(),
                                    niter = 1)
    })
    # fit <- with(imp, dpql(fixed = formula(fm),
    #                       random = formula(rm.dpql),
    #                       family = binomial(),
    #                       niter = 5)
  } else {
    suppressWarnings(for (i in 1:M) {
      dd <- imp$data
      dd[is.na(dd$x1), which(names(imp$nmis) == "x1")] <- imp$imp$x1[, i]
      dd[is.na(dd$x2), which(names(imp$nmis) == "x2")] <- imp$imp$x2[, i]
      fit$analyses[[i]] <- dpql.fit(fixed = formula(fm),
                                    random = formula("~st"),
                                    data = dd,
                                    family = gaussian(),
                                    niter = 1)
    })
  }
  res <- poolres.cts(fit)
  
  # if ("2l.bin.dPQL" %in% method) {list$niter <- }
  
  # # models of x1 and x2
  # fitx1 <- fitx2 <- list(analyses = NULL)
  # fitx1$analyses <- fitx2$analyses <- vector(mode = "list", length = M)
  # for (i in 1:M) {
  #   dd <- imp$data
  #   if (sum(is.na(dd$x1)) > 0) dd[is.na(dd$x1), which(names(imp$nmis) == "x1")] <- imp$imp$x1[, i]
  #   if (sum(is.na(dd$x2)) > 0) dd[is.na(dd$x2), which(names(imp$nmis) == "x2")] <- imp$imp$x2[, i]
  #   
  #   fitx1$analyses[[i]] <- dpql.fit(fixed = formula(x1 ~ x2 + x3),
  #                                   random = formula("~st"),
  #                                   data = dd,
  #                                   family = binomial(),
  #                                   niter = 100)
  # }
  # resx1 <- poolres.binom(fitx1)
  # 
  # for (i in 1:M) {
  #   dd <- imp$data
  #   # if (sum(is.na(dd$x1)) > 0) dd[is.na(dd$x1), which(imp$nmis != 0)] <- imp$imp$x1[, i]
  #   if (sum(is.na(dd$x2)) > 0) dd[is.na(dd$x2), which(names(imp$nmis) == "x2")] <- imp$imp$x2[, i]
  #   fitx2$analyses[[i]] <- dpql.fit(fixed = formula(x2 ~ x3),
  #                                   random = formula("~st"),
  #                                   data = dd,
  #                                   family = gaussian(),
  #                                   niter = 1)
  # }
  # resx2 <- poolres.cts(fitx2)
  
  
  list(fixef = res$fixef, ranef = res$ranef)#, 
       # fixefx1 = resx1$fixef, ranefx1 = resx1$ranef,
       # fixefx2 = resx2$fixef, ranefx2 = resx2$ranef)  
  
}


