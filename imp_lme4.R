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
  dimf <- length(fixef(fit$analyses[[1]]))
  dimr <- length(as.data.frame(VarCorr(fit$analyses[[1]]))[,"vcov"])
  qfix <- matrix(NA, m, dimf)
  ufix <- matrix(0,  dimf, dimf)
  qran <- matrix(NA,  m, dimr)  
  for (i in 1:m){
    qfix[i,] <- fixef(fit$analyses[[i]])
    ufix <- ufix + as.matrix(vcov(fit$analyses[[i]]))
    qran[i,] <- as.data.frame(VarCorr(fit$analyses[[i]]))[,"vcov"]
  }
  Qbarfix <- apply(qfix, 2, mean)
  Ubarfix <- ufix/m
  Bfix <- cov(qfix)
  Tfix <- Ubarfix + Bfix*(1 + 1/m)
  riv <- (1 + 1/m)*diag(Bfix/Ubarfix)
  lambda <- (1 + 1/m)*diag(Bfix/Tfix)
  df <- mice.df(m, lambda, df.residual(fit$analyses[[1]]), method = "smallsample")
  fmi <- (riv + 2/(df + 3))/(riv + 1)
  resfix <- t(rbind(Qbarfix, sqrt(diag(Tfix)), Qbarfix/sqrt(diag(Tfix)), dnorm(abs(Qbarfix/sqrt(diag(Tfix)))),
                    df, riv, fmi))
  rownames(resfix) <- rownames(summary(fit$analyses[[1]])$coefficients)
  colnames(resfix) <- c(colnames(summary(fit$analyses[[1]])$coefficients), "df", "riv", "fmi")
  ranef <- apply(qran, 2, mean)
  list (fixef = resfix, ranef = ranef) # ranef has length 2: random intercept sd and residual sd
}

# Pool the parameters of models for x1 and x2
poolres.cts <- function(fit,...){
  m <- length(fit$analyses)
  dimf <- length(fixef(fit$analyses[[1]]))
  dimr <- length(as.data.frame(VarCorr(fit$analyses[[1]]))[,"vcov"])
  qfix <- matrix(NA, m, dimf)
  ufix <- matrix(0,  dimf, dimf)
  qran <- matrix(NA,  m, dimr)  
  for (i in 1:m){
    qfix[i,] <- fixef(fit$analyses[[i]])
    ufix <- ufix + as.matrix(vcov(fit$analyses[[i]]))
    qran[i,] <- as.data.frame(VarCorr(fit$analyses[[i]]))[,"vcov"][1:dimr]
  }
  Qbarfix <- apply(qfix, 2, mean)
  Ubarfix <- ufix/m
  Bfix <- cov(qfix)
  Tfix <- Ubarfix + Bfix*(1 + 1/m)
  riv <- (1 + 1/m)*diag(Bfix/Ubarfix)
  lambda <- (1 + 1/m)*diag(Bfix/Tfix)
  df <- mice.df(m, lambda, df.residual(fit$analyses[[1]]), method = "smallsample")
  fmi <- (riv + 2/(df + 3))/(riv + 1)
  resfix <- t(rbind(Qbarfix, sqrt(diag(Tfix)), Qbarfix/sqrt(diag(Tfix)), dnorm(abs(Qbarfix/sqrt(diag(Tfix)))),
                    df, riv, fmi))
  rownames(resfix) <- rownames(summary(fit$analyses[[1]])$coefficients)
  colnames(resfix) <- c(colnames(summary(fit$analyses[[1]])$coefficients), "p-value", "df", "riv", "fmi")
  ranef <- apply(qran, 2, mean)
  list (fixef = resfix, ranef = ranef) # ranef has length 2:  random intercept sd and residual sd
}


mire <- function(dat, randmodel, method,
                 binom.out = FALSE,
                 MAX = 10, M = 5, ...){
  imp0 <- mice(dat, print = F, maxit=0)
  pred <- imp0$pred  
  # pred[pred == 1] <- 2 # otherwise, all remaining variables are also included as random effects in imputation models
  # 1: fixed effect, 2: random effect; 0: self; -2: cluster
  pred[pred[,"st"] != 0,"st"] <- -2
  imp <- mice(dat, print = F, pred=pred, method = method, maxit = MAX, m = M)
  if (binom.out) {
    suppressWarnings(fit <- with(imp, glmer(formula(randmodel), family = binomial)))
    res <- poolres.binom(fit)
  } else {
    suppressWarnings(fit <- with(imp, lmer(formula(randmodel))))
    res <- poolres.cts(fit)
  }
  
  # if ("2l.bin.dPQL" %in% method) {list$niter <- }
  
  # # models of x1 and x2
  # suppressWarnings(fitx1 <- with(imp, glmer(x1 ~ x2 + x3 + (1|st), family = binomial)))
  # resx1 <- poolres.binom(fitx1)
  # suppressWarnings(fitx2 <- with(imp, lmer(x2 ~ x3 + (1|st))))
  # resx2 <- poolres.cts(fitx2)
  
  
  list(fixef = res$fixef, ranef = res$ranef)#, 
       # fixefx1 = resx1$fixef, ranefx1 = resx1$ranef,
       # fixefx2 = resx2$fixef, ranefx2 = resx2$ranef)  
}


