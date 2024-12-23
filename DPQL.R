source("lmm.fit.R")

# fixed <- "x1~y+x2+x3"
# random <- "~st"
# data <- dat[complete.cases(dat), ]

dpql.fit <- function(fixed,
                     random,
                     data,
                     family = "binomial",
                     niter = 10L) {
  if (is.character(family)) {
    family <- get(family)
  }
  if (is.function(family)) {
    family <- family()
  }
  
  cluster <- attr(terms(x = formula(random), data = data), "term.labels")
  fixefid <- attr(terms(x = formula(fixed), data = data), "term.labels")
  if (attr(terms(x = formula(fixed), data = data), "intercept") == 1) intercept <- TRUE
  
  id.site.uniq <- unique(data[, cluster])
  K <- length(id.site.uniq)
  qz <- 1
  
  # step 1 initialize at the lead site to get beta0, ui0 = 0
  ranef <- rep(0, K)
  names(ranef) <- id.site.uniq
  
  # use glm to get init working Y from center cluster
  # coef.glm <- matrix(NA, nrow = K,
  #                    ncol = ifelse(intercept, length(fixefid) + 1, length(fixefid)))
  # rownames(coef.glm) <- id.site.uniq
  eta <- numeric(length = nrow(data))
  # for (cl in id.site.uniq) {
  #   data.cl <- subset(data, data[, cluster] == cl)
  #   if (intercept) {
  #     X <- cbind(1, as.matrix(data.cl[, fixefid, drop = FALSE])); colnames(X)[1] <- "Intercept"
  #   } else {
  #     X <- as.matrix(data.cl[, fixefid])
  #   }
  #   Y <- c(data.cl[, all.vars(formula(fixed))[1]])
  #   Z <- rep(1, length(Y))
  #   
  #   ## consider using the average from all sites
  #   fit.cl <- glm(Y ~ X - 1, family = family)
  #   fixef.cl <- fit.cl$coefficients
  #   
  #   coef.glm[cl, ] <- fixef.cl
  #   
  #   eta[data[, cluster] == cl] <- fit.cl$linear.predictors
  # }
  # ##### init fixef
  # fixef <- colMeans(coef.glm)
  fixef <- rep(0, ifelse(intercept, length(fixefid) + 1, length(fixefid)))
  
  # # # use lead site
  # # data1 <- subset(data, data[, cluster] == cluster1)
  # # X1 <- as.matrix(data1[, fixefid])
  # # Y1 <- c(data1[, 1])
  # # Z1 <- data1[, ranefid]
  # 
  # # fixef <- ifelse(is.na(fixef), 0, fixef)
  # # w <- fit0$prior.weights
  # eta <- fit0$linear.predictors
  # zz <- eta + fit0$residuals  
  # wz <- fit0$weights # working weights, final step of IWLS
  
  SiXYZ <- list() #vector(mode = "list", length = K)
  
  for (it in 1:niter) {
    etaold <- eta
    for (cl in id.site.uniq) {
      data.cl <- subset(data, data[, cluster] == cl)
      if (intercept) {
        X <- cbind(1, as.matrix(data.cl[, fixefid, drop = FALSE])); colnames(X)[1] <- "Intercept"
      } else {
        X <- as.matrix(data.cl[, fixefid])
      }
      Y <- c(data.cl[, all.vars(formula(fixed))[1]])
      Z <- rep(1, length(Y))
      
      eta.cl = c(X%*%fixef + ranef[cl])
      eta[data[, cluster] == cl] <- eta.cl
      mu = family$linkinv(eta.cl)
      mu.eta.val = family$mu.eta(eta.cl)
      wz = mu.eta.val^2/family$variance(mu) #* w
      zz = eta.cl + (Y-mu)/mu.eta.val #- offset  ## zz is working outcome
      
      # SiXYZ[[cl]] <- lmm.get.summary(Y = zz, X = X, Z = Z, weights = wz, id.site = id.site)
      SiXYZ[[cl]] <- list(SiX  = t(X*wz) %*% X,
                          SiXY = t(X*wz) %*% zz,
                          SiY  = sum(zz^2*wz),
                          SiXZ = t(X*wz) %*% Z, 
                          SiZ = sum(Z^2*wz),
                          SiZY = sum(Z*zz*wz),
                          ni = length(Y))
    }
    
    names(SiXYZ) <- id.site.uniq
    
    fit.dlmm.tmp <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=T, hessian=T)
    
    fixef <- fit.dlmm.tmp$b # ; fixef
    
    ranef <- unlist(fit.dlmm.tmp$ui)
    names(ranef) <- id.site.uniq
    
    result <- list()
    result$fixef <- c(fixef)
    result$fix.vcov <- solve(fit.dlmm.tmp$XvX)*fit.dlmm.tmp$s2
    result$ranef <- ranef
    if (family$family == "gaussian" ) {
      result$ran.vcov <- fit.dlmm.tmp$V
    } else {
      result$ran.vcov <- fit.dlmm.tmp$V/fit.dlmm.tmp$s2
    }
    result$niter <- it
    result$N <- nrow(data)
    if (family$family == "gaussian") result$sigma <- sqrt(fit.dlmm.tmp$s2)
    if (sum((eta - etaold)^2) < 1e-6 * sum(eta^2)) {
      break
    }
    if (it == niter) {warning("Maximum iteration reached.")}
  }
  return(result)
}