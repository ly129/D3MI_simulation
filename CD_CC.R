source("datagen.R")
set.seed(1112)
# library(mice)

uneven <- F
hetero <- "strong"
mis.var <- 2
out.cts <- T
n.sites <- 5
ss <- 600
n.sim <- 1000

# method <- character(5)
# if (1 %in% mis.var) method[3] <- "2l.bin.dPQL"
# if (2 %in% mis.var) method[4] <- "2l.cts.dPQL"
# 
# cdata <- mdata <- cd <- cc <- vector(mode = "list", length = n.sim)
# for (i in 1:n.sim) {
#   cdata[[i]] <- datagen(n_vec = n_vec, hetero = hetero, out.cts = out.cts)
#   mdata[[i]] <- misgen(cdata[[i]], parr = 0.15, mis.var = mis.var)
# }
# 
# meany <- matrix(nrow = n.sim, ncol = n.sites)
# for (i in 1:n.sim) {
#   meany[i,] <- aggregate(y~st, data = cdata[[i]], FUN = mean)[, 2]
# }
# apply(meany, 2, max)
# apply(meany, 2, min)
# 
# meany.m <- matrix(nrow = n.sim, ncol = n.sites)
# for (i in 1:n.sim) {
#   meany.m[i,] <- aggregate(y~st, data = mdata[[i]][complete.cases(mdata[[i]]), ], FUN = mean)[, 2]
# }
# apply(meany.m, 2, max)
# apply(meany.m, 2, min)


if (identical(mis.var, 1:2)) {
  folder <- "missing12/"
} else if (mis.var == 1) {
  folder <- "missing1/"
} else if (mis.var == 2) {
  folder <- "missing2/"
}

# dPQL with 1 iteration

file.name <-  paste0("simdata/", folder, ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")

load(file.name)

loc <- c(which(colnames(cdata[[1]]) == "st" ) , which(colnames(cdata[[1]]) == "y"))
randmodel <- paste("y ~ ", paste(colnames(cdata[[1]])[-loc], collapse="+"),
                   "+ ( 1",
                   # "+ ( 1 +",
                   # paste(colnames(cdata)[-loc],collapse="+"),
                   "|", colnames(cdata[[1]])[loc[1]], ")")

cd <- cc <- vector(mode = "list", length = n.sim)

if (out.cts) {
  for (i in 1:n.sim) {
    # cat("Iteration", i, "\n")
    cd[[i]] <- lme4::lmer(randmodel,
                          data = as.data.frame(cdata[[i]]))
    
    cc[[i]] <- lme4::lmer(randmodel,
                          data = as.data.frame(mdata[[i]]))
  }
} else {
  for (i in 1:n.sim) {
    # cat("Iteration", i, "\n")
    cd[[i]] <- lme4::glmer(randmodel,
                           data = as.data.frame(cdata[[i]]),
                           family = binomial)
    cc[[i]] <- lme4::glmer(randmodel,
                           data = as.data.frame(mdata[[i]]),
                           family = binomial)
  }
}

cd.fix <- cc.fix <- matrix(data = NA, nrow = 4, ncol = n.sim)

for (i in 1:n.sim) {
  cd.fix[, i] <- lme4::fixef(cd[[i]])
  cc.fix[, i] <- lme4::fixef(cc[[i]])
}





cglm <- vector(mode = "list", length = n.sim)
method <- character(5)
if (1 %in% mis.var) method[3] <- "logreg"
if (2 %in% mis.var) method[4] <- "norm"
cglm.fix <- matrix(data = NA, nrow = 4, ncol = n.sim)

for (i in 1:n.sim) {
  # cat("Iteration", i, "\n")
  imp <- mice(mdata[[i]], m = 5, maxit = 10, method = method, printFlag = FALSE)
  
  if (out.cts) {
    cglm[[i]] <- with(imp, lm(y ~ x1 + x2 + x3))
  } else {
    cglm[[i]] <- with(imp, glm(y ~ x1 + x2 + x3, family = binomial()))
  }
  
  cglm.fix[, i] <- summary(pool(cglm[[i]]))$estimate
}





# bias cd
(alpha.mean <- apply(cd.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
(se <- sqrt(sum((cd.fix - alpha.mean)^2)/n.sim))
# rmse
(rmse <- sqrt(sum((cd.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))

# bias cc
(alpha.mean <- apply(cc.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
(se <- sqrt(sum((cc.fix - alpha.mean)^2)/n.sim))
# rmse
(rmse <- sqrt(sum((cc.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))

# bias
(alpha.mean <- apply(cglm.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
( se <- sqrt( sum( (cglm.fix - alpha.mean)^2 )/n.sim ) )
# rmse
(rmse <- sqrt(sum((cglm.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))

result.name <- paste0("results/", folder, "dPQL_", ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")
# save.image(result.name)
result.name
