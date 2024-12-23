rm(list = ls())

source("imp_lme4.R")

library(mice)
library(lme4)

uneven <- T
hetero <- "strong"
mis.var <- 2
out.cts <- T
n.sites <- 10
ss <- 1200
# n.sim <- 1000

if (identical(mis.var, 1:2)) {
  folder <- "missing12/"
} else if (mis.var == 1) {
  folder <- "missing1/"
} else if (mis.var == 2) {
  folder <- "missing2/"
}


file.name <-  paste0("simdata/",folder, ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")

load(file.name)

# for (i in 1:n.sim)
loc <- c(which(colnames(cdata[[1]]) == "st" ) , which(colnames(cdata[[1]]) == "y"))
randmodel <- paste("y ~ ", paste(colnames(cdata[[1]])[-loc], collapse="+"),
                   "+ ( 1",
                   # "+ ( 1 +",
                   # paste(colnames(cdata)[-loc],collapse="+"),
                   "|", colnames(cdata[[1]])[loc[1]], ")")

method <- character(5)
if (1 %in% mis.var) method[3] <- "2l.bin"
if (2 %in% mis.var) method[4] <- "2l.lmer"

mlmi <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  cat("Iteration", i, "\n")
  mlmi[[i]] <- mire(dat = mdata[[i]],
                    randmodel,
                    method = method,
                    binom.out = !out.cts)
  mlmi[[i]]$fixef
}

mlmi.fix <- matrix(data = NA, nrow = 4, ncol = n.sim)

for (i in 1:n.sim) {
  mlmi.fix[, i] <- mlmi[[i]]$fixef[, 1]
}

# bias
(alpha.mean <- apply(mlmi.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
(se <- sqrt(sum((mlmi.fix - alpha.mean)^2)/n.sim))
# rmse
(rmse <- sqrt(sum((mlmi.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))

result.name <- paste0("results/", folder, "lme4_", ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")
save.image(result.name)
result.name
