rm(list = ls())

source("imp_PQL.R")

library(mice)
library(MASS)


uneven <- F
hetero <- "strong"
mis.var <- 1:2
out.cts <- T
n.sites <- 5
ss <- 600
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

fm <- "y~x1+x2+x3"

# i <- 20

method <- character(5)
if (1 %in% mis.var) method[3] <- "2l.bin.PQL"
if (2 %in% mis.var) method[4] <- "2l.cts.PQL"

pqlmi <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  cat("Iteration", i, ", ")
  pqlmi[[i]] <- mire.pql(dat = mdata[[i]],
                         fm = fm,
                         method = method,
                         binom.out = !out.cts,
                         MAX = 10,
                         M = 5)
  pqlmi[[i]]
}

pql.fix <- matrix(data = NA, nrow = 4, ncol = n.sim)

for (i in 1:n.sim) {
  pql.fix[, i] <- pqlmi[[i]]$fixef[, 1]
}

# bias
(alpha.mean <- apply(pql.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
(se <- sqrt(sum((pql.fix - alpha.mean)^2)/n.sim))
# rmse
(rmse <- sqrt(sum((pql.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))


# apply(pql.fix, MARGIN = 1, FUN = sd)
# save.image(file.name)
result.name <- paste0("results/", folder, "PQL_", ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")
save.image(result.name)
result.name