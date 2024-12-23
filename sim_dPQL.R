rm(list = ls())

source("imp_dPQL.R")

uneven <- T
hetero <- "strong"
mis.var <- 2
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

# dPQL with 1 iteration

file.name <-  paste0("simdata/", folder, ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")

load(file.name)

fm <- "y~x1+x2+x3"

method <- character(5)
if (1 %in% mis.var) method[3] <- "2l.bin.dPQL"
if (2 %in% mis.var) method[4] <- "2l.cts.dPQL"

dpqlmi <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  cat("Iteration", i, ", ")
  dpqlmi[[i]] <- mire.dpql(dat = mdata[[i]],
                           fm = fm,
                           method = method,
                           binom.out = !out.cts,
                           MAX = 10,
                           M = 5)
  dpqlmi[[i]]
}

dpql.fix <- matrix(data = NA, nrow = 4, ncol = n.sim)

for (i in 1:n.sim) {
  dpql.fix[, i] <- dpqlmi[[i]]$fixef[, 1]
}

# bias
(alpha.mean <- apply(dpql.fix, MARGIN = 1, FUN = mean))
(bias <- sqrt( sum( (alpha.mean - c(-0.5, 0.5, 0.5, 0.1))^2 ) ))

# se
(se <- sqrt(sum((dpql.fix - alpha.mean)^2)/n.sim))
# rmse
(rmse <- sqrt(sum((dpql.fix - c(-0.5, 0.5, 0.5, 0.1))^2)/n.sim))

result.name <- paste0("results/", folder, "dPQL_", ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")
save.image(result.name)
result.name
