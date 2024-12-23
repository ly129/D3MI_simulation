rm(list = ls())

source("datagen.R")
library(mice)
library(MASS)
set.seed(1112)

uneven <- F
hetero <- "strong"
mis.var <- 1:2
out.cts <- T
n.sites <- 5
ss <- 600
n.sim <- 1000

# file.name <- paste0("results/PQL_", ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")

if (n.sites == 10) {
  if (uneven) {
    n_vec <- c(150, rep(50, 9))
  } else {
    n_vec <- rep(60, 10)
  }
} else {
  if (uneven) {
    n_vec <- c(200, rep(100, 4))
  } else {
    n_vec <- rep(120, 5)
  }
}

n_vec <- n_vec * ss/600

n_vec

cdata <- mdata <- pqlmi <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  cdata[[i]] <- datagen(n_vec = n_vec, hetero = hetero, out.cts = out.cts)
  mdata[[i]] <- misgen(cdata[[i]], parr = 0.15, mis.var = mis.var)
}

if (identical(mis.var, 1:2)) {
  folder <- "missing12/"
} else if (mis.var == 1) {
  folder <- "missing1/"
} else if (mis.var == 2) {
  folder <- "missing2/"
}


file.name <-  paste0("simdata/", folder, ifelse(uneven, "uneven", "even"), "_", hetero, "_n.sites", n.sites, "_ss", ss, ".Rdata")
save.image(file.name)
