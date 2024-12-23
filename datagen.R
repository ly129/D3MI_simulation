

# data generation 
datagen <- function(n_vec, hetero = "strong", out.cts = TRUE){
  ncl <- length(n_vec)
  
  ##### fixed effect
  alpha <- c(-0.5, 0.5, 0.5, 0.1) # old mub
  
  ##### random effects, reduce to random intercept only
  sig.uy <- switch(hetero, strong = 0.5, moderate = 0.3, weak = 0.1)
  sig.y <- 0.5
  
  dat <- NULL
  for (k in 1:ncl){
    n <- n_vec[k]
    x3 <- rnorm(n, 0, 0.5)
    
    u2 <- rnorm(1, 0, 0.1) # bi12
    x2 <- rnorm(n, u2 + 0.1 * x3, 0.1) # sigmax2 is eij2
    # cor(x2, x3)
    
    # mu <- rmvnorm(1, mean = mux, sigma = eta)
    u1 <- rnorm(1, 0, 0.1)
    # x1 <- rnorm(n, mu[1] + mu[2]*x1, sigmax1)
    x1.logit <- u1 + 0.1 * x2 + 0.1 * x3 #; hist(x1.logit)
    x1 <- rbinom(n, 1, 1/(1 + exp(-x1.logit))) #; table(x1)
    tempx <- cbind(x1, x2, x3)
    
    uy <- rnorm(1, 0, sig.uy)
    if (out.cts) {
      tempy <- cbind(1,tempx)%*%alpha + uy + rnorm(n, 0, sig.y)
    } else {
      logit <- cbind(1, tempx) %*% alpha + uy
      tempy <- rbinom(n, 1, 1/(1 + exp(-logit))) #logit link
    }
    
    dat <- rbind(dat, cbind(k, tempy, tempx))
  }
  colnames(dat) <- c("st", "y", "x1", "x2", "x3")
  summary(dat); by(dat[, "y"], INDICES = dat[, "st"], FUN = summary); # by(dat[, "x1"], INDICES = dat[, "st"], FUN = summary)
  return(dat)
}

misgen <- function(dat, parr = .15, mis.var = 1) {
  # gamma <- matrix(0, nrow = 5, ncol=2) 
  # gamma[1,] <- log(parr/(1 - parr))
  gamma <- matrix(0, nrow = 5, ncol = 2)
  gamma[2,2] <- 1; gamma[3,2] <- 1.5; gamma[5,2] <- -1
  gamma[1,2] <- log(parr/(1 - parr)) - gamma[2,2]*mean(dat[,"y"]) - gamma[3,2]*mean(dat[,"x1"]) - gamma[5, 2]*mean(dat[, "x3"])
  gamma[2,1] <- 1; gamma[4,1] <- 1.5; gamma[5,1] <- -1
  gamma[1,1] <- log(parr/(1 - parr)) - gamma[2,1]*mean(dat[,"y"]) - gamma[4,1]*mean(dat[,"x2"]) - gamma[5, 1]*mean(dat[, "x3"])
  
  mis.id <- matrix(NA, nrow = nrow(dat), ncol = length(mis.var))
  for (j in 1:length(mis.var)) {
    logit <- cbind(1, dat[,-1])%*%gamma[, mis.var]
    p <- 1/(1 + exp(-logit))
    
    mis.id[, j] <- rbinom(nrow(dat), size = 1, prob = p) == 1
  }
  for (j in 1:length(mis.var)) {
    dat[mis.id[, j], mis.var[j] + 2] <- NA # first col is cluster, second col is y
  }
  
  return(dat)
}


# uneven <- TRUE
# hetero <- "strong"
# mis.var <- 1
# out.cts <- TRUE
# 
# if (uneven) {
#   n_vec <- c(150, rep(50, 9))
# } else {
#   n_vec <- rep(60, 10)
# }
# 
# n.sim <- 1000
# cdata <- mdata <- mlmi <- vector(mode = "list", length = n.sim)
# for (i in 1:n.sim) {
#   cdata[[i]] <- datagen(n_vec = n_vec, hetero = hetero, out.cts = out.cts)
#   mdata[[i]] <- misgen(cdata[[i]], parr = 0.2, mis.var = mis.var)
# }

# summary(mdata[[1]])
# md.pattern(mdata[[1]])
# path <- paste0("/Users/yilian/Library/CloudStorage/Box-Box/ITCR/ITCR_dGLMM/code/simulation/simdata", )
# save(file = )


# cdata <- datagen(n_vec, hetero = "strong", out.cts = TRUE)
# mdata <- misgen(cdata, parr = 0.15, mis.var = 1:2); summary(mdata)
