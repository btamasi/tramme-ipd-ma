## ========================================================================
## Supplementary code to demostrate the structure of the main simulation
## exercise of the article
##
## "Individual participant data meta-analysis with mixed-effects
##  transformation models" by Tamasi, Crowther, Puhan, Steyerberg, Hothorn
##
## IMPORTANT: The code does not exactly replicate the simulation results,
## because it is based on a subset of the original 3CIA dataset.
##
## 03/11/2021 by B. Tamasi
## ========================================================================
options(warn = -1)

library("survival")
library("tramME") ## >= 0.1.0
library("lme4")

## === Load the example dataset (see demo("IPD-MA", package = "tramME"))
load("3CIApart.RData")

## --- Calculate the scales for later use
sc <- sapply(data3CIA[, c("age", "fev1pp", "mmrc")], function(x) {
  diff(range(x, na.rm = FALSE))
})

## === Estimate Model 3
fit <- CoxphME(sui | ages + mmrcs + fev1pps ~ 1 + (ages + mmrcs + fev1pps | cohort),
  data = data3CIA, log_first = TRUE, order = 4, support = c(1, 150),
  control = optim_control(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-8))
stopifnot(fit$opt$convergence == 0)
re_str <- attr(fit$param, "re")
pr <- list(beta = coef(fit, with_baseline = TRUE), theta = varcov(fit, as.theta = TRUE))

## === Seeds to be used in the simulation
## NOTE: control the number of draws by setting this vector
seeds <- 1001:1500

## === Helper function to extract results
extract_fns <- function(object, t) {
  nd <- object$data[1, ]
  nd[] <- 1
  nd <- nd[rep(1, length(t)), ]
  nd[[variable.names(object, "response")]] <- t
  mod <- object$model$ctm
  ord <- get("order", envir = environment(object$model$ctm$bases$response))
  X <- model.matrix(mod, data = nd)
  cf <- coef(object, with_baseline = TRUE)
  nm <- names(cf)
  intidx <- grep("^Bs[1-9]*\\(", nm)
  nm <- unique(sub("^Bs.*\\:", "", nm[intidx]))
  idx <- split(intidx, rep(seq_along(nm), each = ord+1))
  out <- lapply(idx, function(i) {
    X_ <- X[, i, drop = FALSE]
    cf_ <- cf[i]
    drop(X_ %*% cf_)
  })
  names(out) <- nm
  out
}

## === Evaulate functions at t0
t0 <- seq(30, 150, by = 30)
tfn <- do.call("cbind", extract_fns(fit, t0))

## === Function to simulate from Model 3 and refit the Model
## -- tramME
simt <- function(seed) {
  ## --- Simulate from model
  sim <- simulate(fit, what = "joint", seed = seed)
  resp <- sim$responses
  re_true <- as.matrix(tramME:::.re_format(sim$ranef,
    re_str$termsize, re_str$names, re_str$blocksize, re_str$levels)[[1]])
  time2 <- system.time({
    data3CIA$ys <- resp
    mm2 <-  try(CoxphME(ys | ages + mmrcs + fev1pps ~ 1 + (ages + mmrcs + fev1pps | cohort),
      data = data3CIA, log_first = TRUE, order = 4, support = c(1, 150), initpar = pr,
      control = optim_control(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-8)))
      fe2 <- do.call("cbind", extract_fns(mm2, t0))
      re2 <- as.matrix(ranef(mm2)[[1]])
      vc2 <- varcov(mm2)[[1]]
    })
  if (inherits(mm2, "try-error") || mm2$opt$convergence != 0) {
    print(paste("Failed:", seed))
    return(list(fe = NA, re = NA, vc = NA, time = time2, seed = seed))
  } else {
    print(paste("Done:", seed))
    return(list(fe = fe2, re = re2, vc = vc2, time = time2, seed = seed))
  }
}

## -- Garcia et al (2019)
simg <- function(seed) {
  ## --- Simulate from model
  sim <- simulate(fit, what = "joint", seed = seed)
  resp <- sim$responses
  re_true <- as.matrix(tramME:::.re_format(sim$ranef,
    re_str$termsize, re_str$names, re_str$blocksize, re_str$levels)[[1]])
  ## --- Set up output
  fe1 <- matrix(NA, nrow = length(t0), ncol = 4)
  re1 <- array(NA, dim = c(dim(re_true), length(t0)))
  vc1 <- array(NA, dim = c(dim(varcov(fit)[[1]]), length(t0)))
  time1 <- system.time({
    start <- NULL
    for (i in seq_along(t0)) {
      data3CIA$yi <- ifelse(!is.na(resp[, 3]) & resp[, 3] < t0[i], 1, 0)
      mm1 <- try(
        glmer(yi ~ ages + mmrcs + fev1pps + (ages + mmrcs + fev1pps | cohort),
          data = data3CIA, family = binomial(link = "cloglog"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e4),
                                 check.conv.singular = "ignore"),
          start = start))
      if (inherits(mm1, "try-error") || mm1@optinfo$conv$opt != 0)
        break
      start <- list(theta = mm1@theta, fixef = fixef(mm1))
      fe1[i, ] <- fixef(mm1)
      re1[, , i] <- as.matrix(ranef(mm1)[[1]])
      vc1[, , i] <- as.matrix(VarCorr(mm1)[[1]])
    }
    colnames(fe1) <- names(fixef(mm1))
  })
  if (any(is.na(fe1))) {
    print(paste("Failed:", seed))
  } else {
    print(paste("Done:", seed))
  }
  return(list(fe = fe1, re = re1, vc = vc1, time = time1, seed = seed))
}

## === Simulations
## NOTE: This takes a while to run. Consider parallelizing these loops
sim_tramME <- lapply(seeds, function(ii) simt(ii))
sim_garcia <- lapply(seeds, function(ii) simg(ii))

## === Exctracting the results, creating figures
## --- helper function
rescale <- function(x, sc) {
  n <- ncol(x)
  nm <- sub("s$", "", colnames(x))
  for (i in 1:n) {
    s <- sc[nm[i]]
    if (is.na(s)) s <- 1
    x[, i] <- x[, i] / s
  }
  x
}

## --- Fixed effects
rres_garcia <- lapply(sim_garcia, function(x) {
  if (any(is.na(x$fe))) return(matrix(NA, ncol = 4, nrow = length(t0)))
  rescale(x$fe, sc)
})

rres_tramME <- lapply(sim_tramME, function(x) {
  if (any(is.na(x$fe))) return(matrix(NA, ncol = 4, nrow = length(t0)))
  rescale(x$fe, sc)
})

rtval <- rescale(tfn, sc)

## --- Random effects
re_str <- attr(fit$param, "re")
rre <- lapply(seq_along(sim_garcia), function(i) {
  seed <- sim_garcia[[i]]$seed
  out <- list()
  re <- simulate(fit, what = "ranef", seed = seed)
  re <- as.matrix(tramME:::.re_format(re[[1]], re_str$termsize, re_str$names,
                                      re_str$blocksize, re_str$levels)[[1]])
  out$true <- rescale(re, sc)
  if (any(is.na(tr <- sim_tramME[[i]]$re))) out$tramME <- NA * re
  else out$tramME <- rescale(tr, sc)
  if (any(is.na(ga <- sim_garcia[[i]]$re)))
    out$garcia <- array(NA, dim = c(dim(re), length(t0)))
  else {
    ga <- lapply(asplit(ga, 3), function(y) {
      colnames(y) <- colnames(re)
      rescale(y, sc)
    })
    out$garcia <- array(unlist(ga), c(dim(out$true), length(t0)))
  }
  out
})

mse <- list()
mse$tramME <- matrix(unlist(lapply(rre, function(x) {
  colMeans((x$true - x$tramME)^2)
})), ncol = length(rre))

mse$garcia <- matrix(unlist(lapply(rre, function(x) {
  ga <- lapply(asplit(x$garcia, 3), function(y) {
    (x$true - y)^2
  })
  ga <- array(unlist(ga), c(dim(x$true), length(t0)))
  apply(ga, 2, mean)
})), ncol = length(rre))

## --- RE covariance matrix
ss <- c(1, 1 / sc[c("age", "mmrc", "fev1pp")])
rvc_true <- diag(ss) %*% varcov(fit)[[1]] %*% diag(ss)

rvc_garcia <- lapply(sim_garcia, function(x) {
  ga <- lapply(asplit(x$vc, 3), function(xx) {
    diag(ss) %*% xx %*% diag(ss)
  })
  array(unlist(ga), c(dim(rvc_true), length(t0)))
})

rvc_tramME <- lapply(sim_tramME, function(x) {
  if (any(is.na(x$vc))) return(matrix(NA, nrow = 4, ncol = 4))
  diag(ss) %*% x$vc %*% diag(ss)
})

## --- Time
ti <- sapply(seq_along(sim_garcia), function(i) {
  c(sim_garcia[[i]]$time[3], sim_tramME[[i]]$time[3])
})

## --- Non-convergent cases
which(sapply(sim_garcia, function(x) any(is.na(x$fe))))
which(sapply(sim_tramME, function(x) any(is.na(x$fe))))

## --- Figure 3
pdf("Fig3.pdf", width = 9, height = 7.5)
n <- ncol(rres_garcia[[1]])
par(mfrow = c(ceiling(n/2), 2), mar = c(4, 5.5, 3, 1), cex = 0.8,
    las = 1, xaxs="i")
ylabs <- c(expression(hat(h)(t)), rep(expression(hat(beta)(t)), 3))
titles <- c("Log-cumulative baseline hazard",
            expression(Age), expression(mMRC), expression(FEV[1]))
plabs <- c("A", "B", "C", "D")
for (i in 1:n) {
  y1 <- sapply(rres_garcia, function(xx) {
    xx[, i]
  })
  y2 <- sapply(rres_tramME, function(xx) {
    xx[, i]
  })
  plot(0, type = "n", xlim = c(0, length(t0)),
       ylim = range(c(y1), c(y2), rtval[, i], na.rm = TRUE),
       xaxt = "n", ylab = ylabs[i], main = titles[i], xlab = "Time (in months)")
  grid(nx = 0, ny = NULL)
  abline(v = 1:length(t0), col = "lightgray", lty = 2)
  for (j in 0:(length(t0)-1)) {
    segments(x0 = j, x1 = j+1, y0 = rtval[j+1, i], lwd = 3)
  }
  c1 <- colorspace::qualitative_hcl(6, "Dark 3")[1]
  c2 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.1)[1]
  c3 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.7)[1]
  boxplot(t(y1), at = 1:length(t0) - 0.75, col = c2,
          medcol=c3, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1,
          lwd = 2, add = TRUE, boxwex = 0.3, xaxt = "n")
  c1 <- colorspace::qualitative_hcl(6, "Dark 3")[3]
  c2 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.1)[3]
  c3 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.7)[3]
  boxplot(t(y2), at = 1:length(t0) - 0.25, col = c2,
          medcol=c3, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1,
          lwd = 2, add = TRUE, boxwex = 0.3, xaxt = "n")
  axis(1, at = 1:length(t0) - 0.5, labels = paste("t =", t0), tck = 0)
  blx <- grconvertX(0.12, from = "nfc", to = "user")
  bly <- grconvertY(0.98, from = "nfc", to = "user")
  text(blx, bly, labels = plabs[i], xpd = TRUE, cex = 1.2)
  if (i == 1) {
    legend("bottomright", c("True value", "Garcia et al.", "Transformation model"),
           col = c(1, colorspace::qualitative_hcl(6, "Dark 3")[c(1, 3)]),
           lwd = c(3, 2, 2), bty = "n", cex = 0.9)
  }
}
dev.off()

## --- Figure S-8
pdf("FigS-8.pdf", width = 11, height = 9)
idx <- which(lower.tri(rvc_true, diag = TRUE), arr.ind = TRUE)
nm <- c("Frailty term", expression(Age), expression(mMRC), expression(FEV[1]))
par(mfrow = c(4, 4))
for (ii in 1:nrow(idx)) {
  mb <- if (idx[ii, ][1] == 4) 4 else 4
  ml <- if (idx[ii, ][2] == 1) 4.5 else 4.5
  par(mfg = idx[ii, ], mar = c(mb, ml, 2, 1), las = 1, cex = 0.6,
      xaxs = "i", xpd = FALSE)
  yy <- rvc_true[idx[ii, ][1], idx[ii, ][2]]
  yt <- sapply(rvc_tramME, function(x) {
    x[idx[ii, ][1], idx[ii, ][2]]
  })
  yg <- sapply(rvc_garcia, function(x) {
    x[idx[ii, ][1], idx[ii, ][2], ]
  })
  ylab <- if (idx[ii, ][2] == 1) nm[idx[ii, ][1]] else ""
  main <- if (idx[ii, ][1] == idx[ii, ][2]) nm[idx[ii, ][1]] else NULL
  plot(0, type = "n", xlim = c(0, length(t0) + 2),
       ylim = range(yt, yg, yy, na.rm = TRUE),
       xaxt = "n", xlab = "", ylab = ylab, main = main, cex.lab = 1.2)
  grid(nx = 0, ny = NULL)
  abline(h = yy, lwd = 3)
  c1 <- rep(colorspace::qualitative_hcl(6, "Dark 3")[c(1, 3)], c(length(t0), 1))
  c2 <- rep(colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.1)[c(1, 3)],
            c(length(t0), 1))
  c3 <- rep(colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.7)[c(1, 3)],
            c(length(t0), 1))
  boxplot(t(rbind(yg, yt)), col = c2,
          medcol=c3, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1,
          lwd = 2, add = TRUE, xaxt = "n")
  axis(1, at = 1:length(t0), labels = FALSE, tck = 0)
  aloc <- par("usr")[3] - diff(grconvertY(0:1, "inches", "user")) * 0.1
  text(1:(length(t0)+1), aloc, srt = 50, adj = 1,
       labels = c(paste("t =", t0), "tramME"), xpd = TRUE)
  par(mfg = c(1, 2))
  legend("topleft", inset = c(-0.2, 0),
         legend = c("True value", "Garcia et al.", "Transformation model"),
         col = c(1, colorspace::qualitative_hcl(6, "Dark 3")[c(1, 3)]),
         lwd = c(3, 2, 2), xpd = TRUE, cex = 1.4)
}
dev.off()

## --- Figure S-9
pdf("FigS-9.pdf", width = 7, height = 4)
par(mfrow = c(1, nrow(mse$tramME)), las = 1, cex = 0.8,
    mar = c(7.4, 4.1, 4.1, 1))
c1 <- colorspace::qualitative_hcl(6, "Dark 3")[c(1, 3)]
c2 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.1)[c(1, 3)]
c3 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.7)[c(1, 3)]
titles <- c("Frailty term",
            expression(Age), expression(mMRC), expression(FEV[1]))
for (ii in 1:nrow(mse$tramME)) {
  ylab <- if (ii == 1) "RMSE" else NULL
  boxplot(sqrt(cbind(mse$garcia[ii, ], mse$tramME[ii, ])), lwd = 2, xaxt = "n",
          col = c2, medcol=c3, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1,
          ylab = ylab, main = titles[ii])
  grid(nx = 0, ny = NULL)
  axis(1, at = 1:2, labels = FALSE)
  aloc <- par("usr")[3] - diff(grconvertY(0:1, "inches", "user")) * 0.1
  text(1:2, aloc, srt = 50, adj = 1,
       labels = c("Garcia et al.", "Transformation model"), xpd = TRUE)
}
dev.off()
