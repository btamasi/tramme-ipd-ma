## ========================================================================
## Supplementary code to demostrate the structure of the secondary
## simulation exercise (the effect of misspecification) of the article
##
## "Individual participant data meta-analysis with mixed-effects
##  transformation models" by Tamasi, Crowther, Puhan, Steyerberg, Hothorn
##
## IMPORTANT: The code does not exactly replicate the simulation results,
## because it is based on a subset of the original 3CIA dataset.
##
## 03/11/2021 by B. Tamasi
## ========================================================================
library("survival")
library("tramME") ## 0.1.2

## === Load the example dataset (see demo("IPD-MA", package = "tramME"))
load("3CIApart.RData")

## --- Calculate the scales for later use
sc <- sapply(data3CIA[, c("age", "fev1pp", "mmrc")], function(x) {
  diff(range(x, na.rm = FALSE))
})


## === Estimate Model 1
fit <- CoxphME(sui ~ ages + mmrcs + fev1pps + (ages + mmrcs + fev1pps | cohort),
  data = data3CIA, log_first = TRUE, order = 4, support = c(1, 150),
  control = optim_control(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-8))
stopifnot(fit$opt$convergence == 0)

## Extract expected probabilistic index (integrating out the random effects)
E_PI <- function(mod, var = variable.names(mod, "shift"), scale = NULL) {
  if (is.null(scale)) {
    scale <- rep(1, length(var))
    names(scale) <- var
  }
  fun <- switch(mod$model$ctm$todistr$name,
                "minimum extreme value" = function(b) 1 / (1 + exp(b)),
                "logistic" = function(b) 1 - tram::PI(logOR = b),
                stop("Unkown model!"))
  fun_w <- function(re, b, sd, sc) fun((b + re) / sc) * dnorm(re, mean = 0, sd = sd)
  out <- sapply(var, function(v) {
    sc <- scale[v]
    sd <- sqrt(VarCorr(mod)[[1]]$var[v])
    b  <- coef(mod)[v]
    res <- integrate(fun_w, lower = -Inf, upper = Inf,
                     b = b, sd = sd, sc = sc)
    if (res$message == "OK") return(res$value)
    return(NA)
  })
  names(out) <- var
  out
}

sc2 <- sc
names(sc2) <- paste0(names(sc), "s")
tpi <- E_PI(fit, scale = sc2)


## === Simulation
simt <- function(seed) {
  set.seed(seed)
  ## --- Simulate from model
  sim <- simulate(fit, what = "response")
  data3CIA$ys <- sim
  dat2 <- data3CIA ## for OOS LL calculations
  dat2$ys <- simulate(fit, what = "response")
  ## --- No misspec
  time1 <- system.time({
    mm1 <-  try(
      CoxphME(ys ~ ages + mmrcs + fev1pps + (ages + mmrcs + fev1pps | cohort),
      data = data3CIA, log_first = TRUE, order = 4, support = c(1, 150), estinit = FALSE,
      control = optim_control(iter.max = 500, eval.max = 500, rel.tol = 1e-8))
    )
    pi1 <- E_PI(mm1, scale = sc2)
    ll1 <- logLik(mm1)
    oosll1 <- logLik(mm1, newdata = dat2)
    })
  if (inherits(mm1, "try-error") || mm1$opt$convergence != 0) {
    print(paste("Failed CoxphME:", seed))
    out <- list(pi1 = NA, ll1 = NA, oosll1 = NA)
  } else {
    out <- list(pi1 = pi1, ll1 = ll1, oosll1 = oosll1)
  }
  out$time1 <- time1
  ## -- Misspec
  time2 <- system.time({
    mm2 <-  try(
      ColrME(ys ~ ages + mmrcs + fev1pps + (ages + mmrcs + fev1pps | cohort),
      data = data3CIA, log_first = TRUE, order = 4, support = c(1, 150), estinit = FALSE,
      control = optim_control(iter.max = 500, eval.max = 500, rel.tol = 1e-8))
    )
    pi2 <- E_PI(mm2, scale = sc2)
    ll2 <- logLik(mm2)
    oosll2 <- logLik(mm2, newdata = dat2)
    })
  if (inherits(mm2, "try-error") || mm2$opt$convergence != 0) {
    print(paste("Failed ColrME:", seed))
    out$pi2 <- NA
    out$ll2 <- NA
    out$oosll2 <- NA
  } else {
    out$pi2 <- pi2
    out$ll2 <- ll2
    out$oosll2 <- oosll2
  }
  out$time2 <- time2
  if ((seed %% 10) == 0) print(paste("Done: ", seed))
  out$seed <- seed
  return(out)
}

## NOTE: control the number of draws by setting this vector
seeds <- 1001:2000
## sim_misspec <- lapply(seeds, function(ii) simt(ii))
sim_misspec <- parallel::mclapply(seeds, function(ii) simt(ii), mc.cores = 8)


## -- Non-convergent/problematic cases
which(sapply(sim_misspec, function(x) any(is.na(x$pi1))))
which(sapply(sim_misspec, function(x) any(is.na(x$pi2))))

## -- Model selection
mean(sapply(sim_misspec, function(x) x$ll1 > x$ll2), na.rm = TRUE)
mean(sapply(sim_misspec, function(x) x$oosll1 > x$oosll2), na.rm = TRUE)

## -- Figure S-10
pi_ph <- sapply(sim_misspec, function(x) {
  if (length(x$pi1) == 1 && is.na(x$pi1)) {
    out <- tpi
    out[] <- NA
    return(out)
  }
  x$pi1
})
pi_po <- sapply(sim_misspec, function(x) {
  if (length(x$pi2) == 1 && is.na(x$pi2)) {
    out <- tpi
    out[] <- NA
    return(out)
  }
  x$pi2
})

pdf("FigS-10.pdf", width = 7, height = 4.5)
par(mfrow = c(1, nrow(pi_ph)), las = 1, cex = 0.8,
    mar = c(7.4, 4.1, 4.1, 1))
c1 <- colorspace::qualitative_hcl(6, "Dark 3")[c(1, 3)]
c2 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.1)[c(1, 3)]
c3 <- colorspace::qualitative_hcl(6, "Dark 3", alpha = 0.7)[c(1, 3)]
## c1 <- gray(0); c2 <- gray(0, alpha = 0.1); c3 <- gray(0, alpha = 0.7)
titles <- c(expression(Age), expression(mMRC), expression(FEV[1]))
for (ii in 1:nrow(pi_ph)) {
  ylab <- if (ii == 1) "PI" else NULL
  boxplot(cbind(pi_ph[ii, ], pi_po[ii, ]), lwd = 2, xaxt = "n",
          col = c2, medcol=c3, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1,
          ylab = ylab, main = titles[ii])
  if (ii == 1)
    legend("bottomright", "True value", lwd = 2, col = 1, bty = "n", cex = 0.9)
  abline(h = tpi[ii], lwd = 2)
  grid(nx = 0, ny = NULL)
  axis(1, at = 1:2, labels = FALSE)
  aloc <- par("usr")[3] - diff(grconvertY(0:1, "inches", "user")) * 0.1
  text(1:2, aloc, srt = 50, adj = 1,
       labels = c("PH", "PO (misspec)"), xpd = TRUE)
}
dev.off()
