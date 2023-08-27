## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  cache = !TRUE, 
  cache.path = "cache/",
  collapse = TRUE,
  comment = "#>",
  dev = "png",
  dev.args = list(type = "cairo-png"),
  fig.width = 7,
  fig.height = 5
)


## ----packages, results = "hide", cache = FALSE--------------------------------
library(ggplot2)
library(patchwork)
library(INLA)
library(INLAspacetime)
library(inlabru)


## ----inlasetup----------------------------------------------------------------
inla.setOption(
  inla.mode = "compact",
  num.threads = "5:-1",
  smtp = "pardiso",
  pardiso.license = "~/.pardiso.lic")


## ----ctrc---------------------------------------------------------------------
ctrc <- list(
  waic = TRUE,
  dic = TRUE,
  cpo = TRUE)


## ----filenames----------------------------------------------------------------
u0 <- paste0(
    "http://inla.r-inla-download.org/",
    "r-inla.org/case-studies/Cameletti2012/")
coofl <- "coordinates.csv"
datafl <- "Piemonte_data_byday.csv"
bordersfl <- "Piemonte_borders.csv"


## ----downloadporders----------------------------------------------------------
### get the domain borders
if(!file.exists(bordersfl))
    download.file(paste0(u0, bordersfl), bordersfl)
dim(pborders <- read.csv(bordersfl))


## ----coords-------------------------------------------------------------------
### get the coordinates
if(!file.exists(coofl))
    download.file(paste0(u0, coofl), coofl)
dim(locs <- read.csv(coofl))


## ----downloaddata-------------------------------------------------------------
### get the dataset
if(!file.exists(datafl))
    download.file(paste0(u0, datafl), datafl)
dim(pdata <- read.csv(datafl))


## ----dinspect-----------------------------------------------------------------
head(pdata)


## ----time---------------------------------------------------------------------
range(pdata$Date <- as.Date(pdata$Date, "%d/%m/%y"))
pdata$time <- as.integer(difftime(
    pdata$Date, min(pdata$Date), units = "days")) + 1


## ----covars-------------------------------------------------------------------
### prepare the covariates
xnames <- c("A", "WS", "TEMP", "HMIX", "PREC", "EMI")
xmean <- colMeans(pdata[, xnames])
xsd <- sapply(pdata[xnames], sd)

### prepare the data (st loc, scale covariates and log PM10)
dataf <- data.frame(pdata[c("UTMX", "UTMY", "time")],
                    scale(pdata[xnames], xmean, xsd),
                    y = log(pdata$PM10))
str(dataf)


## ----tmesh--------------------------------------------------------------------
nt <- max(pdata$time)
h <- 1
tmesh <- inla.mesh.1d(
  loc = seq(1, nt + h/2, h), 
  degree = 1)
tmesh$n


## ----smesh--------------------------------------------------------------------
smesh <- inla.mesh.2d(
    cbind(locs[,2], locs[,3]),
    loc.domain = pborders,
    max.edge = c(50, 300),
    offset = c(10, 140),
    cutoff = 5,
    min.angle = c(26, 21))
smesh$n


## ----smeshvis-----------------------------------------------------------------
par(mfrow = c(1,1), mar = c(0,0,1,0))
plot(smesh, asp = 1)
lines(pborders, lwd = 2, col = "green4")
points(locs[, 2:3], pch = 19, col = "blue")


## ----lprec--------------------------------------------------------------------
lkprec <- list(
    prec = list(prior = "pcprec", param = c(1, 0.05)))


## ----lhood--------------------------------------------------------------------
lhood <- like(
  formula = y ~ .,
  family = "gaussian",
  control.family = list(
    hyper = lkprec),
  data = dataf)


## ----mcomp--------------------------------------------------------------------
M <- ~ -1 + Intercept(1) + A + WS + TEMP + HMIX + PREC + EMI +
    field(list(space = cbind(UTMX, UTMY), 
               time = time),
          model = stmodel)


## ----stmodeldef, echo = FALSE, eval = FALSE-----------------------------------
## stmodel <- stModel.define(
##     smesh, tmesh, model,
##     control.priors = list(
##         prs = c(100, 0.05),
##         prt = c(1, 0.05),
##         psigma = c(2, 0.05)),
##     constr = TRUE)


## ----m102---------------------------------------------------------------------
model <- "102"
stmodel <- stModel.define(
    smesh, tmesh, model,
    control.priors = list(
        prs = c(150, 0.05),
        prt = c(10, 0.05),
        psigma = c(5, 0.05)),
    constr = TRUE)


## ----ini----------------------------------------------------------------------
theta.ini <- c(4, 7, 9, 1)


## ----fitcode, echo = FALSE, eval = FALSE--------------------------------------
##   bru(M,
##       lhood,
##       options = list(
##           control.mode = list(theta = theta.ini, restart = TRUE),
##           control.compute = ctrc))


## ----fit102-------------------------------------------------------------------
fit102 <- 
    bru(M,
        lhood,
        options = list(
            control.mode = list(theta = theta.ini, restart = TRUE),
            verbose = TRUE,
            control.compute = ctrc))


## ----cpu1---------------------------------------------------------------------
fit102$cpu


## ----sfixef-------------------------------------------------------------------
fit102$summary.fixed[, c(1, 2, 3, 5)]


## ----thyper-------------------------------------------------------------------
post.h <- list(
  sigma_e = inla.tmarginal(function(x) exp(-x/2), 
                           fit102$internal.marginals.hyperpar[[1]]),
  range_s = inla.tmarginal(function(x) exp(x), 
                           fit102$internal.marginals.hyperpar[[2]]),
  range_t = inla.tmarginal(function(x) exp(x), 
                           fit102$internal.marginals.hyperpar[[3]]),
  sigma_u = inla.tmarginal(function(x) exp(x), 
                           fit102$internal.marginals.hyperpar[[4]])
)


## ----shyper-------------------------------------------------------------------
shyper <- t(sapply(post.h, function(m) 
  unlist(inla.zmarginal(m, silent = TRUE))))
shyper[, c(1, 2, 3, 7)]


## ----compare1-----------------------------------------------------------------
c(shyper[c(1, 4, 2), 1], 
  a = exp(-h * sqrt(8 * 0.5) / shyper[3, 1]))


## ----m121---------------------------------------------------------------------
model <- "121"
stmodel <- stModel.define(
    smesh, tmesh, model,
    control.priors = list(
        prs = c(200, 0.05),
        prt = c(5, 0.05),
        psigma = c(2, 0.05)),
    constr = TRUE)


## ----fit121-------------------------------------------------------------------
fit121 <- 
    bru(M,
        lhood,
        options = list(
            control.mode = list(theta = theta.ini, restart = TRUE),
            verbose = TRUE,
            control.compute = ctrc))


## ----results------------------------------------------------------------------
results <- list("u102" = fit102, "u121" = fit121)


## ----cpu----------------------------------------------------------------------
sapply(results, function(r) r$cpu)


## ----nfn----------------------------------------------------------------------
sapply(results, function(r) r$misc$nfunc)


## ----pmode--------------------------------------------------------------------
sapply(results, function(r) r$mode$theta)

par(mfrow = c(1, 1), mar = c(3, 3, 0, 0.0), mgp = c(2, 1, 0))
uu.hist <- lapply(results, function(r)
    hist(r$summary.random$field$mean,
         -60:60/20, plot = FALSE))
ylm <- range(uu.hist[[1]]$count, uu.hist[[2]]$count)
plot(uu.hist[[1]], ylim = ylm,
     col = rgb(1, 0.1, 0.1, 1.0), border = FALSE, 
     xlab = "u", main = "")
plot(uu.hist[[2]], add = TRUE, col = rgb(0.1, 0.1, 1, 0.5), border = FALSE)
legend("topleft", c("separable", "non-separable"), 
       fill = rgb(c(1,0.1), 0.1, c(0.1, 1), c(1, 0.5)), 
       border = 'transparent', bty = "n")

## ----h2pmd--------------------------------------------------------------------
posts.h2 <- lapply(1:2, function(m) vector("list", 4L))
for(m in 1:2) {
    posts.h2[[m]]$sigma_e =  
      data.frame(
        parameter = "sigma_e", 
        inla.tmarginal(
          function(x) exp(-x/2), 
          results[[m]]$internal.marginals.hyperpar[[1]]))
    for(p in 2:4) {
      posts.h2[[m]][[p]] <-   
      data.frame(
        parameter = c(NA, "r_s", "r_t", "sigma")[p], 
        inla.tmarginal(
          function(x) exp(x), 
          results[[m]]$internal.marginals.hyperpar[[p]])
      )
    }
}


## ----hpmds2-------------------------------------------------------------------
posts.df <- rbind(
  data.frame(model = "102", do.call(rbind, posts.h2[[1]])),
  data.frame(model = "121", do.call(rbind, posts.h2[[2]]))
)

ggplot(posts.df) +
  geom_line(aes(x = x, y = y, group = model, color = model)) +
  ylab("Density") + xlab("") + 
  facet_wrap(~parameter, scales = "free")


## ----sfits--------------------------------------------------------------------
t(sapply(results, function(r) {
  c(DIC = mean(r$dic$local.dic, na.rm = TRUE),
    WAIC = mean(r$waic$local.waic, na.rm = TRUE),
    LPO = -mean(r$po$po, na.rm = TRUE), 
    LCPO = -mean(r$cpo$cpo, na.rm = TRUE))
}))


## ----g5cv---------------------------------------------------------------------
g5cv <- lapply(
  results, inla.group.cv, num.level.sets = 5, 
  strategy = "posterior", size.max = 50)


## ----group1m102---------------------------------------------------------------
ii <- 1000
g5cv$u102$group[[ii]]


## ----group1m121---------------------------------------------------------------
g5cv$u121$group[[ii]]


## ----obs1---------------------------------------------------------------------
dataf[g5cv$u102$group[[ii]]$idx, ]


## ----nlcv---------------------------------------------------------------------
sapply(g5cv, function(r) -mean(log(r$cv), na.rm = TRUE))

table(compare.cv <- g5cv[[1]]$cv < g5cv[[2]]$cv)

ii2 <- c(which(compare.cv)[ii],
         which(!compare.cv)[ii])
ii2

hh <- tapply(dataf$y, compare.cv, hist, breaks = 50, plot = FALSE)

par(mfrow = c(2, 2), mar = c(3, 3, 0, 0), mgp = c(2, 1, 0))
plot(dataf$y ~ factor(compare.cv))
ggplot(dataf) +
    geom_boxplot(aes(x = time, y = y, group = compare.cv))

par(mfrow = c(1, 2), mar = c(3, 3, 0, 0), mgp = c(2, 1, 0))
for (i in ii2) {
    x0 <- seq(dataf$y[i] - 1, dataf$y[i] + 1, 0.01)
    plot(x0, dnorm(x0, g5cv[[1]]$mean[i], g5cv[[1]]$sd[i]), type = "l", lwd = 2)
    lines(x0, dnorm(x0, g5cv[[2]]$mean[i], g5cv[[2]]$sd[i]), col = 4, lwd = 2)
    abline(v = dataf$y[i], col = 2, lwd = 3)
    c(g5cv[[1]]$cv[i], g5cv[[2]]$cv[i])
}
