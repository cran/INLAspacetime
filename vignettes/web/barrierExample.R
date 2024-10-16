## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png",
  dev.args = list(type = "cairo-png"),
  fig.width = 10,
  fig.height = 7,
  dpi = 300,
  out.width = "99%",
  fig.align = "center"
)
knitr::opts_chunk$set(echo = TRUE)


## ----domain-------------------------------------------------------------------
rxy <- c(10, 6)
r <- mean(rxy)
bb <- rbind(
    c(0, 1)*rxy[1],
    c(0, 1)*rxy[2])
domain.xy <- cbind(
    c(0, 1, 1, 0, 0) * rxy[1],
    c(0, 0, 1, 1, 0) * rxy[2])


## ----bdomain------------------------------------------------------------------
barrier1 <- cbind(
  cos(seq(0, 2*pi, length=50)) * 0.7 + 0.50 * rxy[1], 
  sin(seq(0, 2*pi, length=50)) * 1.0 + 0.75 * rxy[2])
barrier2 <- cbind(
  c(0.50, 1.50, 1.50, 0.50, 0.50) * rxy[1],
  c(0.33, 0.33, 0.40, 0.40, 0.33) * rxy[2])


## ----sfpackage----------------------------------------------------------------
library(sf)
domain <- st_sfc(st_polygon(list(domain.xy)))
barriers <- st_sfc(st_multipolygon(
  list(st_polygon(list(barrier1)), 
       st_polygon(list(barrier2)))))


## ----domainvis----------------------------------------------------------------
library(ggplot2)
gg0 <- ggplot() + 
  xlab("") + ylab("") +
  theme_minimal() + coord_fixed()
gg0 + xlim(bb[1, ]) +
  geom_sf(data = domain, fill = "blue") +
  geom_sf(data = barriers, fill = "gray")


## ----packages-----------------------------------------------------------------
library(INLA)
library(INLAspacetime)
library(inlabru)
library(patchwork)


## ----mesh---------------------------------------------------------------------
mesh <- inla.mesh.2d(
    loc.domain = domain.xy, 
    max.edge = c(0.03, 0.1) * r,
    offset = c(0.1, 0.3) * r,
    cutoff = 0.01 * r)
mesh$n


## ----ibarrier-----------------------------------------------------------------
triBarrier <- unlist(fmesher::fm_contains(
  x = barriers, 
  y = mesh, 
  type = "centroid"))


## ----vmesh--------------------------------------------------------------------
triCenters.xy <- cbind(
  mesh$loc[mesh$graph$tv[,1], 1:2] +
  mesh$loc[mesh$graph$tv[,2], 1:2] +
  mesh$loc[mesh$graph$tv[,3], 1:2])/3
gg0 + 
    gg(mesh) +
    geom_point(aes(
        x = triCenters.xy[triBarrier, 1],
        y = triCenters.xy[triBarrier, 2])) 


## ----rfparams-----------------------------------------------------------------
sigma <- 1
(ranges <- r * c(0.5, 0.05)) 


## ----bfem---------------------------------------------------------------------
bfem <- mesh2fem.barrier(mesh, triBarrier)


## ----Q------------------------------------------------------------------------
Q <- inla.barrier.q(bfem, ranges = ranges, sigma = sigma)


## ----fcorr--------------------------------------------------------------------
localCorrel <- function(locs, mesh, Q) {
  nl <- nrow(locs)
  ii <- sapply(1:nl, function(i)
    which.min(rowSums(sweep(
      mesh$loc[, 1:ncol(locs)], 2, locs[i, ], "-")^2)))
  b <- matrix(0, nrow(Q), nl)
  for(i in 1:nl)
    b[ii[i], i] <- 1
  cc <- inla.qsolve(Q, b)
  s <- sqrt(diag(inla.qinv(Q)))
  for(i in 1:nl)
    cc[, i] <- cc[, i] / (s * s[ii[i]])
  return(drop(cc))
}


## ----lcorrels-----------------------------------------------------------------
locs <- cbind(c(0.4, 0.6, 0.7, 0.5) * rxy[1], 
              c(0.7, 0.6, 0.3, 0.5) * rxy[2])
mcorrels <- localCorrel(locs, mesh, Q)
dim(mcorrels) 


## ----pgrid--------------------------------------------------------------------
pgrid <- inla.mesh.projector(
    mesh,
    xlim = bb[1, ],
    ylim = bb[2, ],
    dims = round(500 * rxy/r)) 
gcorrels <- as.matrix(inla.mesh.project(
  pgrid, field = mcorrels
))


## ----grid.df------------------------------------------------------------------
grid.df <- data.frame(
  x = rep(pgrid$x, times = length(pgrid$y)), 
  y = rep(pgrid$y, each = length(pgrid$x)))


## ----gcorrels-----------------------------------------------------------------
ggcorrels <- do.call(
  rbind, 
  lapply(1:4, function(l) 
    data.frame(grid.df, 
               loc = paste(sprintf("%1.1f", locs[l, 1]), 
                           sprintf("%1.1f", locs[l, 2]), sep = ", "), 
               correlation = gcorrels[, l])))


## ----beval, echo = FALSE, ref.label = 'ggbarriers', eval = TRUE---------------

## ----ggcorrels----------------------------------------------------------------
gg0 + 
  geom_raster(
    data = ggcorrels[ggcorrels$correlation>0.1, ], 
    mapping = aes(x = x, y = y, fill = correlation)) + 
  facet_wrap(~ loc) + 
  add.b0 + gg.add ## look to the appendix for the code for this


## ----sample-------------------------------------------------------------------
u <- inla.qsample(1, Q, seed = 1)[,1]


## ----usimgrid-----------------------------------------------------------------
ugrid.sim <- inla.mesh.project(pgrid, field = u)


## ----uviz---------------------------------------------------------------------
grid.df$u <- as.vector(ugrid.sim)
gg0 + 
  geom_raster(
    data = grid.df, 
    aes(x = x, y = y, fill = u)) + 
  add.b0 + gg.add ## look to the appendix for the code for this


## ----locs0--------------------------------------------------------------------
n0 <- 500
set.seed(2)
xy0 <- cbind(
  runif(n0, bb[1, 1], bb[1, 2]),
  runif(n0, bb[2, 1], bb[2, 2]))


## ----ib-----------------------------------------------------------------------
ii <- which(sapply(st_intersects(
  x = st_as_sf(data.frame(x=xy0[, 1], y = xy0[, 2]), coords = c('x', 'y')),
  y = barriers), length)==0)

## ----locs---------------------------------------------------------------------
dataset <- data.frame(
  x = xy0[ii, 1], y = xy0[ii, 2])
(n <- nrow(dataset))


## ----dataset------------------------------------------------------------------
sigma.e <- 1
set.seed(3)
dataset$outcome <- 
    drop(inla.mesh.project(
        mesh,
        loc = cbind(dataset$x, dataset$y),
        field = u)) + 
  10 + rnorm(n, 0.0, sigma.e)


## ----bmodel-------------------------------------------------------------------
bmodel <- barrierModel.define(
    mesh = mesh, 
    barrier.triangles = triBarrier,
    prior.range = c(1, 0.01), ## P(range < 1) = 0.01
    prior.sigma = c(1, 0.01), ## P(sigma > 1) = 0.01
    range.fraction = 0.1)


## ----mformulae----------------------------------------------------------------
model <- outcome ~ Intercept(1) +
    field(cbind(x, y), model = bmodel)


## ----mfitcode, echo = FALSE, warning = FALSE, eval = FALSE--------------------
## result <- bru(
##     model, dataset, family = "gaussian")

## ----mfit, echo = FALSE-------------------------------------------------------
sresfl <- "result_barrier.R"
if(file.exists(sresfl)) {
  source(sresfl)
} else {
  result <- bru(
      model, dataset, family = "gaussian")
  snams <- c("summary.fixed", "summary.random", 
             "internal.marginals.hyperpar")
  result <- result[snams]
  dump("result", file = sresfl)
}


## ----fit, eval = FALSE--------------------------------------------------------
## result <- bru(
##     model, dataset, family = "gaussian")


## ----sfix---------------------------------------------------------------------
result$summary.fix


## ----pmargs-------------------------------------------------------------------
pmarginals <- 
    list(
        data.frame(
            param = "sigma.e",
            inla.tmarginal(
                function(x) exp(-x/2),
                result$internal.marginals.hyperpar[[1]])),
        data.frame(
            param = "range",
            inla.tmarginal(
                function(x) exp(x),
                result$internal.marginals.hyperpar[[2]])),
        data.frame(
            param = "sigma",
            inla.tmarginal(
                function(x) exp(x),
                result$internal.marginals.hyperpar[[3]]))
    )


## ----psummar------------------------------------------------------------------
rbind(true = c(sigma.e = sigma.e, range = ranges[1], sigma = 1),
      sapply(pmarginals, function(m)
             unlist(inla.zmarginal(m[, -1], TRUE))[1:2]))


## ----pmargsvis, fig.width = 7, fig.height = 4, out.width = "69%"--------------
ggplot(do.call(rbind, pmarginals)) + 
    geom_line(aes(x=x, y=y)) + 
    facet_wrap(~param, scales = "free") +
    theme_minimal()


## ----umean--------------------------------------------------------------------
grid.df$u.mean <- as.vector(
  inla.mesh.project(
    pgrid,
    result$summary.random$field$mean))
grid.df$u.sd <- as.vector(
  inla.mesh.project(
    pgrid,
    result$summary.random$field$sd))


## ----ibarrierpix--------------------------------------------------------------
gInBarrier <- which(sapply(st_intersects(
  x = st_as_sf(grid.df, coords = c('x', 'y')),
  y = barriers), length)==0)


## ----umeanviz-----------------------------------------------------------------
gg0 + 
  geom_raster(
    data = grid.df[!gInBarrier, ], 
    aes(x = x, y = y, fill = u.mean)) + 
  gg.add ## look to the appendix for the code for this


## ----usdviz-------------------------------------------------------------------
gg0 + 
  geom_raster(
    data = grid.df[!gInBarrier, ], 
    aes(x = x, y = y, fill = u.sd)) + 
  geom_point(data = dataset, aes(x = x, y = y)) +
  gg.add ## look to the appendix for the code for this


## ----old2new, echo = FALSE, results = "asis"----------------------------------
old2new <- data.frame(
  Old = c("inla.barrier.fem()", "inla.barrier.pcmatern()"),
  New = c("mesh2fem.barrier()", "barrierModel.define()"))
knitr::kable(old2new)


## ----ggbarriers, echo = TRUE, eval = FALSE------------------------------------
## gg.add <- list(
##   scale_fill_distiller(
##     type = "div",
##     palette = "RdBu",
##     na.value = "transparent")
## )
## add.b0 <- list(
##   geom_polygon(
##     mapping = aes(x, y),
##     data = data.frame(
##       x = pmin(barrier1[,1], bb[1, 2]),
##       y = pmin(barrier1[,2], bb[2, 2])
##     ), colour = "black",
##   fill = gray(0.9, 0.5)), ##"transparent"),
##   geom_polygon(
##     mapping = aes(x, y),
##     data = data.frame(
##       x = pmin(barrier2[,1], bb[1, 2]),
##       y = pmin(barrier2[,2], bb[2, 2])
##     ), colour = "black",
##     fill = gray(0.9, 0.5))##"transparent")
## )

