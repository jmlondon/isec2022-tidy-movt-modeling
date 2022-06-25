library(crawl)
library(ctmcmove)
library(raster)
library(lubridate)
library(tidyverse)
# library(nPacMaps) #devtools::install_github("jmlondon/nPacMaps")
library(viridis)
library(rasterVis)
library(sp)
library(mgcv)
library(foreach)

### Load my custom function for producing
source('ctmc2glm_alt.R')

### LOAD SOME DATA
# npac_poly = nPacMaps::npac() %>% as("Spatial")

npac_hab = readRDS("data/npac_habitat.rds") # See below if not made already

data(northernFurSeal)

pup_frame <- northernFurSeal %>% `coordinates<-`(~long+lat) %>%
  `proj4string<-`(CRS("+init=epsg:4326")) %>% spTransform(CRS(proj4string(npac_hab[[1]])))


### SET UP A CRAWL FIT
## Uses v 2.2.0

set.seed(321)
fit <- crwMLE(
  data=pup_frame,
  err.model=list(x=~loc_class-1),
  drift=TRUE,
  # initial.state = initial,
  Time.name="GMT",
  fixPar=c(log(250), log(500), log(1500), NA, NA, NA,4, NA, NA),
  theta=c(log(2000), log(2000), 5, 0, 0),
  constr=list(
    lower=c(rep(log(1500),2), rep(-Inf,3)), upper=rep(Inf,5)
  ), method='L-BFGS-B'
)

### Create the process imputer
post_simulator = crwSimulator(fit, predTime="10 mins", parIS=0)


### Group raster layers already in vector gradient form
grad_list = list(
  wind = stack(npac_hab$surface_wind_u, npac_hab$surface_wind_v),
  geo_curr = stack(npac_hab$geo_curr_u, npac_hab$geo_curr_v)
) %>% map(~'names<-'(.x, c("grad.x","grad.y")))

raster.list = list(
  stack.static=npac_hab$sst,
  grad.list=grad_list
)

### The SST layer has 'NA' over some land portions, so we need to let the model know about that.
zero.idx = c(
  which(is.na(raster::values(raster.list$stack.static$sst))),
  which(is.na(raster::values(raster.list$grad.list$geo_curr$grad.x)))
) %>% unique() %>% sort()

### A vector of daily times-- for effects prediction, etc
daily = seq(
  ceiling_date(pup_frame@data$GMT[1], "day"),
  tail(floor_date(pup_frame@data$GMT, "day"),1),
  24*60*60
  )

### Create plan for parallel computing-- see 'future' package
library(foreach)
library(doFuture)
library(doRNG)
registerDoFuture()
plan('multisession', workers=6)

### Execute imputation and estimation in parallel (using 25 imputed tracks)
MI_fit = tibble(rep=1:24) %>%
  mutate(
    PI_rep = foreach(i=1:n())%dorng%{
      examplerast = raster.list[[1]][[1]]
      # Generate imputed pathi
      pi_path = crwPostIS(post_simulator, fullPost=F) %>%
      {data.frame(.$alpha.sim, Time=.$Time, locType=.$locType)} %>%
        dplyr::filter(locType=="p") %>% dplyr::select(mu.x, mu.y, Time)
      # Get it in the right form for 'ctmcmove'
      path = list(t=pi_path$Time,xy=as.matrix(pi_path[,c("mu.x","mu.y")]))
      # Discretize the space
      ctmc = ctmcmove::path2ctmc(xy=path$xy,t=path$t, rast=examplerast, zero.idx = zero.idx)
      # Create model data for estimation
      glm_data = ctmc2glm_alt(ctmc,raster.list=raster.list, zero.idx = zero.idx)
      # Fit model
      form = "z ~ crw + s(t,by=sst) + s(t,by=wind_grad) + s(t,by=geo_curr_grad)"
      fit = gam(as.formula(form), select=T, family="poisson", offset=log(tau), data=glm_data)
      # Get the effects info
      Xp = smoothCon(s(t), data=glm_data)[[1]] %>%
        PredictMat(data=data.frame(t=as.numeric(daily)/3600))
      V_b = vcov(fit, unconditional = T)
      b = coef(fit)
      eff_nms = c("wind", "geo_curr", "sst")
      eff_list = foreach(k=1:3)%do%{
        idx = grep(eff_nms[k], names(b))
        mu_eff = Xp%*%b[idx]
        V_eff = Xp%*%V_b[idx,idx]%*%t(Xp)
        list(mu_eff, V_eff)
      } %>% `names<-`(eff_nms)
      list(fit=fit, eff_list=eff_list)
    }
    )
plan('sequential')

### Make some effects plots
library(mvtnorm)
library(coda)
library(cowplot)

### Draw a 'posterior' sample
post_sample = MI_fit %>%
  mutate(
    wind_posterior = map(PI_rep, ~data.frame(rmvnorm(5000, .x$eff_list$wind[[1]], .x$eff_list$wind[[2]]))),
    geo_curr_posterior = map(PI_rep, ~data.frame(rmvnorm(5000, .x$eff_list$geo_curr[[1]], .x$eff_list$geo_curr[[2]]))),
    sst_posterior = map(PI_rep, ~data.frame(rmvnorm(5000, .x$eff_list$sst[[1]], .x$eff_list$sst[[2]])))
  ) %>% select(rep, contains("posterior"))

### Plot wind effects
wind_post = post_sample %>% select(wind_posterior) %>% unnest(cols = c(wind_posterior))
wind_summ = data.frame(time=daily, est=colMeans(wind_post), hpd=HPDinterval(mcmc(wind_post),0.9))
wind_plt = ggplot(data=wind_summ)+geom_path(aes(x=time, y=est), col='blue') +
  geom_ribbon(aes(x=time, ymin=hpd.lower, ymax=hpd.upper), alpha=0.2, fill='blue') +
  geom_abline(slope=0, intercept = 0)+
  ylab("Wind effect") + xlab("Date")
ggsave(wind_plt, filename="wind_eff.png", width=8, height=6, units="in", dpi=200)

### Plot current effects
curr_post = post_sample %>% select(geo_curr_posterior) %>% unnest(cols = c(geo_curr_posterior))
curr_summ = data.frame(time=daily, est=colMeans(curr_post), hpd=HPDinterval(mcmc(curr_post),0.9))
curr_plt = ggplot(data=curr_summ)+geom_path(aes(x=time, y=est), col='blue') +
  geom_ribbon(aes(x=time, ymin=hpd.lower, ymax=hpd.upper), alpha=0.2, fill='blue') +
  geom_abline(slope=0, intercept = 0)+
  ylab("Geostraphic current effect") + xlab("Date")
ggsave(curr_plt, file="curr_eff.png", width=8, height=6, units="in", dpi=200)


### Plot SST effects
sst_post = post_sample %>% select(sst_posterior) %>% unnest(cols = c(sst_posterior))
sst_summ = data.frame(time=daily, est=colMeans(sst_post), hpd=HPDinterval(mcmc(sst_post),0.9))
sst_plt = ggplot(data=sst_summ)+geom_path(aes(x=time, y=est), col='blue') +
  geom_ribbon(aes(x=time, ymin=hpd.lower, ymax=hpd.upper), alpha=0.2, fill='blue') +
  geom_abline(slope=0, intercept = 0)+
  ylab("SST effect") + xlab("Date")
ggsave(sst_plt, filename="sst_eff.png", width=8, height=6, units="in", dpi=200)
