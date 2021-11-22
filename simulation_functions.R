# Version 2 of the simulations
# (after reviewer comments from Epidemiology)


##' @name measure_intervals
##' @description Creates a set of observation times based on a simulation structure
##' 
measure_intervals <- function(n,
                              interval_length=1,
                              interval_gap=0,
                              J=1,
                              random_initial_only=FALSE,
                              start=0,
                              start_interval_length=NULL) {
  if (J==1){
    times <- runif(n, 0, interval_length)  
  } else if (J>1) {
    if(random_initial_only) {
      t1 <- runif(n, 0, interval_length)  
      times <- rep(t1, each=J) + interval_length*rep(0:(J-1), times=n)
    } else {
      times <- runif(n*J, 0, interval_length)  
      if (length(interval_gap)==1) {
        times <- times + (interval_length+interval_gap)*rep(0:(J-1), times=n)  
      } else if (length(interval_gap)==J){
        times <- times + rep(interval_gap, times=n)        
      } else {
        stop()
      }
    }
  } else {
    stop("J must be >= 1.")
  }
  if (!is.null(start_interval_length)){
    start <- start + rep(runif(n, min=0, max=start_interval_length), each=J)
  }
  times <- times + start
  times
}



##' @name create_expsim_dataset
##' @description Creates a simulated HAP dataset
##' @param ngroups Number of groups in the study
##' @param study_type Indicator of study type. Currently only "parallel" and "crossover" supported.
##' @param J Number of observations per group
##' @param J1 Number of observations in the first condition for the first group of a crossover study
##' @param groupmeans The group mean HAP values (often on log scale).
##' @param sigE Observation standard error
##' @param sigA Standard deviation of household random effects
##' @param timefn Function for temporal variation. Was previously assumed to have mean zero; that is no longer required.
##' @param times Times of each measurement. If equal to "sample", then `measure_intervals()` is called.
##' @param truth_timerange Time range of values to compute the (true) LTA over
##' @details This code is based on functions in the `bercs` package.
##' @import bercs
create_expsim_dataset <- function(ngroups=2,
                                  study_type="parallel",
                                  n=200, # per group
                                  J=2,
                                  J1=J,
                                  J2=J,
                                  firstgroup=1,
                                  groupmeans=c(3, 4),
                                  sigE=1,
                                  sigA=1,
                                  verbose=FALSE,
                                  timefn=function(w) {0}, ## no longer assumed to have mean zero.
                                  times=NULL,
                                  truth_timerange=range(times), ## added this now.
                                  return_structure=FALSE,
                                  ...) {
  

  if (study_type=="parallel"){
    
  
  es <- create_exposure_simulation_skeleton_parallel(ngroups=ngroups,
                                                     nclusters=1,
                                                     nunits=n,
                                                     nobs=J,
                                                     verbose=verbose)
  } else if (study_type=="crossover"){

        es <- create_exposure_simulation_skeleton_crossover(ngroups=ngroups,
                                                        nclusters=1,
                                                        nunits=n,
                                                        nobs1=rep(c(J1, J2), each=n/2),
                                                        nobs2=rep(c(J2, J1), each=n/2),
                                                        firstgroup=firstgroup,
                                                        verbose=verbose)
    J <- J1 + J2
  } else {
    stop("'study_type' not supported.")
  }
  
  # Set group means
  es <- expsim_update_parameter(es,
                                level="group",
                                type="mean",
                                value=groupmeans)
  # Set unit random effect standard deviation
  es <- expsim_update_parameter(es,
                                level="unit",
                                type="sd",
                                value=sigA)
  # Set observation standard deviation
  es <- expsim_update_parameter(es,
                                level="observation",
                                type="sd",
                                value=sigE)
  # Add time function
  if (!is.null(timefn)){
    es <- expsim_update_parameter(es, level="time", type="mean", value=timefn)
  }
  # Add times
  if (!is.null(times)){
    if (times[1]=="sample") {
      if (study_type=="parallel"){
      times <- measure_intervals(n=n*2, J=J,...)
      } else if (study_type=="crossover"){
        times <- measure_intervals(n=n, J=J1 + J2,...)
      } 
    } 
    es <- sim_update_times(es, time=times)
  }
  
  
  # Sample unit random effects
  es <- expsim_update_parameter(es,
                                level="unit",
                                type="re")
  
  # Create start/stop times for each stove type (by person)

  
  # Compute (conditional) observation means and sample observations
  expsim_demo <- expsim_sample_observations(es) 
  dat <- with(expsim_demo$standata,
              data.frame(group=group_of_obs,
                         unit=unit_of_obs,
                         obs_num=rep(1:J, times=n),
                         cluster=cluster_of_obs,
                         x=w,
                         time=time,
                         mu_wtime=expsim_demo$structure$meanW))
  
  
  dat <- dat %>%
    group_by(unitstove=interaction(group, unit)) %>%
  mutate(stove_stop=max(time)) %>%
    ungroup()

  if (study_type=="crossover"){
    # Assumes all start with group 1
    dat <- dat %>%
      group_by(unit) %>%
      mutate(stove_start=case_when(group==1~min(time),
                                   group==2~max(time[group==1]))) %>%
      ungroup()
  } else {
    dat <- dat %>%
      group_by(unitstove) %>%
      mutate(stove_start=min(time)) %>%
      ungroup()
  }
  

myintegrate <- function(fn, start, stop){
  out <- numeric(length(start))
  for (i in 1:length(out)){
   out[i] <- integrate(fn, start[i], stop[i])$value/(stop[i] - start[i])
  }
  out
}
  dat <- dat %>%
    mutate(mu_notime=mu_wtime - timefn(time),
           mu_timerange=mu_wtime - timefn(time)  + integrate(timefn, truth_timerange[1], truth_timerange[2])$value/(truth_timerange[2] - truth_timerange[1]))
  
  
  
    dat$mu_obsrange <- dat$mu_wtime - timefn(dat$time) + myintegrate(timefn, dat$stove_start, dat$stove_stop)
                         
  if (return_structure) {
    dat <- list(dat=dat,
                structure=expsim_demo)
  }
  dat
}

predict_dropFEs <- function(model,
                            inds_keep=0, ...){
  if (inherits(model, "merMod")) {
    if (length(inds_keep)==1 && inds_keep==0){
      out <- predict(model, ...) 
    } else {
      out <-  predict(model, ...) - model.matrix(model)[, -inds_keep, drop=FALSE]%*% fixef(model)[-inds_keep] 
    }
  } else {
    if (length(inds_keep)==1 && inds_keep==0){
      out <- predict(model, ...) 
    } else {
      out <-  predict(model, ...) - model.matrix(model)[, -inds_keep, drop=FALSE]%*% coef(model)[-inds_keep] 
    }
  }
  out
} 

# Gets prediction of the LTA for a given time range.
predict_avgtime <- function(model, dat, group_var,timerange_min, timerange_max,seq_length=200,...){
  if (length(timerange_min)==1 && length(timerange_max)==1){
    time <- seq(timerange_min, timerange_max, length=seq_length)
    dat <- dat %>%
      select(-c(time)) %>%
      distinct({{group_var}}, .keep_all = TRUE) %>%
      expand_grid(time)
  } else {
    stop("not supported yet")
  }
  if (!is.null(dat$stove_start)){
    dat <- dat %>%
      filter(time >=stove_start)
  }
  if (!is.null(dat$stove_stop)){
    dat <- dat %>%
      filter(time <=stove_stop)
  }
  out <- predict_avg(model=model,
              dat=dat,
              group_var={{group_var}},
              ...)
  out
}

predict_avg <- function(model,
                            dat,
                            group_var,
                            ...){
  out <- predict(model,newdata=dat,...) 
  out <- cbind(dat, pred=out)
  if (!missing(group_var)) {
    out <- out  %>%
      group_by({{group_var}}) %>%
      summarize(pred=mean(pred))
  }
  out
} 

extract_mod_LTA <- function(mod, dat, group_var, timerange_min, timerange_max){
  predout <- predict_avgtime(mod, dat=dat, group_var={{group_var}}, timerange_min=timerange_min, timerange_max=timerange_max)
  out <- list(lta=predout$pred,
              unit_id=predout[[1]])
  out
}


extract_mod_values <- function(mod, dat){
 if (inherits(mod, "merMod")){
    vcovmod <- as.data.frame(VarCorr(mod))
    out <- list(fixed=fixef(mod),
                unit_re=ranef(mod)[["unit"]]$`(Intercept)`,
                unit_sd=vcovmod$sdcor[vcovmod$grp=="unit"],
                resid_sd=vcovmod$sdcor[vcovmod$grp=="Residual"])
  } else if (inherits(mod, "lm")) {
    out <- list(fixed=coef(mod),
                unit_re=NA,
                unit_sd=NA,
                resid_sd=summary(mod)$sigma)
  }
  out
}





# LET THIS ACCOUNT FOR TIME!
extract_mod_values_v0 <- function(mod, inds_keep=0){
  if (inherits(mod, "merMod")){
    vcovmod <- as.data.frame(VarCorr(mod))
    out <- list(fixed=fixef(mod),
                unit_re=ranef(mod)[["unit"]]$`(Intercept)`,
                # unit_est_lta=predict(mod, type="response"),
                unit_est_lta=predict_dropFEs(mod, type="response", inds_keep=inds_keep),
                unit_sd=vcovmod$sdcor[vcovmod$grp=="unit"],
                resid_sd=vcovmod$sdcor[vcovmod$grp=="Residual"])
  } else if (inherits(mod, "lm")) {
    out <- list(fixed=coef(mod),
                unit_re=NA,
                # unit_est_lta=fitted(mod),
                unit_est_lta=predict_dropFEs(mod, type="response", inds_keep=inds_keep),
                unit_sd=NA,
                resid_sd=summary(mod)$sigma)
  }
  out
}


##' @name get_lta_metrics
##' @description Calculates MSE/bias metrics for a model
##' @param mod Model object to predict from
##' @param true_obsrange Calculate metrics for true LTA from mu_obsrange?
##' @param true_timerange Calculate metrics for true LTA from mu_timerange?
get_lta_metrics <- function(mod,
                            dat,
                            group_var=unitid,
                            timerange=c(0, 1),
                            true_obsrange=TRUE,
                            true_timerange=TRUE,
                            pred_obsrange=TRUE,
                            pred_timerange=TRUE){
  MSE <- numeric(4)
  names(MSE) <- c("TobsPobs", "TobsPtime", "TtimePobs", "TtimePtime")
  bias <- MSE
  dat_nodup <- dat %>%
    distinct({{group_var}}, .keep_all=TRUE)
  # Calculate LTA using timerange
  if (pred_timerange){
    lta_timerange <- extract_mod_LTA(mod,
                                     dat=dat,
                                     group_var={{group_var}},
                                     timerange_min=timerange[1],
                                     timerange_max=timerange[2])
  
    if (true_obsrange){
      MSE["TobsPtime"] <- get_lta_mse(dat_pred=lta_timerange,
                                      dat_true=dat_nodup,
                                      var_pred="lta",
                                      var_true="mu_obsrange")
      bias["TobsPtime"] <- get_lta_bias(dat_pred=lta_timerange,
                                        dat_true=dat_nodup,
                                        var_pred="lta",
                                      var_true="mu_obsrange")
      
      
    }
    if (true_timerange){
      MSE["TtimePtime"] <- get_lta_mse(dat_pred=lta_timerange,
                                       dat_true=dat_nodup,
                                       var_pred="lta",
                                      var_true="mu_timerange")
      bias["TtimePtime"] <- get_lta_bias(dat_pred=lta_timerange,
                                         dat_true=dat_nodup,
                                         var_pred="lta",
                                       var_true="mu_timerange")
    }
  }
  
  
  # Calculate LTA using obsrange
  if (pred_obsrange){
    lta_obsrange <- extract_mod_LTA(mod,
                                     dat=dat,
                                     group_var={{group_var}},
                                    timerange_min=min(dat$time),
                                    timerange_max=max(dat$time))
    
    if (true_obsrange){
      MSE["TobsPobs"] <- get_lta_mse(dat_pred=lta_obsrange,
                                     dat_true=dat_nodup,
                                     var_pred="lta",
                                      var_true="mu_obsrange")
      bias["TobsPobs"] <- get_lta_bias(dat_pred=lta_obsrange,
                                       dat_true=dat_nodup,
                                       var_pred="lta",
                                         var_true="mu_obsrange")
    }
    if (true_timerange){
      MSE["TtimePobs"] <- get_lta_mse(dat_pred=lta_obsrange,
                                      dat_true=dat_nodup,
                                      var_pred="lta",
                                       var_true="mu_timerange")
      bias["TtimePobs"] <- get_lta_bias(dat_pred=lta_obsrange,
                                        dat_true=dat_nodup,
                                        var_pred="lta",
                                       var_true="mu_timerange")
    }
  }

  
  out <- list(MSE=MSE,
              bias=bias)
  out
}
  
  

##' @name get_lta_mse
##' @description Calculate MSE of a predicted LTA
get_lta_mse <- function(dat_pred,
                        dat_true,
                        var_pred="lta",
                        var_true="mu"){
  mse <- get_mse(truth=dat_true[[var_true]],
                 pred=dat_pred[[var_pred]])
  mse
}

get_mse <- function(truth,
                    pred){
  mse <- mean((pred- truth )^2)
  mse
}

get_lta_mse0 <- function(res,
                        dat,
                        muvar="mu"){
  mse <- mean((res$unit_est_lta[!duplicated(dat$unit)] - dat[[muvar]][!duplicated(dat$unit)])^2)
  mse
}

##' @rdname get_lta_mse
##' @description Calculate MSE of a predicted LTA
get_lta_bias <- function(dat_pred,
                        dat_true,
                        var_pred="lta",
                        var_true="mu"){
  mse <- get_bias(truth=dat_true[[var_true]],
                 pred=dat_pred[[var_pred]])
  mse
}
get_bias <- function(truth,
                    pred){
  bias <- mean((pred- truth ))
  bias
}


get_lta_bias_v0 <- function(res,
                         dat){
  bias <- mean((res$unit_est_lta[!duplicated(dat$unit)] - dat$mu[!duplicated(dat$unit)]))
  bias
}


run_one_sim <- function(ngroups=2,
                        study_type="parallel",
                        n=200, # per group
                        J=2,
                        J1=J,
                        J2=J,
                        firstgroup=1,
                        groupmeans=c(3, 4),
                        sigE=1,
                        sigA=1,
                        verbose=FALSE,
                        timefn=NULL, ## **ASSUMED TO HAVE MEAN ZERO**
                        times=NULL,
                        time_spline_df=NULL,
                        time_spline_df2=NULL,
                        truth_timerange_min=min(times),
                        truth_timerange_max=max(times),
                        ...){
  dat <- create_expsim_dataset(ngroups=ngroups,
                               study_type=study_type,
                               groupmeans=groupmeans,
                               n=n,
                               J=J,
                               J1=J1,
                               J2=J2,
                               firstgroup=firstgroup,
                               sigE=sigE,
                               sigA=sigA,
                               verbose=verbose,
                               timefn=timefn,
                               times=times,
                               truth_timerange=c(truth_timerange_min,
                                                 truth_timerange_max),
                               ...)

    # HH Average Model
  mod0 <- lm(x~0 + factor(unit), data=dat)
  
  res0 <- extract_mod_values(mod0,dat=dat)
  metrics0 <- get_lta_metrics(mod=mod0,
                              dat=dat,
                              group_var=unitstove,
                              timerange=c(truth_timerange_min,truth_timerange_max))

  # LMM: Group + HH-RE
  mod1 <- lmer(x~factor(group) + (1|unit), data=dat)
  res_wgroup <- extract_mod_values(mod1, dat=dat)
  metrics_wgroup <- get_lta_metrics(mod=mod1,
                              dat=dat,
                              group_var=unitstove,
                              timerange=c(truth_timerange_min,truth_timerange_max))
  
  # LMM: HH-RE only
  mod_nogroup <- lmer(x~(1|unit), data=dat)
  res_nogroup <- extract_mod_values(mod_nogroup, dat=dat)
  metrics_nogroup <- get_lta_metrics(mod=mod_nogroup,
                                    dat=dat,
                                    group_var=unitstove,
                                    timerange=c(truth_timerange_min,truth_timerange_max))
  
  # LMM: Factor Time + HH-RE
  dat$obs_num_contr_sum <- factor(dat$obs_num)
  contrasts(dat$obs_num_contr_sum) <- contr.sum(n=nlevels(dat$obs_num_contr_sum))
  mod_factime <- lmer(x~  factor(group) + obs_num_contr_sum   +  (1|unit), data=dat)
  res_factime <- extract_mod_values(mod_factime, dat=dat)
  metrics_factime <- get_lta_metrics(mod=mod_factime,
                                     dat=dat,
                                     group_var=unitstove,
                                     timerange=c(truth_timerange_min,truth_timerange_max))

  if (!is.null(times)){
    # # LMM: Linear Time
    # mod_lintime <- lmer(x~factor(group) + scale(time)  + (1|unit), data=dat)
    # res_lintime <- extract_mod_values(mod_lintime, inds_keep=1:(ngroups))
    # mse_lintime <- get_lta_mse(res_lintime, dat)
    # if (length(unique(dat$time))<10){

    # } else {
    #   mse_factime <- NA
    #   bias_factime <- NA
    #   res_factime <- list(unit_sd=NA,
    #                       resid_sd=NA)
    # }
    # LMM: Spline Time + HH-RE
    if (!is.null(time_spline_df)){
      # print(time_spline_df)
      if (time_spline_df>length(unique(dat$time))) {
        time_spline_df <- length(unique(dat$time))
      }
      
    
      mod_splinetime <- lmer(x~factor(group) + ns(time, df=time_spline_df)  + (1|unit), data=dat)

      res_splinetime <- extract_mod_values(mod_splinetime, dat=dat)
      metrics_splinetime <- get_lta_metrics(mod=mod_splinetime,
                                            dat=dat,
                                            group_var=unitstove,
                                            timerange=c(truth_timerange_min,truth_timerange_max))
      
      

      if (time_spline_df2>length(unique(dat$time))) {
        time_spline_df2 <- length(unique(dat$time))
      }
      mod_splinetime2 <- lmer(x~factor(group) + scale(ns(time, df=time_spline_df2))  + (1|unit), data=dat)
      res_splinetime2 <- extract_mod_values(mod_splinetime2, dat=dat)
      metrics_splinetime2 <- get_lta_metrics(mod=mod_splinetime2,
                                            dat=dat,
                                            group_var=unitstove,
                                            timerange=c(truth_timerange_min,truth_timerange_max))
    }
  } 
  # if (is.null(times)){
  #   mse_factime <- NA
  #   bias_factime <- NA
  #   res_factime <- list(unit_sd=NA,
  #                       resid_sd=NA)
  #   # mse_lintime <- NA
  #   # res_lintime <- list(unit_sd=NA,
  #   # resid_sd=NA)
  # }
  if (is.null(time_spline_df)){
    mse_splinetime_fullscale <- NA
    bias_splinetime_fullscale <- NA
    res_splinetime_fullscale <- list(unit_sd=NA,
                           resid_sd=NA)
    mse_splinetime <- NA
    bias_splinetime <- NA
    res_splinetime <- list(unit_sd=NA,
                                    resid_sd=NA)
  }
  if (is.null(time_spline_df2)){
    mse_splinetime2 <- NA
    bias_splinetime2 <- NA
    res_splinetime2 <- list(unit_sd=NA,
                           resid_sd=NA)
  }
  
  
  out <- list(MSE=c(hh_avg=metrics0$MSE,
                    wgroup=metrics_wgroup$MSE,
                    nogroup=metrics_nogroup$MSE,
                    factime=metrics_factime$MSE,
                    splinetime=metrics_splinetime$MSE,
                    splinetime2=metrics_splinetime2$MSE),
              bias=c(hh_avg=metrics0$bias,
                     wgroup=metrics_wgroup$bias,
                     nogroup=metrics_nogroup$bias,
                     factime=metrics_factime$bias,
                     splinetime=metrics_splinetime$bias,
                     splinetime2=metrics_splinetime2$bias),
              sigA=c(hh_avg=res0$unit_sd,
                     wgroup=res_wgroup$unit_sd,
                     nogroup=res_nogroup$unit_sd,
                     # lintime=res_lintime$unit_sd,
                     factime=res_factime$unit_sd,
                     splinetime=res_splinetime$unit_sd,
                     splinetime2=res_splinetime2$unit_sd),
              sigE=c(hh_avg=res0$resid_sd,
                     wgroup=res_wgroup$resid_sd,
                     nogroup=res_nogroup$resid_sd,
                     # lintime=res_lintime$resid_sd,
                     factime=res_factime$resid_sd,
                     splinetime=res_splinetime$resid_sd,
                     splinetime2=res_splinetime2$resid_sd))
  unlist(out)


}

run_many <- function(B=1,
                     ngroups=2,
                     study_type="parallel",
                     n=200, # per group
                     J=2,
                     J1=J,
                     J2=J,
                     firstgroup=1,
                     groupmeans=c(3, 4),
                     sigE=1,
                     sigA=1,
                     verbose=FALSE,
                     timefn=NULL, ## **ASSUMED TO HAVE MEAN ZERO**
                     times=NULL,
                     time_spline_df=NULL,
                     time_spline_df2=NULL,
                     truth_timerange_min=min(times),
                     truth_timerange_max=max(times),
                     cores=1,
                     ...
) {
  res <- mclapply(1:B, function(w) run_one_sim(ngroups=ngroups,
                                               study_type=study_type,
                                               groupmeans=groupmeans,
                                               n=n,
                                               J=J,
                                               J1=J1,
                                               J2=J2,
                                               firstgroup=firstgroup,
                                               sigE=sigE,
                                               sigA=sigA,
                                               verbose=verbose,
                                               timefn=timefn,
                                               times=times,
                                               time_spline_df=time_spline_df,
                                               time_spline_df2=time_spline_df2,
                                               truth_timerange_min=truth_timerange_min,
                                               truth_timerange_max=truth_timerange_max,
                                               ...
  ),
  mc.cores=cores)
  res <- simplify2array(res)
  as.data.frame(t(rowMeans(res)))
}




