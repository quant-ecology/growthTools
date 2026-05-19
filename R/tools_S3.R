# ============================================================
# Growth Rate Modeling Framework (S3 OOP Version)
# ============================================================


#' Fitting methods for growth models
#' 
#' A suite of growth models is available which can be fit to time series of 
#' ln(abundance) over time and used to estimate population exponential growth
#' rate. The generic method is fit_model(), which dispatches to growth model
#' specific fitting functions, ultimately these are to specific get.gr() functions.
#' 
#' @param model Growth model type
#' @param x Time variable
#' @param y ln(abundance)
#' @param \dots Additional arguments passed to fitting function (not used?)
#' 
#' @export
fit_model <- function(model, x, y, ...) {
  UseMethod("fit_model")
}

#' @rdname fit_model
#' @export
fit_model.linear_model <- function(model, x, y, ...) {
  get.gr.linear(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.lag_model <- function(model, x, y, ...) {
  get.gr.lag(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.sat_model <- function(model, x, y, ...) {
  get.gr.sat(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.flr_model <- function(model, x, y, ...) {
  get.gr.flr(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.lagsat_model <- function(model, x, y, ...) {
  get.gr.lagsat(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.satdecay_model <- function(model, x, y, ...) {
  get.gr.satdecay(x, y, ...)
}

#' @rdname fit_model
#' @export
fit_model.satdecay_ode_model <- function(model, x, y, ...) {
  get.gr.satdecay.ode(x, y, ...)
}



#' Evaluation methods for growth models
#' 
#' Generic for a suite of methods applied to different growth models to extract or 
#' calculate key metrics for evaluating the quality of model fits. This includes 
#' the exponential growth rate estimate itself, as well as measures for the 
#' fit in different portions of the model (e.g., initial exponential phase vs.
#' after saturating). Input follows a common structure, but internals vary by 
#' model type.
#' 
#' @param model Growth model type
#' @param fit Actual fit of growth model
#' @param x Time variable
#' @param y ln(abundance)
#' @param \dots Additional arguments (not used)
#' 
#' @export
evaluate_model <- function(model, fit, x, y, ...) {
  UseMethod("evaluate_model")
}

#' @rdname evaluate_model
#' @export
evaluate_model.linear_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  results<-list(
    coef = coef(fit),
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)[2]),
      se = unname(sqrt(diag(vcov(fit)))[2]),
      slope_n = length(x),
      slope_r2 = get.R2(preds, y),
      pre_n = NA,
      pre_r2 = NA,
      post_n = NA,
      post_r2 = NA,
      tmax = NA,
      nmax = NA
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.lag_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  b1 <- coef(fit)["B1"]
  
  exp_idx <- x >= (b1 - 0.1)
  pre_idx <- x <= b1

  results<-list(
    coef = coef(fit),
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)["b"]),
      se = unname(sqrt(diag(vcov(fit)))["b"]),
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = sum(pre_idx),
      pre_r2 = get.R2(preds[pre_idx], y[pre_idx]),
      post_n = NA,
      post_r2 = NA,
      tmax = NA,
      nmax = NA
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.sat_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  b2 <- coef(fit)["B2"]
  
  exp_idx <- x <= (b2 + 0.1)
  post_idx <- x >= b2
  
  results<-list(
    coef = coef(fit),
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)["b"]),
      se = unname(sqrt(diag(vcov(fit)))["b"]),
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = NA,
      pre_r2 = NA,
      post_n = sum(post_idx),
      post_r2 = get.R2(preds[post_idx], y[post_idx]),
      tmax = NA,
      nmax = NA
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.lagsat_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  b1 <- coef(fit)["B1"]
  b2 <- coef(fit)["B2"]
  
  exp_idx <- x >= (b1 - 0.1) & x <= (b2 + 0.1)
  pre_idx <- x <= b1
  post_idx <- x >= b2

  results<-list(
    coef = coef(fit),
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)["b"]),
      se = unname(sqrt(diag(vcov(fit)))["b"]),
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = sum(pre_idx),
      pre_r2 = get.R2(preds[pre_idx], y[pre_idx]),
      post_n = sum(post_idx),
      post_r2 = get.R2(preds[post_idx], y[post_idx]),
      tmax = NA,
      nmax = NA
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.flr_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  b2 <- coef(fit)["B2"]
  
  exp_idx <- x <= (b2 + 0.1)
  post_idx <- x >= b2
  
  results<-list(
    coef = coef(fit),
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)["b"]),
      se = unname(sqrt(diag(vcov(fit)))["b"]),
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = NA,
      pre_r2 = NA,
      post_n = sum(post_idx),
      post_r2 = get.R2(preds[post_idx], y[post_idx]),
      tmax = NA,
      nmax = NA
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.satdecay_model <- function(model, fit, x, y, ...) {
  
  preds <- predict(fit)
  
  cfs<-coef(fit)
  
  #b2 <- coef(fit)["B2"]
  b2 <- cfs$B2
  
  exp_idx <- x <= (b2 + 0.1)
  post_idx <- x >= b2
  
  results<-list(
    coef = cfs,
    preds = preds,
    metrics = list(
      slope = unname(coef(fit)["b"]),
      se = unname(sqrt(diag(vcov(fit)))["b"]),
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = NA,
      pre_r2 = NA,
      post_n = sum(post_idx),
      post_r2 = get.R2(preds[post_idx], y[post_idx]),
      tmax = cfs$B2,
      nmax = cfs$a + cfs$b*cfs$B2
    )
  )
  return(results)
}

#' @rdname evaluate_model
#' @export
evaluate_model.satdecay_ode_model <- function(model, fit, x, y, ...) {
  
  # to generate predictions from an mle2 model with implicit NLL calc,
  # need to compute from scratch.
  
  # extract and back transform coefficients
  tcoef<-function(cfs){
    cfs<-as.list(cfs)
    vec<-c(alpha=exp(cfs$log_alpha), vmax=exp(cfs$log_vmax), c=1 / (1 + exp(-cfs$theta_c)),
           d=exp(cfs$log_d), n0=cfs$n0, sigma=exp(cfs$log_sigma))
    as.list(vec)
  }
  cfs <- tcoef(coef(fit))

  # calculate predicted values using these coefficients from the fit:
  preds<-satdecay.ode(x,cfs$alpha,cfs$vmax,cfs$c,cfs$d,10,cfs$n0)
  
  # Numerical estimate of peak abundance
  derived <- derive.satdecay.stats(cfs,r0 = 10)
  #tmax.est<-satdecay.ode.peak.time(cfs)
  
  exp_idx <- x <= (derived$tmax + 0.1)
  post_idx <- x >= derived$tmax
  
  results<-list(
    coef = cfs,
    preds = preds,
    metrics = list(
      slope = cfs$vmax * (1 - cfs$d), 
      se = NA, #need to figure out calculation of se for this composite parameter
      slope_n = sum(exp_idx),
      slope_r2 = get.R2(preds[exp_idx], y[exp_idx]),
      pre_n = NA,
      pre_r2 = NA,
      post_n = sum(post_idx),
      post_r2 = get.R2(preds[post_idx], y[post_idx]),
      tmax = derived$tmax,
      nmax = derived$nmax
    )
  )
  return(results)
}


#' Plotting methods for growth models
#' 
#' Generic for a suite of methods used to plot the fits of different growth 
#' models. Specific methods have yet to be written.
#' 
#' @param model Growth model type
#' @param fit Actual fit of growth model
#' @param x Time variable
#' @param y ln(abundance)
#' @param \dots Additional arguments (not used)
#' 
#' @export
plot_model <- function(model, fit, x, y, ...) {
  UseMethod("plot_model")
}



#' Defining growth model classes
#' 
#' Each growth model exists as both a 'growth_model' object and carries a specific
#' class determining which type of growth model it is. These functions establish
#' these class definitions.
#' 
#' @export
linear_model <- function() {
  structure(
    list(name = "linear"),
    class = c("linear_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
lag_model <- function() {
  structure(
    list(name = "lag"),
    class = c("lag_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
sat_model <- function() {
  structure(
    list(name = "sat"),
    class = c("sat_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
flr_model <- function() {
  structure(
    list(name = "flr"),
    class = c("flr_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
lagsat_model <- function() {
  structure(
    list(name = "lagsat"),
    class = c("lagsat_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
satdecay_model <- function() {
  structure(
    list(name = "satdecay"),
    class = c("satdecay_model", "growth_model")
  )
}

#' @rdname linear_model
#' @export
satdecay_ode_model <- function() {
  # these models will, for now, require Julia. Make sure it's available when
  # setting up a model, before any fitting is even attempted.
  ensure_julia()
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("JuliaCall is required for fitting satdecay_ode model")
  }
  
  structure(
    list(name = "satdecayode"),
    class = c("satdecay_ode_model", "growth_model")
  )
}




#' Growth fit class creation
#' 
#' @param model Growth model type
#' @param fit Actual fit of growth model
#' @param data List containing the x (time) and y (ln(abundance) data used for fitting
#' @param results Growth model results, including diagnostics/metrics from evaluating fit
#' @param status Status of growth model fit and evaluation
#' @param error Error(s) arising when model fit or evaluation failed
#' @param growth.rate.valid Does the fit meet conditions for robust growth rate estimate?
#' 
#' @export
new_growth_fit <- function(model,fit,data,results,status='ok',error=NA,growth.rate.valid = FALSE){
  
  structure(
    list(
      model = model,
      fit = fit,
      data = data,
      results = results,
      status = status,
      error = error,
      growth.rate.valid = growth.rate.valid
    ),
    class = "growth_fit"
  )
}

#' Print method for growth_fit object
#' 
#' @param object Object of class growth_fit
#' @param \dots Additional arguments (not used)
#' 
#' @export
print.growth_fit <- function(object, ...){
  
  cat("Model:", object$model$name, "\n")
  cat("Growth rate valid:", object$growth.rate.valid, "\n")
  
  if (object$growth.rate.valid) {
    cat("Slope:", object$results$metrics$slope, "\n")
    cat("SE:", object$results$metrics$se, "\n")
    cat("Slope R2:", object$results$metrics$slope_r2, "\n")
  }
}

#' Generic predict method for growth_fit object
#' 
#' @param object Object of class growth_fit
#' @param newdata New data (if any) to use for making predictions
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @method predict growth_fit
predict.growth_fit <- function(object, newdata = NULL, ...) {
  
  # if not requesting predictions for new data (on time steps)
  if (is.null(newdata)) {
    return(object$results$preds)
  }
  
  # require data frame
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame (in predict.growth_fit())")
  }
  
  # check for viable fit:
  if (is.null(object$fit)) {
    return(rep(NA, nrow(newdata))) # if none, return NAs
  }
  
  # require x column
  if (!"x" %in% names(newdata)) {
    stop("newdata must contain column 'x' (in predict.growth_fit())")
  }
  
  # otherwise use model to make new predictions, given newdata
  model <- object$model
  predict_model(model, object, newdata, ...)
}


#' Prediction methods for growth models
#' 
#' A suite of growth models is available which can be fit to time series of 
#' ln(abundance) over time and used to estimate population exponential growth
#' rate. The generic method for generating predictions from these growth_fit
#' objects, predict_model(), dispatches to growth model specific prediction
#' functions. All of these are ultimately accessed through predict.growth_fit().
#' 
#' @param model Growth model type
#' @param object growth_fit object (produced by run_growth_model, must include results of evaluate_model)
#' @param newdata New data (if any) used for making predictions; must contain column named 'x' for time values
#' @param \dots Additional arguments passed to fitting function (not used?)
#' 
#' @export
predict_model <- function(model, object, newdata, ...) {
  UseMethod("predict_model")
}

#' @rdname predict_model
#' @export
predict_model.linear <- function(model, object, newdata, ...) {
  predict(object$fit, newdata = newdata)
}

#' @rdname predict_model
#' @export
predict_model.lag <- function(model, object, newdata, ...) {
  cfs <- object$results$coef
  #cfs <- as.list(coef(fit))
  lag(newdata$x, a  = cfs$a, b  = cfs$b, B1 = cfs$B1, s  = 1E-10)
}

#' @rdname predict_model
#' @export
predict_model.sat <- function(model, object, newdata, ...) {
  cfs <- object$results$coef
  #cfs <- as.list(coef(fit))
  sat(newdata$x, a  = cfs$a, b  = cfs$b, B2 = cfs$B2, s  = 1E-10)
}

#' @rdname predict_model
#' @export
predict_model.lagsat <- function(model, object, newdata, ...) {
  cfs <- object$results$coef
  #cfs <- as.list(coef(fit))
  lagsat(newdata$x, a  = cfs$a, b  = cfs$b, B1 = cfs$B1, B2 = cfs$B2, s  = 1E-10)
}

#' @rdname predict_model
#' @export
predict_model.flr <- function(model, object, newdata, ...) {
  cfs <- object$results$coef
  #cfs <- as.list(coef(fit))
  flr(newdata$x, a  = cfs$a, b  = cfs$b, B2 = cfs$B2, s  = 1E-10)
}

#' @rdname predict_model
#' @export
predict_model.satdecay <- function(model, object, newdata, ...) {
  cfs <- object$results$coef
  #cfs <- as.list(coef(fit))
  satdecay(newdata$x, a  = cfs$a, b  = cfs$b, b2 = cfs$b2, B2 = cfs$B2, s  = 1E-10)
}

#' @rdname predict_model
#' @export
predict_model.satdecay_ode_model <- function(model, object, newdata, ...) {
  cfs <- object$results$coef #extract_satdecay_coefs(fit)
  satdecay.ode(newdata$x, cfs$alpha, cfs$vmax, cfs$c, cfs$d, 10, cfs$n0)
}


#' Generic plotting method for growth_fit object
#' 
#' @param object Object of class growth_fit
#' @param newdata New data (if any) to use for making predictions
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @method plot growth_fit
plot.growth_fit <- function(object, main = NULL, savepath = NULL, ngrid = 200, ...){
  
  obs.x<-object$data$x
  obs.y<-object$data$y

  xr <- range(obs.x, na.rm = TRUE)
  
  grid <- seq(xr[1], xr[2], length.out = ngrid)
  
  preds <- predict(object, newdata = data.frame(x=grid))
  print(length(preds))
  print(preds)
  
  if (!is.null(savepath)) {
    grDevices::pdf(savepath)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  graphics::plot(obs.y ~ obs.x,
    xlab = "Time", ylab = "ln(abundance)", main = main)
  
  graphics::lines(preds~grid, col = "blue", lwd = 2)

  invisible(list(grid = grid, preds = preds))
}


#' Fit and evaluate specific growth rate model
#' 
#' This generic function fits a specific user-selected growth rate model to time
#' series data on ln(abundance). In addition to providing the final model fit, 
#' derived metrics are calculated after the model is fit (the evaluation step). 
#' 
#' @param model Desired model type
#' @param x Time steps
#' @param y ln(abundance)
#' @param min.exp.obs control parameter specifying the minimum number of observations that must fall within the estimated exponential phase in order to consider lag, sat, lagsat, and flr models; defaults to 3.
#' @param internal.r2.cutoff control parameter specifying the R2 criteria that may be applied to drop fits where the number of observations in the exponential portion is equal to 3. The default value of zero permits all fits of 3 obs to be considered.
#' 
#' @return This function returns a growth_fit object
#' 
#' @export
run_growth_model <- function(model,x,y,min.exp.obs = 3,internal.r2.cutoff = 0){
  
  # helper function in case of fit failure:
  fail_fit <- function(stage, err = NA, x=x, y=y){
    n <- length(x)

    new_growth_fit(
      model = model,
      fit = NA,
      data = list(x=x,y=y),
      results = list(
        coef = NA,
        preds = rep(NA,n),
        metrics = list(
          slope = NA,
          se = NA,
          slope_n = NA,
          slope_r2 = NA,
          pre_n = NA,
          pre_r2 = NA,
          post_n = NA,
          post_r2 = NA,
          tmax = NA,
          nmax = NA
        )
      ),
      status = stage,
      error = err,
      growth.rate.valid = FALSE
    )
  }
  
  # attempt fit:
  fit <- try(fit_model(model, x, y), silent = TRUE)
  
  if(inherits(fit, "try-error")) {
    return(fail_fit(stage="fit_failed", err=attr(fit, "condition"),x=x,y=y))
  }
  
  # If fit succeeded, can evaluate the model:
  results <- try(evaluate_model(model, fit, x, y), silent = TRUE)
  
  # Note: this is probably too severe? what if only some metrics can't be calculated?
  # this would blank all of them out in the event of any error...
  if(inherits(results, "try-error")) {
    return(fail_fit(stage = "evaluation_failed",err = attr(results, "condition"),x=x,y=y))
  }
  
  # check whether fit has properties that are sufficient for the growth rate 
  # estimate to be robust/valid; e.g., based on enough observations
  metrics <- results$metrics
  growth.rate.valid <- !(
    is.na(metrics$slope_n) ||
      metrics$slope_n < min.exp.obs ||
      (metrics$slope_n == min.exp.obs &&
         !is.na(metrics$slope_r2) &&
         metrics$slope_r2 < internal.r2.cutoff)
  )
  
  # return standard object.
  # will indicate that fit is ok, but growth rate may not be valid; see above
  new_growth_fit(
    model = model,
    fit = fit,
    data = list(x=x,y=y),
    results = results,
    status = "ok",
    error = NA,
    growth.rate.valid = growth.rate.valid
  )
}

#' Growth rate result class creation
#' 
#' @param models Set of growth rate models employed
#' @param successful Which model(s) fit successfully
#' @param best Growth rate model selected as the best model
#' @param ictab Model comparison results (information criteria table)
#' 
#' @export
new_growth_rate_result <- function(models,successful,best,ictab){
  
  structure(
    list(
      models = models,
      successful = successful,
      best = best,
      ictab = ictab
    ),
    class = "growth_rate_result"
  )
}

#' Print method for best model in the suite of growth rate results
#' 
#' @param object Object of class growth_rate_result
#' @param \dots Additional arguments (not used)
#' 
#' @export
print.growth_rate_result <- function(object, ...) {
  
  cat("====================================\n")
  cat("Growth Rate Analysis\n")
  cat("====================================\n\n")
  
  cat("Best Model:", object$best$model$name, "\n")
  cat("Slope:", object$best$results$metrics$slope, "\n")
  cat("SE:", object$best$results$metrics$se, "\n")
  cat("Slope R2:", object$best$results$metrics$slope_r2, "\n")
  cat("\n")
  
  print(object$ictab)
}

#' Summary method for suite of successful growth rate results
#' 
#' @param object Object of class growth_rate_result
#' @param \dots Additional arguments (not used)
#' 
#' @export
summary.growth_rate_result <- function(object, ...) {
  
  data.frame(
    model = sapply(object$successful, function(x) x$model$name),
    slope = sapply(object$successful, function(x) x$slope),
    se = sapply(object$successful, function(x) x$se),
    slope_n = sapply(object$successful, function(x) x$slope_n),
    slope_r2 = sapply(object$successful, function(x) x$slope_r2)
  )
}

#' Glance method for suite of successful growth rate results
#' 
#' @param object Object of class growth_rate_result
#' @param \dots Additional arguments (not used)
#' 
#' @export
glance.growth_rate_result <- function(object, ...) {
  
  tibble::tibble(
    best_model = object$best$model$name,
    slope = object$best$results$metrics$slope,
    se = object$best$results$metrics$se,
    slope_n = object$best$results$metrics$slope_n,
    slope_r2 = object$best$results$metrics$slope_r2,
    nmax = object$best$results$metrics$nmax
  )
}

#' Augment method for growth model fit
#' 
#' Provides observation-level predictions and residuals based on a model fit
#' 
#' @param object Object of class growth_fit
#' @param \dots Additional arguments (not used)
#' 
#' @export
augment.growth_fit <- function(object, ...) {
  
  tibble::tibble(
    x = object$x,
    y = object$y,
    fitted = predict(object$fit),
    residual = residuals(object$fit)
  )
}

#' Empty growth model construction
#' 
#' Given the (character string) name of a model, generate a model specification
#' object, e.g. "linear" becomes "linear_model()", the infrastructure for holding
#' a model of this type. This empty infrastructure is later populated when we
#' actually ask R to fit this model using the fit_model() command.
#' 
#' @param name Name of a growth model type
#' 
#' @export
make_model <- function(name) {
  
  switch(
    name,
    
    linear = linear_model(),
    lag = lag_model(),
    sat = sat_model(),
    flr = flr_model(),
    lagsat = lagsat_model(),
    satdecay = satdecay_model(),
    satdecayode = satdecay_ode_model(),
    
    stop(paste("Unknown model:", name))
  )
}



#' Extract exponential growth rate from a time series of ln(population abundance)
#' 
#' This meta-function takes a time series of abundance, and attempts to extract an 
#' estimate of exponential growth rate, using one or more of a suite of possible methods.
#' These methods allow for the possibility that exponential growth may lag or saturate,
#' or both, over the course of the time series. All selected methods are used to fit 
#' models to the time series. Subsequently, model comparison (based on AIC) is used to 
#' determine which model best fits the focal data.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param id Label corresponding to the population/strain/species of interest; used to determine the title and file name of saved plot, if any.
#' @param methods Must be a character vector containing one or more of \code{'linear'}, \code{'lag'}, \code{'sat'}, \code{'flr'}, or \code{'lagsat'}
#' @param model.selection control parameter to specify which IC metric to use in model selection; default is AICc, which corrects for small sample sizes and converges asymptotically on AIC.
#' @param min.exp.obs control parameter specifying the minimum number of observations that must fall within the estimated exponential phase in order to consider lag, sat, lagsat, and flr models; defaults to 3.
#' @param internal.r2.cutoff control parameter specifying the R2 criteria that may be applied to drop fits where the number of observations in the exponential portion is equal to 3. The default value of zero permits all fits of 3 obs to be considered.
#' @param plot.best.Q logical; should the best fitting model be plotted?
#' @param fpath character; if best model is to be plotted, provide the file path for saving the plot
#' @param zero.time if TRUE, shift time axis so that each time series starts at time = 0
#' 
#' @return A data frame containing the identity of the best model, the content of the best model, the estimated slopes of the increasing linear portion of the regressions (ie, exponential growth rate), the standard errors associated with these slopes, the IC table used to determine the best model, and the full list of all models fit. See vignette for details.
#' 
#' @examples
#' sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10),rep('D',10)),
#'                 dtime=rep(seq(1,10),4),
#'                 ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),
#'                            c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,
#'                            c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1),
#'                            c(5.5,4.5,3.5,2.5,1.5,0,0.2,-0.1,0,-0.1)))
#'                            
#' # for single replicate                            
#' sdat2<-sdat[sdat$trt=='A',]
#' 
#' # calculate growth rate using all available methods:
#' res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = TRUE,id = 'Population A')
#' res$best.model
#' res$best.slope
#' 
#' @export
get.growth.rate <- function(x,y,id,
    methods = c("linear","lag","sat","flr","lagsat","satdecay","satdecayode"),
    model.selection = "AICc",min.exp.obs = 3,internal.r2.cutoff = 0,
    plot.best.Q = TRUE, fpath = NA, zero.time = TRUE){
  
  # Clean data
  keep <- !is.na(y)
  
  x <- x[keep]
  y <- y[keep]
  
  if (zero.time) {
    x <- x - min(x, na.rm = TRUE)
  }
  
  if (length(unique(x)) < 2) {
    stop("Fewer than two unique time points")
  }
  

  # Set up models
  model_objects <- lapply(methods,make_model)
  

  # Fit all models
  results <- lapply(
    model_objects,
    run_growth_model,
    x = x,
    y = y,
    min.exp.obs = min.exp.obs,
    internal.r2.cutoff = internal.r2.cutoff
  )
  
  names(results) <- methods
  
  # which worked?
  successful <- Filter(
    function(x) {
      !is.null(x) &&
        isTRUE(x$growth.rate.valid)
    },
    results
  )
  
  if (length(successful) == 0) {
    stop("All requested models failed")
  }
  
  # Model comparison, using only models that successfully fit
  mod.list <- lapply(successful, function(x) x$fit)

  mod.names <- sapply(
    successful,
    function(x) x$model$name
  )
  
  ictab <- switch(
    model.selection,
    
    AIC = bbmle::AICtab(
      mod.list,
      mnames = mod.names
    ),
    
    AICc = bbmle::AICctab(
      mod.list,
      mnames = mod.names,
      nobs=length(x)
    ),
    
    BIC = bbmle::BICtab(
      mod.list,
      mnames = mod.names,
      nobs=length(x)
    ),
    
    stop("Invalid model.selection")
  )
  
  # extract best model
  best.name <- attr(ictab, "row.names")[1]
  
  best <- successful[
    sapply(
      successful,
      function(x) x$model$name == best.name
    )
  ][[1]]
  

  # Return results as object
  new_growth_rate_result(
    models = results,
    successful = successful,
    best = best,
    ictab = ictab
  )
}


# ============================================================
# EXAMPLE USAGE
# ============================================================

# result <- get.growth.rate(x, y)

# print(result)

# summary(result)

# result$best$slope
# result$best$model$name
# result$ictab
# result$successful