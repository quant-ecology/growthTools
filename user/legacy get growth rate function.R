

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
#' @param plot.best.Q logical; should the best fitting model be plotted?
#' @param fpath character; if best model is to be plotted, provide the file path for saving the plot
#' @param methods Must be a character vector containing one or more of \code{'linear'}, \code{'lag'}, \code{'sat'}, \code{'flr'}, or \code{'lagsat'}
#' @param id Label corresponding to the population/strain/species of interest; used to determine the title and file name of saved plot, if any.
#' @param model.selection control parameter to specify which IC metric to use in model selection; default is AICc, which corrects for small sample sizes and converges asymptotically on AIC.
#' @param min.exp.obs control parameter specifying the minimum number of observations that must fall within the estimated exponential phase in order to consider lag, sat, lagsat, and flr models; defaults to 3.
#' @param internal.r2.cutoff control parameter specifying the R2 criteria that may be applied to drop fits where the number of observations in the exponential portion is equal to 3. The default value of zero permits all fits of 3 obs to be considered.
#' @param verbose logical; display or suppress diagnostics
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
get.growth.rate<-function(x,y,id,plot.best.Q=F,fpath=NA,methods=c('linear','lag','sat','flr','lagsat','satdecay'),model.selection=c('AICc'),min.exp.obs=3,internal.r2.cutoff=0,verbose=FALSE,zero.time=TRUE){
  
  # thin vectors if abundance measure is NA
  x<-x[!is.na(y)]
  y<-y[!is.na(y)]
  
  if(zero.time){
    x <- x-min(x,na.rm = TRUE)
  }
  
  if(sum(methods %in% c('linear','lag','sat','flr','lagsat','satdecay'))==0){
    print('Error! None of the specified methods matched a currently implemented approach')
  }
  
  if(verbose){print(paste('data set id = ',id))}
  
  # Initialize empty data structures
  modlist<-list(gr=NA,gr.lag=NA,gr.sat=NA,gr.flr=NA,gr.lagsat=NA,gr.satdecay=NA)
  class(modlist$gr)<-class(modlist$gr.lag)<-class(modlist$gr.sat)<-class(modlist$gr.flr)<-class(modlist$gr.lagsat)<-class(modlist$gr.satdecay)<-'try-error'
  gr<-gr.lag<-gr.sat<-gr.flr<-gr.lagsat<-gr.satdecay<-NA
  slope.gr<-slope.gr.lag<-slope.gr.sat<-slope.gr.flr<-slope.gr.lagsat<-slope.gr.satdecay<-NA
  se.gr<-se.gr.lag<-se.gr.sat<-se.gr.flr<-se.gr.lagsat<-se.gr.satdecay<-NA
  
  slope.n.gr<-slope.n.gr.lag<-slope.n.gr.sat<-slope.n.gr.flr<-slope.n.gr.lagsat<-slope.n.gr.satdecay<-NA
  slope.r2.gr<-slope.r2.gr.lag<-slope.r2.gr.sat<-slope.r2.gr.flr<-slope.r2.gr.lagsat<-slope.r2.gr.satdecay<-NA
  
  pre.n.gr<-pre.n.gr.lag<-pre.n.gr.sat<-pre.n.gr.flr<-pre.n.gr.lagsat<-pre.n.gr.satdecay<-NA
  pre.r2.gr<-pre.r2.gr.lag<-pre.r2.gr.sat<-pre.r2.gr.flr<-pre.r2.gr.lagsat<-pre.r2.gr.satdecay<-NA
  
  post.n.gr<-post.n.gr.lag<-post.n.gr.sat<-post.n.gr.flr<-post.n.gr.lagsat<-post.n.gr.satdecay<-NA
  post.r2.gr<-post.r2.gr.lag<-post.r2.gr.sat<-post.r2.gr.flr<-post.r2.gr.lagsat<-post.r2.gr.satdecay<-NA
  
  if(length(unique(x))==2){
    print('Caution: only two unique time points, high risk of over-fitting. Methods other than "linear" are likely to fail')
  }
  
  # if there are two or more unique time points with data:
  if(length(unique(x))>=2){
    
    # Fill data structures for each given method, if requested by user:
    if('linear' %in% methods){
      gr<-get.gr(x,y) 
      slope.gr<-coef(gr)[[2]]
      slope.n.gr<-length(x)
      slope.r2.gr<-get.R2(predict(gr),y)
      se.gr<-sqrt(diag(stats::vcov(gr)))['x']
      modlist$gr<-gr
    }
    if('lag' %in% methods){
      gr.lag<-try(get.gr.lag(x,y))  
      if(prod(class(gr.lag)!='try-error')){
        b1.cutoff <- coef(gr.lag)[1]-0.1  # where does exponential phase begin?
        pds.lag <- predict(gr.lag) # predicted values
        obs.lag <- y # observed values
        
        slope.gr.lag <- unname(coef(gr.lag)['b'])
        se.gr.lag <- sqrt(diag(stats::vcov(gr.lag)))['b']
        
        slope.n.gr.lag <- length(x[x>=b1.cutoff])  # how many observations above cutoff
        slope.r2.gr.lag <- get.R2(pds.lag[x>=b1.cutoff],obs.lag[x>=b1.cutoff])
        
        pre.n.gr.lag<-length(x[x<=coef(gr.lag)['B1']])
        pre.r2.gr.lag <- get.R2(pds.lag[x<=coef(gr.lag)['B1']],obs.lag[x<=coef(gr.lag)['B1']])
        
        # if exponential portion is based on fewer than min.exp.obs observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.lag < min.exp.obs  | (slope.n.gr.lag==min.exp.obs & slope.r2.gr.lag<internal.r2.cutoff)){
          class(gr.lag)<-'try-error'          
        }
      }
      modlist$gr.lag<-gr.lag
    }
    if('sat' %in% methods){
      gr.sat<-try(get.gr.sat(x,y))
      if(prod(class(gr.sat)!='try-error')){
        b2.cutoff <- coef(gr.sat)[1]+0.1  # where does exponential phase end?
        pds.sat <- predict(gr.sat) # predicted values
        obs.sat <- y # observed values 
        
        slope.gr.sat <- unname(coef(gr.sat)['b'])
        se.gr.sat <- sqrt(diag(stats::vcov(gr.sat)))['b']
        
        slope.n.gr.sat <- length(x[x<=b2.cutoff])  # how many observations below cutoff
        slope.r2.gr.sat <- get.R2(pds.sat[x<=b2.cutoff],obs.sat[x<=b2.cutoff])
        
        post.n.gr.sat <- length(x[x>=coef(gr.sat)['B2']])
        post.r2.gr.sat <- get.R2(pds.sat[x>=coef(gr.sat)['B2']],obs.sat[x>=coef(gr.sat)['B2']])
        
        # if exponential portion is based on fewer than min.exp.obs observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.sat < min.exp.obs | (slope.n.gr.sat==min.exp.obs & slope.r2.gr.sat < internal.r2.cutoff)){
          class(gr.sat)<-'try-error'          
        }
      }
      modlist$gr.sat<-gr.sat
    }
    
    if('satdecay' %in% methods){
      gr.satdecay<-try(get.gr.satdecay(x,y))
      if(prod(class(gr.satdecay)!='try-error')){
        b2.cutoff <- coef(gr.satdecay)[1]+0.1  # where does exponential phase end?
        pds.satdecay <- predict(gr.satdecay) # predicted values
        obs.satdecay <- y # observed values 
        
        slope.gr.satdecay <- unname(coef(gr.satdecay)['b'])
        se.gr.satdecay <- sqrt(diag(stats::vcov(gr.satdecay)))['b']
        
        slope.n.gr.satdecay <- length(x[x<=b2.cutoff])  # how many observations below cutoff
        slope.r2.gr.satdecay <- get.R2(pds.satdecay[x<=b2.cutoff],obs.satdecay[x<=b2.cutoff])
        
        post.n.gr.satdecay <- length(x[x>=coef(gr.satdecay)['B2']])
        post.r2.gr.satdecay <- get.R2(pds.satdecay[x>=coef(gr.satdecay)['B2']],obs.satdecay[x>=coef(gr.satdecay)['B2']])
        
        # if exponential portion is based on fewer than min.exp.obs observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.satdecay < min.exp.obs | (slope.n.gr.satdecay==min.exp.obs & slope.r2.gr.satdecay < internal.r2.cutoff)){
          class(gr.satdecay)<-'try-error'          
        }
      }
      modlist$gr.satdecay<-gr.satdecay
    }
    
    if('flr' %in% methods){
      gr.flr<-try(get.gr.flr(x,y))
      if(prod(class(gr.flr)!='try-error')){
        b2.cutoff <- coef(gr.flr)['B2']+0.1  # where does exponential phase end?
        pds.flr <- predict(gr.flr) # predicted values
        obs.flr <- y # observed values
        
        slope.gr.flr <- unname(coef(gr.flr)['b'])
        se.gr.flr <- sqrt(diag(stats::vcov(gr.flr)))['b']
        
        slope.n.gr.flr <- length(x[x<=b2.cutoff])  # how many observations below cutoff
        slope.r2.gr.flr <- get.R2(pds.flr[x<=b2.cutoff],obs.flr[x<=b2.cutoff])
        
        post.n.gr.flr<-length(x[x>=coef(gr.flr)['B2']])
        post.r2.gr.flr<-get.R2(pds.flr[x>=coef(gr.flr)['B2']],obs.flr[x>=coef(gr.flr)['B2']])
        
        # if exponential portion is based on fewer than min.exp.obs observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.flr < min.exp.obs | (slope.n.gr.flr==min.exp.obs & slope.r2.gr.flr < internal.r2.cutoff)){
          class(gr.flr)<-'try-error'          
        }
      }
      modlist$gr.flr<-gr.flr
    }
    
    if('lagsat' %in% methods){
      gr.lagsat<-try(get.gr.lagsat(x,y))
      if(prod(class(gr.lagsat)!='try-error')){
        b1.cutoff <- coef(gr.lagsat)[1]-0.1  # where does exponential phase begin?
        b2.cutoff <- coef(gr.lagsat)[2]+0.1  # where does exponential phase end?
        pds.lagsat <- predict(gr.lagsat) # predictions 
        obs.lagsat <- y # observed values
        
        slope.gr.lagsat <- unname(coef(gr.lagsat)['b'])
        se.gr.lagsat <- sqrt(diag(stats::vcov(gr.lagsat)))['b']
        
        slope.n.gr.lagsat <- length(x[x<=b2.cutoff & x >=b1.cutoff])  # how many obs btwn cutoffs
        slope.r2.gr.lagsat <- get.R2(pds.lagsat[x<=b2.cutoff & x >=b1.cutoff],
                                     obs.lagsat[x<=b2.cutoff & x >=b1.cutoff])
        
        pre.n.gr.lagsat<-length(x[x<=coef(gr.lagsat)['B1']])
        pre.r2.gr.lagsat <- get.R2(pds.lagsat[x<=coef(gr.lagsat)['B1']],
                                   obs.lagsat[x<=coef(gr.lagsat)['B1']])
        
        post.n.gr.lagsat<-length(x[x>=coef(gr.lagsat)['B2']])
        post.r2.gr.lagsat <- get.R2(pds.lagsat[x>=coef(gr.lagsat)['B2']],
                                    obs.lagsat[x>=coef(gr.lagsat)['B2']])
        
        # if exponential portion is based on fewer than min.exp.obs observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.lagsat < min.exp.obs  | (slope.n.gr.lagsat==min.exp.obs & slope.r2.gr.lagsat < internal.r2.cutoff)){
          class(gr.lagsat)<-'try-error'          
        }
        
      }
      modlist$gr.lagsat<-gr.lagsat
    }
    
    # determine which fits occured and were successful
    successful.fits<-sapply(modlist,detect)
    if(verbose) print(successful.fits)
    
    if(sum(successful.fits)==0){
      print('Error! All results for requested methods failed!')
      stop('Error! All results for requested methods failed!')
    }
    
    # assemble model names, contents, and slopes, but only for successful fits
    mod.names<-c('gr','gr.lag','gr.sat','gr.flr','gr.lagsat','gr.satdecay')[successful.fits]
    mod.list<-list(gr,gr.lag,gr.sat,gr.flr,gr.lagsat,gr.satdecay)[successful.fits]
    slope.ests<-unname(c(slope.gr,slope.gr.lag,slope.gr.sat,slope.gr.flr,slope.gr.lagsat,slope.gr.satdecay)[successful.fits])
    se.ests<-unname(c(se.gr,se.gr.lag,se.gr.sat,se.gr.flr,se.gr.lagsat,se.gr.satdecay)[successful.fits])
    slope.n.vals<-c(slope.n.gr,slope.n.gr.lag,slope.n.gr.sat,slope.n.gr.flr,slope.n.gr.lagsat,slope.n.gr.satdecay)[successful.fits]
    slope.r2.vals<-c(slope.r2.gr,slope.r2.gr.lag,slope.r2.gr.sat,slope.r2.gr.flr,slope.r2.gr.lagsat,slope.r2.gr.satdecay)[successful.fits]
    pre.n.vals<-c(pre.n.gr,pre.n.gr.lag,pre.n.gr.sat,pre.n.gr.flr,pre.n.gr.lagsat,pre.n.gr.satdecay)[successful.fits]
    pre.r2.vals<-c(pre.r2.gr,pre.r2.gr.lag,pre.r2.gr.sat,pre.r2.gr.flr,pre.r2.gr.lagsat,pre.r2.gr.satdecay)[successful.fits]
    post.n.vals<-c(post.n.gr,post.n.gr.lag,post.n.gr.sat,post.n.gr.flr,post.n.gr.lagsat,post.n.gr.satdecay)[successful.fits]
    post.r2.vals<-c(post.r2.gr,post.r2.gr.lag,post.r2.gr.sat,post.r2.gr.flr,post.r2.gr.lagsat,post.r2.gr.satdecay)[successful.fits]
    
    # compare successful models
    switch(model.selection,
           AIC={ictab <- bbmle::AICtab(mod.list,mnames = mod.names)},
           AICc={ictab <- bbmle::AICctab(mod.list,mnames = mod.names)},
           BIC={ictab <- bbmle::BICtab(mod.list,mnames = mod.names)},
           print("Error! Invalid IC method selected in get.nbcurve()"))
    
    best.mod.id<-which(mod.names==attr(ictab,"row.names")[1])
    
    # impose QC based on slope.n and slope.r2 here? or outside of function...
    
    # match model names to outputs:
    names(slope.ests)<-mod.names
    names(se.ests)<-mod.names
    names(slope.n.vals)<-mod.names
    names(slope.r2.vals)<-mod.names
    names(pre.n.vals)<-mod.names
    names(pre.r2.vals)<-mod.names
    names(post.n.vals)<-mod.names
    names(post.r2.vals)<-mod.names
    
    # format output:
    result<-list(best.slope=slope.ests[[best.mod.id]],
                 best.se=se.ests[[best.mod.id]],
                 best.model=as.character(mod.names[[best.mod.id]]),
                 best.model.rsqr=get.R2(predict(mod.list[[best.mod.id]]),y),
                 best.model.slope.n=slope.n.vals[[best.mod.id]],
                 best.model.slope.r2=slope.r2.vals[[best.mod.id]],
                 best.model.pre.ns=pre.n.vals[[best.mod.id]],
                 best.model.pre.rs=pre.r2.vals[[best.mod.id]],
                 best.model.post.ns=post.n.vals[[best.mod.id]],
                 best.model.post.rs=post.r2.vals[[best.mod.id]],
                 best.model.contents=list(mod.list[[best.mod.id]]),
                 slopes=slope.ests,
                 ses=se.ests,
                 slope.ns=slope.n.vals,
                 slope.rs=slope.r2.vals,
                 ictab=ictab,
                 models=list(gr=gr,gr.lag=gr.lag,gr.sat=gr.sat,gr.flr=gr.flr,gr.lagsat=gr.lagsat,gr.satdecay=gr.satdecay))
    #print(result)
    
    if(plot.best.Q){
      if(!is.na(fpath)){
        fpath<-paste(fpath,id[1],'.pdf',sep='')
      }
      
      # want to show the best model in the requested model set... given methods options.
      # no model can end up in the model set if not requested, so this should be o.k. as written
      gigo<-switch(result$best.model,
                   gr=get.gr(x,y,plotQ=TRUE,fpath=fpath,id=id[1]),
                   gr.lag=get.gr.lag(x,y,plotQ=TRUE,fpath=fpath,id=id[1]),
                   gr.sat=get.gr.sat(x,y,plotQ=TRUE,fpath=fpath,id=id[1]),
                   gr.flr=get.gr.flr(x,y,plotQ=TRUE,fpath=fpath,id=id[1]),
                   gr.lagsat=get.gr.lagsat(x,y,plotQ=TRUE,fpath=fpath,id=id[1]),
                   gr.satdecay=get.gr.satdecay(x,y,plotQ=TRUE,fpath=fpath,id=id[1]))
    }
  }else{
    print("Warning: fewer than two unique time points provided")
    
    result<-list(best.slope=NA,
                 best.se=NA,
                 best.model="NA",
                 best.model.rsqr=NA,
                 best.model.slope.n=NA,
                 best.model.slope.r2=NA,
                 best.model.pre.ns=NA,
                 best.model.pre.rs=NA,
                 best.model.post.ns=NA,
                 best.model.post.rs=NA,
                 best.model.contents=list(NA),
                 slopes=NA,
                 ses=NA,
                 slope.ns=NA,
                 slope.rs=NA,
                 models=list(gr=NA,gr.lag=NA,gr.sat=NA,gr.flr=NA,gr.lagsat=NA,gr.satdecay=NA))
  }
  
  return(result)
}

