# Write Methods and Classes for Calibration.

# Class 'calibrates.MixedTS'


setClass("calibrates", 
         representation(coef = "numeric",
                        vcov = "matrix",
                        obj = "numeric"
         ),
)

setClass("calibrates.MixedTS", 
         contains = c("calibrates","param.MixedTS")
)

# Method 'calibrate'

setGeneric("calibrate",
           function(model, Data, QMeasure = FALSE, basis = 360, method = "RMSE")
             standardGeneric("calibrate")
)

setMethod("calibrate","param.MixedTS",
          function(model, Data, QMeasure = FALSE,  basis = 360, method = "RMSE"){
            call <- match.call()
            res <- calibrateaux(mu0 = model@mu0, mu = model@mu, sig = model@sigma, a = model@a, 
                    alpha = model@alpha, lambda_p = model@lambda_p, lambda_m = model@lambda_m,
                    Price = Data@PriceOpt, S0 = Data@UnderPrice, Strike = Data@Strike, 
                    TimeToMat = Data@TimeToMat, yield = Data@qyield , basis = basis, 
                    ret = Data@rate, QMeasure = QMeasure, type = Data@Type, 
                    dateobs = Data@Date, method = method, call=call, 
                    Parametrization=model@Parametrization)
          }
)

calibrateaux <- function(mu0, mu, sig, a, alpha, lambda_p, lambda_m,
                          Price, S0, Strike, TimeToMat, yield, basis, 
                          ret, QMeasure, type, dateobs, method, call, Parametrization){

  if(method %in% c("RMSE","RMSPE","AAE","ARPE", "APE","APPE","MPE1")){
    param <- c(mu0, mu, sig, a, alpha, lambda_p, lambda_m)
    names(param) <- c("mu0", "mu", "sigma", "a", "alpha", "lambda_p", "lambda_m")
  }
  
  if(method %in% c("MPE2")){
    sig_err <- 1
    param <- c(mu0, mu, sig, a, alpha, lambda_p, lambda_m, sig_err)
    names(param) <- c("mu0", "mu", "sigma", "a", 
      "alpha", "lambda_p", "lambda_m", "sig_err")
  }
  
  if(method %in% c("AGNE","APGNE","AMGNE")){
    beta_err <- 2
    gamma_err <- 1
    param <- c(mu0, mu, sig, a, alpha, lambda_p, lambda_m, beta_err, gamma_err)
    names(param) <- c("mu0", "mu", "sigma", "a", 
                      "alpha", "lambda_p", "lambda_m", 
                      "beta_err", "gamma_err")
  }
  
  if(method %in% c("AMGNE")){
    beta_err <- 2
    gamma_err <- 1
    rho_err <- 0
    param <- c(mu0, mu, sig, a, alpha, lambda_p, lambda_m, beta_err, gamma_err, rho_err)
    names(param) <- c("mu0", "mu", "sigma", "a", 
                      "alpha", "lambda_p", "lambda_m", 
                      "beta_err", "gamma_err", "rho_err")
  }
  
  Maturity<-unique(TimeToMat)
  
  S0list <- list()
  Pricelist <- list()
  Strikelist <- list()
  TimeToMatlist <- list()
  yieldlist <- list()
  basislist <- list()
  retlist <- list()
  QMeasurelist <- list()
  typelist <- list()
  for(i in 1:length(Maturity)){ 
    cond <- (TimeToMat== Maturity[i])
    # we build a list for each maturity
    S0list[[i]] <- S0[cond]
    Pricelist[[i]] <- Price[cond]
    Strikelist[[i]] <- Strike[cond]
    TimeToMatlist[[i]] <- TimeToMat[cond]
    yieldlist[[i]] <- yield[cond] 
    basislist[[i]] <- basis
    retlist[[i]] <- ret[cond]
    QMeasurelist[[i]] <- QMeasure
    typelist[[i]] <- type[cond]
  }
  
  
  my.envir <- new.env()
  assign("S0", S0list, envir = my.envir)
  assign("Price", Pricelist, envir = my.envir)
  assign("Strike", Strikelist, envir = my.envir)
  assign("TimeToMat", TimeToMatlist, envir = my.envir)
  assign("yield", yieldlist, envir = my.envir)
  assign("basis", basislist, envir = my.envir)
  assign("ret", retlist, envir = my.envir)
  assign("QMeasure", QMeasurelist, envir = my.envir)
  assign("type", typelist, envir = my.envir)
  assign("method", method, envir = my.envir)
  assign("Parametrization", Parametrization, envir = my.envir)
  
  #   if(method=="mle")
  #     assign("method", "Standard", envir = my.envir)
  
  #   res.partial <- optim(par=param, fn=objectfn, method = "L-BFGS-B", 
  #                      lower = c(-1000, -1000, 0.0001, 0.0001, 0.001, 0.001, 0.001 ),
  #                      upper = c(1000, 1000, 100, 100, 2, 100, 100), my.envir=my.envir)
  
    ui<-matrix(c(0,0,1,0,0,0,0),1,7)
    ui<-rbind(ui,
              c(0,0,0,1,0,0,0),
              c(0,0,0,0,1,0,0),
              c(0,0,0,0,-1,0,0),
              c(0,0,0,0,0,1,0),
              c(0,0,0,0,0,0,1))
    ci <- c(0.0001, 0.0001, 0.0001, -2, 0.0001, 0.0001)
  
    if(method %in% c("MPE2")){
      ui <- cbind(ui, rep(0,dim(ui)[1]))
      ui <- rbind(ui, c(0,0,0,0,0,0,0,1))
      ui <- rbind(ui, c(0,0,0,0,0,0,0,-1))
      ci <- c(ci,0.000001,-10)
  #     param <- c(param, sig_eps=1)
    }
  
  if(method %in% c("AGNE","APGNE","AMGNE")){
    ui <- cbind(ui, rep(0,dim(ui)[1]))
    ui <- cbind(ui, rep(0,dim(ui)[1]))
    ui <- rbind(ui, c(0,0,0,0,0,0,0,1,0))
    ci <- c(ci,0)
    ui <- rbind(ui, c(0,0,0,0,0,0,0,0,1))
    ci <- c(ci,0)
    if(method %in% c("AMGNE")){
      ui <- cbind(ui, rep(0,dim(ui)[1]))
      ui <- rbind(ui, c(0,0,0,0,0,0,0,0,0,1))
      ui <- rbind(ui, c(0,0,0,0,0,0,0,0,0,-1))
      ci <- c(ci,0-10^(-10))
      ci <- c(ci,-1-10^(-10)) 
    }
  }
  
  respartial <-constrOptim(theta=param, f = objectfn, 
                            grad =NULL, ui = ui, ci = ci, 
                            my.envir=my.envir)

  
  if(method %in% c("APE","APPE","MPE1","MPE2","AGNE","APGNE","AMGNE")){
     call <- call
     fullcoef<-coef <- respartial$par
     assign("method", method, envir = my.envir)
     HessMatr <- NA 
     ApproachLog <- TRUE
     HessMatr <- tryCatch(optimHess(par = respartial$par, fn = loglikCalib, my.envir=my.envir, ApproachLog=ApproachLog), error = function(par) NA )
#      if(unique(unlist(my.envir$QMeasure))){
#       dummyHessMatr1 <- HessMatr[-1,]
#       dummyHessMatr <- dummyHessMatr1[,-1]
#       #vcov <- matrix(NA, 7,7)
#       vcov[2:7,2:7] <- tryCatch(solve(dummyHessMatr), error = function(par) matrix(NaN,6,6) )
#     }else{
      vcov <- tryCatch(solve(HessMatr), error = function(par) matrix(NaN,dim(HessMatr)[1],dim(HessMatr)[1]) )
#     }
#     assign("method", "mle2", envir = my.envir)
     ApproachLog <- FALSE
     dummy <- loglikCalib(param = respartial$par, my.envir = my.envir, ApproachLog = ApproachLog)
     min <- dummy$res.mle 
     if(method %in% c("APE","APPE","MPE1","MPE2")){
        fullcoef<-c(fullcoef, mu_eps=dummy$meanErr, sd_eps = dummy$stdErr)
     }else{
       if(method %in% c("AGNE","APGNE")){
         fullcoef<-c(fullcoef, beta_err=dummy$beta_err, gamma_err = dummy$gamma_err)
       }
       if(method %in% c("AMGNE")){
         fullcoef<-c(fullcoef, beta_err=dummy$beta_err, gamma_err = dummy$gamma_err, rho_err = dummy$rho_err)
       }
     }
  
      # "AGNE","APGNE","AMGNE"  

    minuslogl <- objectfn
    nobs <- length(Price)
    method <- "L-BFGS-B"
    details <- list()

   
    respartialA <- new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef), 
                      vcov = vcov, min = min, details = details, minuslogl = minuslogl, 
                      nobs = nobs, method = method)
    return(respartialA)
 }
  
  return(respartial)
  
#   res <- auxcalibrate(mu0, mu, sig, a, alpha, lambda_p, lambda_m,
#                            Price, S0, Strike, TimeToMat, yield, basis, 
#                            ret, QMeasure, type, dateobs, method, call, Parametrization)
}


# auxcalibrate <- function(mu0, mu, sig, a, alpha, lambda_p, lambda_m,
#                          Price, S0, Strike, TimeToMat, yield, basis, 
#                          ret, QMeasure, type, dateobs, method, call, Parametrization){
#   
#   
#   if(method=="mle"){
#     call <- call
#     fullcoef<-coef <- res.partial$par
#     assign("method", "mle", envir = my.envir)
#     HessMatr <- NA 
#     HessMatr <- tryCatch(optimHess(par = res.partial$par, fn = objectfn, my.envir=my.envir), error = function(par) NA )
#     if(unique(unlist(my.envir$QMeasure))){
#       dummyHessMatr1 <- HessMatr[-1,]
#       dummyHessMatr <- dummyHessMatr1[,-1]
#       vcov <- matrix(NA, 7,7)
#       vcov[2:7,2:7] <- tryCatch(solve(dummyHessMatr), error = function(par) matrix(NaN,6,6) )
#     }else{
#       vcov <- tryCatch(solve(HessMatr), error = function(par) matrix(NaN,7,7) )
#     }
#     assign("method", "mle2", envir = my.envir)
#     dummy <- objectfn(param = res.partial$par, my.envir = my.envir)
#     min <- dummy$res.mle 
#     fullcoef<-c(fullcoef, mu_eps=dummy$meanErr, sd_eps = dummy$stdErr)
#     minuslogl <- objectfn
#     nobs <- length(Price)
#     method <- "L-BFGS-B"
#     details <- list()
#     
#     res <- new("mle", call = call, 
#                coef = coef, fullcoef = unlist(fullcoef), 
#                vcov = vcov, min = min, 
#                details = details, minuslogl = minuslogl, 
#                nobs = nobs, method = method)
#   }else{
#     res <- res.partial 
#   }
#   return(res)
#   
# }








objectfn <- function(param, my.envir, ApproachLog){
  
  if(my.envir$QMeasure==TRUE){
    
    res1 <- OptionPrice.aux(mu0 = param["mu0"], mu = param["mu"], sig = param["sigma"], 
                           a = param["a"], alpha = param["alpha"], lambda_p = param["lambda_p"],
                           lambda_m = param["lambda_m"], S0 = my.envir$S0, Strike = my.envir$Strike, 
                           TimeToMat = my.envir$TimeToMat, yield = my.envir$yield, basis = my.envir$basis,
                           ret = my.envir$ret, QMeasure = my.envir$QMeasure, type = my.envir$type, 
                           Parametrization = my.envir$Parametrization)
  }else{
    
    object <- Qparam.MixedTS.aux(mu0 = param["mu0"], mu = param["mu"], sig = param["sigma"], 
                                 a = param["a"], alpha = param["alpha"], lambda_p = param["lambda_p"],
                                 lambda_m = param["lambda_m"], ret = my.envir$ret[1],
                                 Parametrization = my.envir$Parametrization)
    NumbofPrice <- length(my.envir$Strike)
   res1<-numeric(length=NumbofPrice)
#     for(i in c(1:NumbofPrice)){
#       res1[i] <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
#                               a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
#                               lambda_m = object@lambda_m, S0 = my.envir$S0[i], Strike = my.envir$Strike[i], 
#                               TimeToMat = my.envir$TimeToMat[i], yield = my.envir$yield[i], basis = my.envir$basis,
#                               ret = my.envir$ret[i], QMeasure = my.envir$QMeasure, type = my.envir$type[i])
#     }
   
    res1 <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
      a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
      lambda_m = object@lambda_m, S0 = my.envir$S0, Strike = my.envir$Strike, 
      TimeToMat = my.envir$TimeToMat, yield = my.envir$yield, 
      basis = my.envir$basis, ret = my.envir$ret, 
      QMeasure = my.envir$QMeasure, type = my.envir$type, 
      Parametrization = my.envir$Parametrization)
    
    
  }
  
  if(my.envir$method == "RMSE" || my.envir$method == "APE"){
    res.mle <- sum((unlist(res1)- unlist(my.envir$Price))^2, na.rm = TRUE)
  }

  
  if(my.envir$method == "RMSPE" || my.envir$method == "APPE"){
    res.mle <- sum(((unlist(res1) - unlist(my.envir$Price))/unlist(my.envir$Price))^2, na.rm = TRUE)
  }

  if(my.envir$method == "ARPE"){
    res.mle <- sum(abs((unlist(res1) - unlist(my.envir$Price))/unlist(my.envir$Price)), na.rm = TRUE)
  }

  if(my.envir$method == "AAE"){
    res.mle <- sum(abs((unlist(res1) - unlist(my.envir$Price))), na.rm = TRUE)
  }

#   if(my.envir$method == "APE"){
#     error <- unlist(res1)-unlist(my.envir$Price)
#     stderr <- sd(x=error, na.rm = TRUE)
#     param["sig_err"] <- stderr
#   #  meanErr<- mean(x=error, na.rm=TRUE)
#     if(param["sig_err"]!=0){
#     res.mle <- -sum(na.omit(dnorm(x = error,
#       mean = 0, sd = stderr, log = TRUE)))
#     }else{res.mle <- -Inf}
#   }

#   if(my.envir$method == "APPE"){
#     sig_err <- param["sig_err"]
#     PriceTrue <- unlist(my.envir$Price)
#     error <- unlist(res1) - PriceTrue
#    # if(stdErr!=0){
#      res.mle <- -sum(na.omit(dnorm(x = error,
#             mean = 0, sd = PriceTrue*sig_err, 
#             log = TRUE)))
#     #}else{res.mle <- Inf}
#   }

#   if(my.envir$method == "MPE1"){
#     sig_err <- param["sig_err"]
#     PriceTrue <- log(unlist(my.envir$Price))
#     error <- log(unlist(res1)) - PriceTrue
#     res.mle <- -sum(na.omit(dnorm(x = error,
#        mean = 0, sd = sig_err, log = TRUE)))
#   }

if(my.envir$method == "MPE1"){
  res.mle <- sum((log(unlist(res1))- log(unlist(my.envir$Price)))^2, na.rm = TRUE)
}


  if(my.envir$method == "MPE2"){
    sig_err <- param["sig_err"]
    PriceTrue <- log(unlist(my.envir$Price))
    error <- log(unlist(res1)) - PriceTrue
    res.mle <- -sum(na.omit(dnorm(x = error,
        mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
  }

# "AGNE","APGNE","AMGNE"  

if(my.envir$method == "AGNE"){
  beta_err <- param["beta_err"]
  gamma_err <- param["gamma_err"]
  PriceTrue <- unlist(my.envir$Price)
  error <- unlist(res1) - PriceTrue
  #res.mle <- -sum(na.omit(dnorm(x = error,
  # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
  Numb <- length(PriceTrue)
  res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)-log(gamma(1/beta_err)))
               -sum((abs(error)/gamma_err)^beta_err))
}

if(my.envir$method == "APGNE"){
  beta_err <- param["beta_err"]
  gamma_err <- param["gamma_err"]
  PriceTrue <- unlist(my.envir$Price)
  
  
  error <- (unlist(res1)-PriceTrue)/PriceTrue
  #res.mle <- -sum(na.omit(dnorm(x = error,
  # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
  Numb <- length(PriceTrue)
  res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)
                     -log(gamma(1/beta_err)))-sum(log(PriceTrue), na.rm =TRUE)-sum((abs(error)/gamma_err)^beta_err, na.rm =TRUE))
}

if(my.envir$method == "AMGNE"){
  beta_err <- param["beta_err"]
  gamma_err <- param["gamma_err"]
  rho_err <- param["rho_err"]
  PriceTrue <- unlist(my.envir$Price)
  error <- (unlist(res1)-PriceTrue)/(PriceTrue^rho_err)
  #res.mle <- -sum(na.omit(dnorm(x = error,
  # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
  Numb <- length(PriceTrue)
  res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)-log(gamma(1/beta_err)))-sum(log(PriceTrue), na.rm =TRUE)*rho_err-sum((abs(error)/gamma_err)^beta_err, na.rm =TRUE))
}
#   if(my.envir$method == "mle2"){
#     error <- unlist(res1)-unlist(my.envir$Price)
#     stdErr <- sd(x=error, na.rm = TRUE)
#     meanErr<- mean(x=error, na.rm=TRUE)
#     if(stdErr!=0){
#       res.mle <- -sum(na.omit(dnorm(x = error,
#                                     mean = meanErr, sd = stdErr, log = TRUE)))
#     }else{res.mle <- Inf}
#     res.mle<-list(res.mle=res.mle, meanErr=meanErr, stdErr=stdErr)
#   }
  
  return(res.mle)
}

loglikCalib <-function(param, my.envir, ApproachLog){
  if(my.envir$QMeasure==TRUE){
    
    res1 <- OptionPrice.aux(mu0 = param["mu0"], mu = param["mu"], sig = param["sigma"], 
                            a = param["a"], alpha = param["alpha"], lambda_p = param["lambda_p"],
                            lambda_m = param["lambda_m"], S0 = my.envir$S0, Strike = my.envir$Strike, 
                            TimeToMat = my.envir$TimeToMat, yield = my.envir$yield, basis = my.envir$basis,
                            ret = my.envir$ret, QMeasure = my.envir$QMeasure, type = my.envir$type, 
                            Parametrization = my.envir$Parametrization)
  }else{
    
    object <- Qparam.MixedTS.aux(mu0 = param["mu0"], mu = param["mu"], sig = param["sigma"], 
                                 a = param["a"], alpha = param["alpha"], lambda_p = param["lambda_p"],
                                 lambda_m = param["lambda_m"], ret = my.envir$ret[1],
                                 Parametrization = my.envir$Parametrization)
    NumbofPrice <- length(my.envir$Strike)
    res1 <- numeric(length=NumbofPrice)
    #     for(i in c(1:NumbofPrice)){
    #       res1[i] <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
    #                               a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
    #                               lambda_m = object@lambda_m, S0 = my.envir$S0[i], Strike = my.envir$Strike[i], 
    #                               TimeToMat = my.envir$TimeToMat[i], yield = my.envir$yield[i], basis = my.envir$basis,
    #                               ret = my.envir$ret[i], QMeasure = my.envir$QMeasure, type = my.envir$type[i])
    #     }
    
    res1 <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                            a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                            lambda_m = object@lambda_m, S0 = my.envir$S0, Strike = my.envir$Strike, 
                            TimeToMat = my.envir$TimeToMat, yield = my.envir$yield, 
                            basis = my.envir$basis, ret = my.envir$ret, 
                            QMeasure = my.envir$QMeasure, type = my.envir$type, 
                            Parametrization = my.envir$Parametrization)
    
    
  }
  
    if(my.envir$method == "APE"){
      error <- unlist(res1)-unlist(my.envir$Price)
      stdErr <- sd(x=error, na.rm = TRUE)
      meanErr<- mean(x=error, na.rm=TRUE)
      if(stdErr!=0){
#         res.mle <- -sum(na.omit(dnorm(x = error,
#                                       mean = 0, sd = stdErr, log = TRUE)))
        Numb<-length(error)
        gamma_err <- stdErr*sqrt(2)
        res.mle <- -(Numb*(log(2)-log(2*gamma_err)
                           -log(gamma(1/2)))-sum((error/gamma_err)^2, na.rm =TRUE))
        
      }else{res.mle <- -Inf}
      if(ApproachLog==TRUE){
        res.mle <- res.mle
      }else{
        res.mle<-list(res.mle=res.mle, meanErr=meanErr, stdErr=stdErr)
      }
    }
  
    if(my.envir$method == "APPE"){
      realPrice <- unlist(my.envir$Price)
      error <- (unlist(res1)/realPrice - 1)
      stdErr <- sd(x=error, na.rm = TRUE)
      meanErr<- mean(x=error, na.rm=TRUE)
      if(stdErr!=0){
        res.mle <- -sum(na.omit(dnorm(x = error,
        mean = 0, sd = stdErr*realPrice, log = TRUE)))
      }else{res.mle <- -Inf}
      if(ApproachLog==TRUE){
      res.mle <- res.mle
      }else{
        res.mle<-list(res.mle=res.mle, meanErr=meanErr, stdErr=stdErr)
      }
    }

  if(my.envir$method == "MPE1"){
    error <- log(unlist(res1))-log(unlist(my.envir$Price))
    stdErr <- sd(x=error, na.rm = TRUE)
    meanErr<- mean(x=error, na.rm=TRUE)
    if(stdErr!=0){
      res.mle <- -sum(na.omit(dnorm(x = error,
                                    mean = 0, sd = stdErr, log = TRUE)))
    }else{res.mle <- -Inf}
    if(ApproachLog==TRUE){
      res.mle <- res.mle
    }else{
      res.mle<-list(res.mle=res.mle, meanErr=meanErr, stdErr=stdErr)
    }
  }
  if(my.envir$method == "MPE2"){
    sig_err <- param["sig_err"]
    PriceTrue <- log(unlist(my.envir$Price))
    error <- log(unlist(res1)) - PriceTrue
    res.mle <- -sum(na.omit(dnorm(x = error,
                                  mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
    if(ApproachLog==TRUE){
      res.mle <- res.mle
    }else{
      res.mle<-list(res.mle=res.mle, meanErr=-0.5*sig_err^2, stdErr=sig_err)
    }
  }

  if(my.envir$method == "AGNE"){
    beta_err <- param["beta_err"]
    gamma_err <- param["gamma_err"]
    PriceTrue <- unlist(my.envir$Price)
    error <- unlist(res1) - PriceTrue
    #res.mle <- -sum(na.omit(dnorm(x = error,
    # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
    Numb <- length(PriceTrue)
    res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)-log(gamma(1/beta_err)))-sum((abs(error)/gamma_err)^beta_err, na.rm =TRUE))
    if(ApproachLog==TRUE){
      res.mle <- res.mle
    }else{
      res.mle<-list(res.mle=res.mle, beta_err=beta_err, gamma_err=gamma_err)
    }
    
  }
  
  if(my.envir$method == "APGNE"){
    beta_err <- param["beta_err"]
    gamma_err <- param["gamma_err"]
    PriceTrue <- unlist(my.envir$Price)
    error <- (unlist(res1)-PriceTrue)/PriceTrue
    #res.mle <- -sum(na.omit(dnorm(x = error,
    # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
    Numb <- length(PriceTrue)
    res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)
                       -log(gamma(1/beta_err)))-sum(log(PriceTrue), na.rm =TRUE)-sum((abs(error)/gamma_err)^beta_err, na.rm =TRUE))
    if(ApproachLog==TRUE){
      res.mle <- res.mle
    }else{
      res.mle<-list(res.mle=res.mle, beta_err=beta_err, gamma_err=gamma_err)
    }
  }
  
  if(my.envir$method == "AMGNE"){
    beta_err <- param["beta_err"]
    gamma_err <- param["gamma_err"]
    rho_err <- param["rho_err"]
    PriceTrue <- unlist(my.envir$Price)
    error <- (unlist(res1)-PriceTrue)/PriceTrue^rho_err
    #res.mle <- -sum(na.omit(dnorm(x = error,
    # mean = -0.5*sig_err^2, sd = sig_err, log = TRUE)))
    Numb <- length(PriceTrue)
    res.mle <- -(Numb*(log(beta_err)-log(2*gamma_err)
                       -log(gamma(1/beta_err)))-sum(log(PriceTrue)*rho_err, na.rm =TRUE)-sum((abs(error)/gamma_err)^beta_err, na.rm =TRUE))
    if(ApproachLog==TRUE){
      res.mle <- res.mle
    }else{
      names(res.mle)<-""
      res.mle<-list(res.mle=res.mle, 
                    beta_err=beta_err, 
                    gamma_err=gamma_err, 
                    rho_err=rho_err)
    }
  }
  
  return(res.mle)
}
