setGeneric("OptionPrice",
           function(object, S0, Strike, TimeToMat, ret = 0, 
                    yield = 0, basis = 360, QMeasure = TRUE, 
                    type = "CALL")
             standardGeneric("OptionPrice")
)

#
setMethod("OptionPrice","param.MixedTS",
          function(object, S0, Strike, TimeToMat, ret = 0, 
                   yield = 0, basis = 360, QMeasure =TRUE, type = "CALL"){
            if(object@Mixing=="Gamma"){
              if(QMeasure==TRUE){
                   res <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                          a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                          lambda_m = object@lambda_m, S0, Strike = Strike, TimeToMat =TimeToMat, 
                          yield = yield, basis = basis, ret = ret, QMeasure = QMeasure, type = type)
              }else{
                  object<-Qparam.MixedTS.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                          a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                          lambda_m = object@lambda_m, ret = ret)
                  
                  res <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                         a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                         lambda_m = object@lambda_m, S0, Strike = Strike, TimeToMat =TimeToMat, 
                         yield = yield, basis = basis, ret = ret, QMeasure = QMeasure, type = type)
              }
            }else{
              warning("The Q parameters are available only for gamma mixing density")
              return(NULL)
            }            
            return(res)
          }
)

OptionPrice.aux <- function(mu0, mu, sig, a, alpha, lambda_p, lambda_m,
                   S0, Strike, TimeToMat, yield, basis, ret, 
                   QMeasure, type){
  mu0Q <- mu0
  if(QMeasure==TRUE){
    # We choose mu0Q in order to ensure the martingale condition
    term1 <- MixedTS:::charact.MTSgam(t= 1i*1,
                             mu0 = mu0, mu = mu, sig = sig, a = a, 
                             alpha = alpha, lambda_p = lambda_p,
                             lambda_m = lambda_m)
    mu0Q <- ret - a*log(mu+term1)
  }
  # we have identify the model with mu0Q, mu, sig, a, alpha, lambda_p and lambda_m
  # INSERT FOURIER TRANSFORM
  dtvar <- 0.01
  tvar <- seq(-1,1, by = dtvar)
  logret <- tvar/(1-tvar^2)
  correct <- (1 + tvar^2)/(1 - tvar^2)^2
  # remove inf
  cond <- is.finite(logret)
  
  LogRet <- logret[cond]
  Correct <- correct[cond]
  Tvar <- tvar[cond]
  dens <- MixedTS:::dMixedTS.aux(x = LogRet, mu0 = mu0Q, mu = mu, sig = sig, a = a,
          alpha = alpha, lambda_p = lambda_p, lambda_m = lambda_m)
  if(type=="CALL"){
    payoff <- exp(-ret*TimeToMat/basis)*pmax(S0*exp(LogRet)-Strike,0)*Correct
  }else{
    payoff <- exp(-ret*TimeToMat/basis)*pmax(Strike-S0*exp(LogRet),0)*Correct
  }
  integr <- payoff*dens 
  
  length.int <- length(integr)
  integr.RX <- na.omit(integr[seq(1,length.int, by=3)+2])
  integr.mid <- na.omit(integr[seq(1,length.int, by=3)+1])
  integr.LX <- na.omit(integr[seq(1,length.int, by=3)])
  
  new.int <- integr[c(1:length(integr.RX))]
  
  elem.price <- dtvar*2/6*(integr.LX+integr.RX+4*integr.mid) 
  price <- sum(elem.price)
  return(price)
}