# methods for an object of class mixed param.MixedTS
setGeneric("Qparam.MixedTS",
           function(object,ret)
             standardGeneric("Qparam.MixedTS")
)

setMethod("Qparam.MixedTS","param.MixedTS",
          function(object,ret){
            if(object@Mixing=="Gamma"){
              res<-Qparam.MixedTS.aux(mu0=object@mu0,
                                 mu=object@mu, sig=object@sigma, a=object@a,
                                 alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                 ret=ret, Parametrization=object@Parametrization)
            }else{
              warning("The Q parameters are available only for gamma mixing density")
              return(NULL)
            }            
            return(res)
          }
)


Qparam.MixedTS.aux<-function(mu0, mu, sig, a, alpha, lambda_p, lambda_m, ret, Parametrization){
  # Find c^{*} that ensures the underlyng price process is a martingale
  dum.env<-new.env()
  assign("mu0", mu0, envir = dum.env)
  assign("mu", mu, envir = dum.env)
  assign("sig", sig, envir = dum.env)
  assign("a", a, envir = dum.env)
  assign("alpha", alpha, envir = dum.env)
  assign("lambda_p", lambda_p, envir = dum.env)
  assign("lambda_m", lambda_m, envir = dum.env)
  assign("r", ret, envir = dum.env)
  assign("Parametrization", Parametrization, envir = dum.env)
  
  dum.f<-function(par,env){
    cparam <- par
    mgf_c <- as.numeric(MixedTS:::charact.MTSgam(t= -1i*cparam,
             mu0 = env$mu0, mu = env$mu, sig = env$sig, a = env$a, 
             alpha = env$alpha, lambda_p = env$lambda_p,
             lambda_m = env$lambda_m, Parametrization=env$Parametrization))
    mgf_c1 <- as.numeric(MixedTS:::charact.MTSgam(t= -1i*(cparam+1),
              mu0 = env$mu0, mu = env$mu, sig = env$sig, a = env$a, 
              alpha = env$alpha, lambda_p = env$lambda_p,
              lambda_m = env$lambda_m, Parametrization=env$Parametrization))
    res <- ( log(mgf_c1) - log(mgf_c)-unique(env$r[[1]]))^2 
  }
  c_start <- 0
  res <- optim(fn = dum.f, par = c_start,lower=-(lambda_m-10^(-12)), upper=lambda_p-10^(-12),
               method ="L-BFGS-B",env = dum.env)
  c_star <- res$par
  lambda_pq <- lambda_p - c_star
  lambda_mq <- lambda_m + c_star
  mu_q <- (mu+((lambda_p^(alpha-1)-lambda_m^(alpha-1))-(lambda_pq^(alpha-1)-lambda_mq^(alpha-1)))
           /((alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2))))
  
#   mgf <- as.numeric(MixedTS:::charact.MTSgam(t= -1i*c_star,
#          mu0 =0, mu = mu, sig = sig, a = a, 
#          alpha = alpha, lambda_p = lambda_p,
#          lambda_m = lambda_m))
#  charctexp <- log(mgf)
  t <- -1i
  
  charctexp <- as.numeric(((lambda_p-1i*t)^alpha-lambda_p^alpha+(lambda_m+1i*t)^alpha-lambda_m^alpha)/(alpha*(alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
+(1i*t*(lambda_p^(alpha-1)-lambda_m^(alpha-1)))/((alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2))))
  
  sig_2q <- sig^2*(lambda_pq^(alpha-2)+lambda_mq^(alpha-2))/((1-c_star*mu-sig^2*charctexp)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
  sig_q <- sqrt(sig_2q)
  mu0_q <- mu0
  mu_q <- mu_q*(lambda_p^(alpha-2)+lambda_m^(alpha-2))/(lambda_pq^(alpha-2)+lambda_mq^(alpha-2))
  if(is.nan(sig_q)){sig_q <- 10^(-8)}
  res   <- setMixedTS.param(mu0 = mu0_q, mu = mu_q, 
           sigma = sig_q, a=a, alpha = alpha, 
           lambda_p = lambda_pq, lambda_m = lambda_mq)
  
  return(res)
}

# Class for data options

setClass("OptionData", 
         representation(UnderPrice = "numeric",
                        PriceOpt = "numeric",
                        ImpliedVol = "numeric",
                        TimeToMat = "numeric",
                        Strike = "numeric",
                        Type = "character",
                        rate = "numeric",
                        qyield = "numeric",
                        dateobs = "Date"),
                      
         prototype(UnderPrice = numeric(),
                   PriceOpt = numeric(),
                   ImpliedVol = numeric(),
                   Strike = numeric(),
                   Type = character(),
                   rate = numeric(),
                   TimeToMat = numeric(),
                   qyield = numeric(),
                   a<-Sys.Date())
         
)
DataOpt<-function(UnderPrice, PriceOpt, ImpliedVol, TimeToMat, 
                  Strike, type = "CALL", rate = 0, qyield = 0, Date = Sys.Date(),
                  YearBasis = 360){
  
  if((missing(PriceOpt) && missing(ImpliedVol))|| missing(Strike) || missing(TimeToMat)){
    warning("Option Data or Strike values or Time to Maturities are necessary")
    return(NULL)
  }
  if(!missing(ImpliedVol) && missing(PriceOpt) && missing(UnderPrice)){
    NumData <- length(ImpliedVol)
    UnderPrice <- rep(1, NumData)
  }
  if(missing(ImpliedVol) && !missing(PriceOpt) && missing(UnderPrice)){
    warning("Underlying price is necessary")
    return(NULL)
  }else{
    NumData <- length(PriceOpt)
  }
  if(length(Date==1)){
    Date <- rep(Date,NumData)
  }
  if(length(type)==1){
    type<-rep(type,NumData)
  }  
  if(length(rate)==1){
    rate <- rep(rate,NumData)
  }
  if(length(qyield)==1){
    qyield <- rep(qyield,NumData)
  }
  if(missing(PriceOpt)){
    PriceOpt <- OptBSprice(S = UnderPrice , X = Strike,
                           r = rate, expTime = TimeToMat/YearBasis,
                           sig = ImpliedVol, y = qyield, type = type)
  }
  if(missing(ImpliedVol)){
    ImpliedVol <- numeric()
  } 
  res<-new("OptionData",UnderPrice = UnderPrice, PriceOpt = PriceOpt,
           ImpliedVol = ImpliedVol, TimeToMat = TimeToMat,
           Strike = Strike, Type = type, rate = rate, qyield = qyield,
           dateobs = Date)
  
  return(res)
}

OptBSprice <- function(S, X, r, expTime, sig, y, type) {
  NumData<-length(type)
  res<-numeric(NumData)
  rf <- r-y
  d1 <- (log(S/X)+(rf+sig^2/2)*expTime)/sig*sqrt(expTime)
  d2 <- d1 - sig * sqrt(expTime)
  for(i in c(1:NumData)){
    if(type[i]=="CALL"){
      res[i] <- S[i]*pnorm(d1[i]) - X[i]*exp(-rf[i]*expTime[i])*pnorm(d2[i])
    }else{
      res[i] <- X[i]*exp(-rf[i]*expTime[i]) * pnorm(-d2[i]) - S[i]*pnorm(-d1[i])
    }
    return(res)
  }
}
