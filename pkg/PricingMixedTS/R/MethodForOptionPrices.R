setGeneric("OptionPrice",
           function(object, S0, Strike, TimeToMat, ret = 0, 
                    yield = 0, basis = 360, QMeasure = FALSE, 
                    type = "CALL")
             standardGeneric("OptionPrice")
)

#
setMethod("OptionPrice","param.MixedTS",
          function(object, S0, Strike, TimeToMat, ret = 0, 
                   yield = 0, basis = 360, QMeasure = FALSE, type = "CALL"){
            if(object@Mixing=="Gamma"){
              if(QMeasure==TRUE){
                   res <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                          a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                          lambda_m = object@lambda_m, S0, Strike = Strike, TimeToMat =TimeToMat, 
                          yield = yield, basis = basis, ret = ret, QMeasure = QMeasure, type = type,
                          Parametrization=object@Parametrization)
              }else{
                  object<-Qparam.MixedTS.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                          a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                          lambda_m = object@lambda_m, ret = ret, Parametrization=object@Parametrization)
                  
                  res <- OptionPrice.aux(mu0 = object@mu0, mu = object@mu, sig = object@sigma, 
                         a = object@a, alpha = object@alpha, lambda_p = object@lambda_p,
                         lambda_m = object@lambda_m, S0, Strike = Strike, TimeToMat =TimeToMat, 
                         yield = yield, basis = basis, ret = ret, QMeasure = QMeasure, type = type,
                         Parametrization=object@Parametrization)
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
                   QMeasure, type, Fourier=TRUE, Parametrization){
  mu0Q <- mu0
  if(QMeasure==TRUE){
    # We choose mu0Q in order to ensure the martingale condition
    term1 <- MixedTS:::charact.MTSgam(t= -1i*1,
                             mu0 = 0, mu = mu, sig = sig, a = a, 
                             alpha = alpha, lambda_p = lambda_p,
                             lambda_m = lambda_m, Parametrization=Parametrization)
    if(is.list(ret)){
      mu0Q <- lapply(ret,"-",Re(log(as.numeric(term1))))# ret - log(as.numeric(term1))
    }else{
      mu0Q <- ret - Re(log(as.numeric(term1)))
    }
  }

  
  # we have identify the model with mu0Q, mu, sig, a, alpha, lambda_p and lambda_m
  
  # INSERT FOURIER TRANSFORM
  if(Fourier==FALSE){
    dtvar <- 0.01
    tvar <- seq(-5,5, by = dtvar)
  #   logret <- tvar/(1-tvar^2)
  #   correct <- (1 + tvar^2)/(1 - tvar^2)^2
    # remove inf
    logret <- tvar
    cond <- is.finite(logret)
    
    LogRet <- logret[cond]
    #Correct <- correct[cond]
    Tvar <- tvar[cond]
  #   dens <- MixedTS:::dMixedTS.aux(x = LogRet, mu0 = mu0Q*TimeToMat/basis, mu = mu, sig = sig, a = a*TimeToMat/basis,
  #           alpha = alpha, lambda_p = lambda_p, lambda_m = lambda_m)
    
  price=rep(0,length(S0))
    
  for(i in c(1:length(S0))){  
    Param <- setMixedTS.param(mu0 = mu0Q*TimeToMat[i]/basis[1], mu = mu, sigma = sig,
                              a = a*TimeToMat[i]/basis[1],
                              alpha = alpha, lambda_p = lambda_p,
                              lambda_m = lambda_m, Parametrization=Parametrization)
    
    dum <- dMixedTS(object=Param,x=LogRet)
    
    dens <- dum@dens
    if(type[i]=="CALL"){
  #    payoff <- exp(-ret*TimeToMat/basis)*pmax(S0*exp(LogRet)-Strike,0)*Correct
      payoff <- exp(-ret[i]*TimeToMat[i]/basis[1])*pmax(S0[i]*exp(LogRet)-Strike[i],0)
    }else{
  #   payoff <- exp(-ret*TimeToMat/basis)*pmax(Strike-S0*exp(LogRet),0)*Correct
      payoff <- exp(-ret[i]*TimeToMat[i]/basis[1])*pmax(Strike[i]-S0[i]*exp(LogRet),0)
    }
    integr <- payoff*dens 
    price[i] <- Re(sum(integr*dtvar))
    }
  }else{
    # INSERT HERE Fourier Transform
    envMixedTS <- new.env()
    assign("mu0Q", mu0Q[[1]][[1]], envir = envMixedTS)
    assign("mu", mu, envir = envMixedTS)
    assign("sig", sig, envir = envMixedTS)
    assign("a", a, envir = envMixedTS)
    assign("alpha", alpha, envir = envMixedTS) 
    assign("lambda_p", lambda_p, envir = envMixedTS)
    assign("lambda_m", lambda_m, envir = envMixedTS)
    assign("Parametrization", Parametrization, envir = envMixedTS)
    
    phiMixedTS <- function(u, env){ 
      MixedTS:::charact.MTSgam(t= u, mu0 = env$mu0Q, mu = env$mu, sig = env$sig, 
                               a = env$a, alpha = env$alpha, lambda_p = env$lambda_p,
                               lambda_m = env$lambda_m, Parametrization=env$Parametrization)
    }
    if(is.list(S0)){
      numblist <- length(S0)
      price <-list()
      for(i in c(1:numblist)){
        price[[i]] <- rep(NA,length(S0[[i]]))   
        pCall <- aux.FFTpriceCall(phi=phiMixedTS, env = envMixedTS, 
                                  S0[[i]], K=Strike[[i]], rate=unique(ret[[i]]), 
                                  Time=unique(TimeToMat[[i]]/basis[[i]]), 
                                  alp = 1, N=2^12, eta=0.25)
        condType <- (type[[i]]==rep("CALL", length( type[[i]][[1]] )))
        price[[i]][condType] <- pCall[condType]
        price[[i]][condType==FALSE] <- Strike[[i]][condType]*exp(-ret[[i]][condType]*
                                                                   TimeToMat[[i]][condType]/basis[[i]])- S0[[i]][condType]+pCall[[i]][condType]
      }
    }else{
          price <- rep(0,length(S0))
          
          pCall<-price
          for(i in c(1:length(S0))){
            pCall[i] <- aux.FFTpriceCall(phi=phiMixedTS, env = envMixedTS, S0[i], K=Strike[i], 
                       rate=ret[i], Time=TimeToMat[i]/basis[1], 
                       alp = 1, N=2^12, eta=0.25)
            if(type[i]=="CALL")
              price[i] <- pCall[i]
            else
              price[i] <- Strike[i]*exp(-ret*TimeToMat[i]/basis[i])- S0[i]+pCall[i]
          }
    }
    
  }
  #   length.int <- length(integr)
#   integr.RX <- na.omit(integr[seq(1,length.int, by=3)+2])
#   integr.mid <- na.omit(integr[seq(1,length.int, by=3)+1])
#   integr.LX <- na.omit(integr[seq(1,length.int, by=3)])
#   
#   new.int <- integr[c(1:length(integr.RX))]
#   
#   elem.price <- dtvar*2/6*(integr.LX+integr.RX+4*integr.mid) 
#   price <- sum(elem.price)
  return(price)
}

aux.FFTpriceCall<-function(phi, env, S0, K, rate, Time, alp = 1, N=2^12, eta=0.25){
  # phi is the charateristic function.
  # The code select an EMM with the mean correcting approach.
 # r <- rate
#  m <- rate - log(phi(-(0+1i))) # correcting martingale
#  phi.tilde <- function(u,env){(phi(u,env)*exp(1i*u*m))^(Time)} #eqn pag 321 Iacus Book 
  phi.tilde <- function(u,env){(phi(u,env))^(Time)} #eqn pag 321 Iacus Book 
  # corrected Characteristic function
  psi <-function(v,env){
    dum.den <- alp^2 + alp - v^2 + 1i*(2*alp+1)*v
    res<-exp(-rate*Time)*phi.tilde(v-(alp+1)*1i,env)/(dum.den)
    return(res)
  } #////eqn pag. 319 Iacus Book
  # Characteristic Function with Dumped parameter according Carr Madan Cheng paper 
  # We use the simpson role to increase the precision of the integral in the pricing 
  # formula see page 320 Iacus Book
  lamb <- (2*pi)/(eta*N)
  b <- 0.5*N*lamb
  u_1 <- 0:(N - 1)
  k_u <- -b + lamb * u_1
  v_j <- eta * u_1
  arg_ft <- exp(1i * b * v_j) * psi(v=v_j, env=env) * eta/3 * ( 3 + (-1)^(1:N) - (u_1 == 0) ) 
  res.term <-  exp(-alp * k_u) * fft(arg_ft) / pi
  interp <- approx(x = k_u, y = Re(res.term), xout = log(K/S0))
  res<-  S0 * interp$y 
  return(res)
}



setGeneric("GenerOptMrk",
           function(object, S0, Strike, TimeToMat, ret = 0, 
                    yield = 0, basis = 360, QMeasure = TRUE, 
                    type = "CALL", error = "rnorm(1,mean = 0, sd =0.0005)")
             standardGeneric("GenerOptMrk")
)

setMethod("GenerOptMrk","param.MixedTS",
          function(object, S0, Strike, TimeToMat, ret = 0, 
                   yield = 0, basis = 360,
                   QMeasure = TRUE, type = "CALL",
                   error = "rnorm(1,mean = 0, sd =0.0005)"){
            
            res <- aux.GenerOptMrk(object, S0, Strike, TimeToMat, ret = ret, 
                                   yield = yield, basis = basis,
                                   QMeasure = QMeasure, type = type,
                                   error = error)
            return(res)
          }
          )
aux.GenerOptMrk<-function(object, S0, Strike, TimeToMat, ret = ret, 
                yield = yield, basis = basis,
                QMeasure = QMeasure, type = type,
                error = error){
  Modpar <- object
  Call <- matrix(NA, length(TimeToMat), length(Strike))
  StrikeVect <- TimeToMatVect <- PriceOptVect <- numeric()
  YieldVect <- PriceVect <- rateVect <- numeric()
  for(t in c(1 : length(TimeToMat))){
    for(k in c(1:length(Strike))){
      PriceOptVect[t+(k-1)*dim(Call)[1]] <- OptionPrice(Modpar,S0 = S0, Strike[k],
                                                        TimeToMat[t], ret = ret, yield = yield, basis = basis, 
                                                        QMeasure = QMeasure) + eval(parse(text=error))
      TimeToMatVect[t+(k-1)*dim(Call)[1]] <- TimeToMat[t]
      StrikeVect[t+(k-1)*dim(Call)[1]] <- Strike[k]
      rateVect[t+(k-1)*dim(Call)[1]] <- ret
      PriceVect[t+(k-1)*dim(Call)[1]] <- S0
      YieldVect[t+(k-1)*dim(Call)[1]] <- yield        
    }
  }
  DatasetOpt <- DataOpt(UnderPrice = PriceVect[1:(dim(Call)[1]*dim(Call)[2])], 
                        PriceOpt = PriceOptVect[1:(dim(Call)[1]*dim(Call)[2])], 
                        TimeToMat = TimeToMatVect[1:(dim(Call)[1]*dim(Call)[2])], 
                        Strike = StrikeVect[1:(dim(Call)[1]*dim(Call)[2])], 
                        rate = rateVect[1:(dim(Call)[1]*dim(Call)[2])],
                        qyield = YieldVect[1:(dim(Call)[1]*dim(Call)[2])], 
                        YearBasis = basis, type = type)       
  return(DatasetOpt)
}