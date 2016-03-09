\name{OptionPrice}
\alias{OptionPrice}
\alias{OptionPrice,generic}
\alias{OptionPrice,param.MixedTS-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{\code{OptionPrice}: Method for Evaluation of a CALL/PUT Option.
}
\description{The \code{OptionPrice} computes the no-arbitrage price for a CALL or PUT option}
\usage{
\S4method{OptionPrice}{param.MixedTS}(object, S0, Strike, TimeToMat, ret = 0, yield = 0,
  basis = 360, QMeasure = FALSE,  type = "CALL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{param.MixedTS-class} that contains the model parameters under real measure.}
  \item{S0}{Current price of the underlying asset.}
  \item{Strike}{Level of the Strike price.}
  \item{TimeToMat}{Time to Maturity.}
  \item{ret}{a numeric object that contains the costant risk free rate.}
  \item{yield}{convenience yield of the underlying asset.}
  \item{basis}{Identify how to counter the time interval.}
  \item{QMeasure}{a logical variable. If \code{Qmeasure = FALSE} the equivalent martingale
is selected according to the Esscher Transform otherwise th mean correcting
martingale measure is used when \code{Qmeasure = TRUE.}.}
  \item{type}{Identify the nature of the Option. Default \code{CALL}.}
}
%\details{}
\value{The method returns a object containing option prices.}

\references{Available as soon as possible}
\author{
PricingMixedTS team.
}

\examples{
## Available as soon as possible

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
