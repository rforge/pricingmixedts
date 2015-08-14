\name{DataOpt}
\alias{DataOpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Options Information available on the markets.}
\description{The constructor \code{DataOpt} builds an object of class \code{OptionData-class}. This object contains option prices and volatilities available on the markets.}
\usage{
DataOpt(UnderPrice, PriceOpt, ImpliedVol, TimeToMat, Strike,
  type = "CALL", rate = 0, qyield = 0, Date = Sys.Date(), YearBasis = 360)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UnderPrice}{ a numeric object that contains the actual price of the underlying asset.}
  \item{PriceOpt}{ a numeric object that contains the price of  call or put options.}
  \item{ImpliedVol}{ a numeric object that contains the implied volatility of  call or put options.}
  \item{TimeToMat}{ a numeric object that contains the time to maturities of the call/put options.}
  \item{Strike}{ a numeric object that contains the strike prices.}
  \item{type}{ a character object that indicates the nature of the available options. Default \code{call}.}
  \item{rate}{ a numeric object that is the risk-free rate.}
  \item{qyield}{ a numeric object that is the dividend yield rate. }
  \item{Date}{ a Date object that indicates the date that the observation occurs.}
  \item{YearBasis}{ Day counter basis.}
}
\value{ The function returns an object of \code{OptionData-class}.}
\references{ Available as soon as possible}
\author{PricingMixedTS Team}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.