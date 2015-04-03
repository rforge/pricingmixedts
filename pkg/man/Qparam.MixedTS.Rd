\name{Qparam.MixedTS}
\alias{Qparam.MixedTS}
\alias{Qparam.MixedTS,generic}
\alias{Qparam.MixedTS,param.MixedTS-method}
\alias{EsscherMartingaleMeasure}
\alias{MartingaleMeasure}
\alias{PricingMeasure}
\alias{IncompleteMarket}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{\code{Qparam.MixedTS}: MixedTS parameters under Esscher Martingale Measure.
}
\description{The \code{Qparam.MixedTS} converts the parameters of an exponential \code{MixedTS} process under the real measure into the corresponding parameters under the equivalent martingale measure obtained by means of Esscher Transform.}
\usage{
\S4method{Qparam.MixedTS}{param.MixedTS}(object, ret)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{param.MixedTS-class} that contains the model parameters under real measure.}
  \item{ret}{a numeric object that contains the costant risk free rate.}
}
\details{
The exponential \code{MixedTS} Levy process is an incomplete market, then  equivalent martingale measure (EMM) for pricing is not unique. The method \code{Qparam.MixedTS} uses the Esscher Transform that leaves invariant the Levy measure under the EMM.  
}
\value{The method returns an object of class \code{param.MixedTS-class}.
}
\references{Available as soon as possible}
\author{
PricingMixedTS team.
}

\examples{
## Available as soon as possible

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
