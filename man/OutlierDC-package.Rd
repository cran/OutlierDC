\name{OutlierDC-package}
\alias{OutlierDC-package}
\docType{package}
\title{
	Functions for Dectecting Outlying Observations using Quantile Regression in Censored Data
}
\description{
  This package offers three outlier detection algorithms based on censored quantile regression. 
}
\details{
\tabular{ll}{
Package: \tab OutlierDC\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2013-03-16\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab no\cr
}
The package OutlierDC is the implementation of the outlier detection algorithms proposed by Eo, Hong and Cho (2013+). Three outlier detection algorithms using quantile regression are proposed, which are modified versions of existing algorithms for uncensored data. We illustrated with real data from the SEER database, which contains a number of data sets related to various cancers.
}
\author{
Soo-Heang Eo Seung-Mo Hong and HyungJun Cho 
Maintainer: Soo-Heang Eo <hanansh@korea.ac.kr>
}
\references{
Eo SH, Hong SM and Cho H (2013+). OutlierDC: Outlier Detection Algorithms using Quantile Regression for Censored Survival Data. \emph{submitted}.
}
\seealso{
    \code{\link{odc}}, \code{\link{plot}}, \code{\link{coef}}, \code{\link{show}}, \code{\link[quantreg:crq]{quantreg}}
}
\keyword{ package }