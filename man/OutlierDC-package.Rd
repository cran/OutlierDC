\name{OutlierDC-package}
\alias{OutlierDC}
\alias{OutlierDC-package}
\docType{package}
\title{
	Functions for detecting outlying observations for censored data using quantile regression
}
\description{
  This package offers three outlier detection algorithms for censored data using quantile regression. 
}
\details{
\tabular{ll}{
Package: \tab OutlierDC\cr
Type: \tab Package\cr
Version: \tab 0.3-0\cr
Date: \tab 2014-03-23\cr
License: \tab GPL version 3\cr
LazyLoad: \tab no\cr
}
}
\author{
Soo-Heang Eo and HyungJun Cho 
Maintainer: Soo-Heang Eo <eo.sooheang@gmail.com>
}
\references{
Eo S-H, Hong S-M Hong, Cho H (2014). Identification of outlying observations with quantile regression for censored data, \emph{Submitted}.

Wang HJ, Wang L (2009) Locally Weighted Censored Quantile Regression. JASA 104:1117--1128. doi: 10.1198/jasa.2009.tm08230
}
\note{
  We would like to thank Huxia Judy Wang and Lan Wang for permission to use their LCRQ functions.
}
\seealso{
    \code{\link{odc}}, \code{\link{plot}}, \code{\link{coef}}, \code{\link{show}}, \code{\link[quantreg:crq]{quantreg}}
}
\examples{
  \dontrun{
    require(OutlierDC)
    # Toy example 
    data(ebd)
    # The data consists of 402 observations with 6 variables. 
    dim(ebd)
    # To show the first six observations of the dataset,
    head(ebd)
    
    #scoring algorithm
    fit <- odc(Surv(log(time), status) ~ meta, data = ebd)
    fit
    coef(fit)
    plot(fit)

    # Add upper bound for the selection of outleirs
    fit1 <- update(fit, k_s = 4)
    fit1
    plot(fit1)

    # residual-based algorithm
    fit2 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "residual", k_r = 1.5)
    fit2
    plot(fit2)
    
    # To display all of outlying observations in the fitted object
    fit2@outlier.data
    
    # boxplot algorithm
    fit3 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "boxplot", k_b = 1.5)
    fit3
    plot(fit3, ylab = "log survival times", xlab = "metastasis lymph nodes")
}
}
\keyword{ package }