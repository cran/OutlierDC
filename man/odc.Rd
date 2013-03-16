\name{odc}
\alias{odc}
\title{Outlier detection using quantile regression for censored data}
\description{
	outlier detection algorithms using quantile regression for censored data 
}
\usage{
     odc(formula, data, method = c("score", "boxplot","residual"), 
      rq.model = c("Wang", "PengHuang", "Portnoy"), k = 1.5, h = .05 )
}

\arguments{
   \item{formula}{a type of \code{Formula} object with a \code{survival} object on the left-hand side of the ~ operator and covariate terms on the right-hand side. The survival object with survival time and its censoring status is constructed by the \code{\link[survival]{Surv}} function in \code{survival} package.}
   \item{data}{a data frame with variables used in the \code{formula}. It needs at least three variables, including survival time, censoring status, and covariates. }
  \item{method}{the outlier detection method to be used. The options \code{"socre", "boxplot"}, and \code{"residual"} conduct the scoring, boxplot, and residual-based algorithm, respectively. The default algorithm is \code{"score"}.
  }
  \item{rq.model}{the type of censored quantile regression to be used for fitting. The options \code{"Wang", "Portnoy"}, and \code{"PengHuang"} conduct Wang and Wang's, Portnoy's, and Peng and Huang's censored quantile regression approaches, respectively. The default is \code{"Wang"}.
  }
  \item{k}{a constant to control the tightness of cut-offs for residual-based and boxplot algorithms with a default value of 1.5.
   }
   \item{h}{bandwidth for locally weighted censored quantile regression with a default value of 0.05.
   }
}
\details{
	The \code{odc} function conducts three outlier detection algorithms on the basis of censored quantile regression. Three outlier detection algorithms were implemented: residual-based, boxplot, and scoring algorithms. The residual-based algorithm detects outlying observations using constant scale estimates; however, it does not account for the heterogeneity of variablity. When the data is extremely heterogeneous, the boxplot algorithm with censored quantile regression is more effective. The residual-based and boxplot algorithms produce cut-offs to determine whether observations are outliers. In contrast, the scoring algorithm provides the outlying magnitude or deviation of each point from the distribution of observations. Outlier detection is achieved by visualizing the scores. 
}
\value{
	an object of the S4 class "OutlierDC" with the following slots:\cr
	call: evaluated function call\cr
  formula: formula to be used\cr
  raw.data: data to be used for model fitting\cr
  refined.data: the data set after removing outliers\cr
  coefficients: the estimated censored quantile regression coefficient matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr
  fitted.mat: the censored quantile regression fitted value matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr 
  score: outlying scores (scoring algorithm) or residuals (residual-based algorithm)\cr
  cutoff: estimated scale parameter for the residual-based algoritm\cr
  lower: lower fence vector used for the boxplot and scoring algorithms\cr
  upper: upper fence vector used for the boxplot and scoring algorithms\cr
  outliers: logical vector to determine which observations are outliers\cr
  n.outliers: number of outliers detected\cr
  method: outlier detection method to be used\cr
  rq.model: censored quantile regression to be used\cr
  k: constant value to be used for the tightness of cut-offs in the residual-based and boxplot algorithms\cr
  UB: upper fence used for the scoring algorithm with the \code{\link{update}} function\cr
  LB: lower fence used for the scoring algorithm with the \code{\link{update}} function\cr
}
\author{
   Soo-Heang Eo <hanansh@korea.ac.kr> \cr
   HyungJun Cho <hj4cho@korea.ac.kr>
}
\references{
Eo SH, Hong SM and Cho H (2013+). OutlierDC: Outlier Detection Algorithms using Quantile Regression for Censored Survival Data. \emph{submitted}.
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}
}

\examples{
  \dontrun{
    require(OutlierDC)
    data(ebd)
    
    #scoring algorithm
    fit <- odc(Surv(log(time), status) ~ meta, data = ebd)
    fit
    coef(fit)
    plot(fit)
    fit1 <- update(fit, UB = 4, LB = -4)
    fit1
    plot(fit1)

    fit2 <- update(fit, UB = NA, LB = -4)
    fit3 <- update(fit, UB = 4, LB = NA)
 
    # residual-based algorithm
    fit4 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "residual", k = 3)
    fit4
    plot(fit4)
    
    # boxplot algorithm
    fit5 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "boxplot", k = 1.5)
    fit5
    plot(fit5, ylab = "log survival times", xlab = "metastasis lymph nodes")
}
}
\keyword{models}


