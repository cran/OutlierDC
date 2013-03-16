# MODIFIED ON   :   Monday, February 23, 2009 at 11:23
# AUTHOR        :   HUIXIA JUDY WANG and LAN WANG
# AFFILIATION   :   DEPARTMENT OF STATISTICS, NORTH CAROLINA STATE UNIVERSITY AND DEPARTMENT OF STATISTICS, UNIVERSITY OF MINNESOTA
# EMAIL         :   WANG@STAT.NCSU.EDU 
# FUNCTION      :   LOCALLY WEIGHTED CENSORED QUANTILE REGRESSION ESTIMATOR

###################################################################################
## tauhat.func: Generalized Kaplan-Meier estimator (Gonzalez-Manteiga and Cadarso-Suarez, 1994) 
## Note: this is equiavlent to Beran's estimator (1981)
## Z=min(Y,C), delta=I(Y<=C)
## we are estimating P(Y<=y0|x=x0)
####################################################################################
Bnk.func = function(x0, x, h)
{
    # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
    # returns a vector
    # h is the bandwidth
    xx<-(x-x0)/h  
    xx[abs(xx)>=1]<-1
    w<-15*(1-xx^2)^2/16  #biquadratic kernel 
    w<-w/sum(w)
    return(w)
}
tauhat.func <- function(y0, x0, z, x, delta,h)
{
# tau0(y0, x0) = F(T<y0|x0); so y0 is the C_i, and x0 is the xi in the paper
# z is observed vector of response variable
# x is the observed covariate
# delta is the censoring indicator function
# h is the bandwidth
    n<-length(z)
    ###kernel weights#########################
    Bn = Bnk.func(x0, x, h)
    if (y0<max(z))
    {
        # sort the data z, and the delta, Bn correspondingly to the order of sorted y
        z2 = sort(z)
        Order = order(z) # so z[Order] = z2
        Bn2 = Bn[Order]
        delta2 = delta[Order]
        eta = which(delta2==1 & z2<=y0) # the index of those observations satisfying delta2==1 & z2<=y0
        Bn3 = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
        tmp = 1- Bn2 /cumsum(Bn3)[n:1]  
        out = 1-prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
    } 
 else out<-1 
    return(out)
}

LCRQ <- function(y, x, delta, tau, h)
{
# x is one dimensional, not including the intercept yet
# y is the observed survival time = min(T, C)
# delta is the censoring indicator function with 1 standing for uncensored, and 0 censored
# tau is the quantile level of interest
# h is the handwidth used in calculating tauhat for each censored observation
    n = length(y)
    ind = which(delta==0)
    w = rep(1, n)  # the weight
    if(length(ind)>=1)
    {
        for(i in 1:length(ind))
        {
            x0 = x[ind[i]]
            y0 = y[ind[i]]
            tau.star = tauhat.func(y0,x0, y, x, delta,h=h)
            if (tau>tau.star) w[ind[i]] = (tau-tau.star)/(1-tau.star)
        }
    # pseudo observations
    ind2 = which(w!=1)
    y.pse = rep(max(y)+100, length(ind2))
    x.pse = x[ind2]
    yy = c(y, y.pse)
    xx = c(x, x.pse)
    ww = c(w, 1-w[ind2])
    }
    else
    {
        yy=y; xx=x; ww=w
    }
    rq1 = rq(yy~xx, weights=ww, tau=tau)
    result<- rq1$coeff
    result    
}


# a simulated data set (generated under the global linearity assumption) for demonstration
#library(quantreg)
#set.seed(123456)
#tau=0.5; n=200
#x=runif(n,0,1)
#b0=3
#b1=5
#ystar = b0 + x*b1 + rnorm(n)-qnorm(tau)
#ct=runif(n,0,26)
#y = pmin(ystar, ct)
#delta = 1*(y==ystar)
#LCRQ(y, x, delta, tau=tau, h=0.05)
#(naive = rq(y~x, tau=tau)$coef)
#fit <- crq(Surv(y, delta)~x, method = "Por")
#summary(fit, c(tau, tau))[[1]]$coef
