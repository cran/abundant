\name{fit.pfc}
\alias{fit.pfc}
\title{
Fit a high-dimensional principal fitted components model 
using the method of Cook, Forzani, and Rothman (2012).
}
\description{
Let \eqn{(x_1, y_1), \ldots, (x_n, y_n)} denote the \eqn{n} measurements of
the predictor and response, where \eqn{x_i\in R^p} and \eqn{y_i \in R}.
The model assumes that these measurements 
are a realization
of \eqn{n} independent copies of 
the random vector \eqn{(X,Y)'}, where
\deqn{
X = \mu_X + \Gamma \beta\{f(Y) - \mu_f\}+ \epsilon,
}
\eqn{\mu_X\in R^p}; \eqn{\Gamma\in R^{p\times d}} with rank \eqn{d}; 
\eqn{\beta \in R^{d\times r}} with rank \eqn{d}; \eqn{f: R \rightarrow R^r} is a known
vector valued function; \eqn{\mu_f = E\{f(Y)\}};  
\eqn{\epsilon \sim N_p(0, \Delta)}; and \eqn{Y} is independent of \eqn{\epsilon}.
The central subspace is \eqn{\Delta^{-1} {\rm span}(\Gamma)}.
 
 
This function computes estimates of these model parameters
by imposing constraints for identifiability.
The mean parameters \eqn{\mu_X} and \eqn{\mu_f}
are estimated with \eqn{\bar x = n^{-1}\sum_{i=1}^n x_i} and
\eqn{\bar f = n^{-1} \sum_{i=1}^n f(y_i)}.
Let \eqn{\widehat\Phi = n^{-1}\sum_{i=1}^{n} \{f(y_i) - \bar f\}\{f(y_i) - \bar f\}'},
which we require to be positive definite.
Given a user-specified weight matrix \eqn{\widehat W}, 
let
\deqn{
(\widehat\Gamma, \widehat\beta) = 
\arg\min_{G\in R^{p\times d}, B \in R^{d\times r}}  \sum_{i=1}^n [x_i - \bar x - GB\{f(y_i) - \bar f\}]'\widehat W [x_i - \bar x - GB\{f(y_i) - \bar f\}],
} 
subject to the constraints that \eqn{G'\widehat W G} is diagonal and 
\eqn{B \widehat\Phi B' = I}.  The sufficient reduction estimate \eqn{\widehat R: R^p \rightarrow R^d}
is defined by
\deqn{
\widehat R(x) =  (\widehat\Gamma'\widehat W \widehat\Gamma)^{-1} \widehat\Gamma' \widehat W(x - \bar x).
}

}
\usage{
fit.pfc(X, y, r=4, d=NULL, F.user=NULL, weight.type=c("sample", "diag", "L1"), 
        lam.vec=NULL, kfold=5, silent=TRUE, qrtol=1e-10, cov.tol=1e-4, 
        cov.maxit=1e3, NPERM=1e3, level=0.01)
}

\arguments{
  \item{X}{
The predictor matrix with \eqn{n} rows and \eqn{p} columns. The \eqn{i}th row is \eqn{x_i} defined above.
}
  \item{y}{
The vector of measured responses with \eqn{n} entries.  The \eqn{i}th entry is \eqn{y_i} defined above. 
}
  \item{r}{
When polynomial basis functions are used (which is the case when \code{F.user=NULL}), \code{r} is the polynomial 
order, i.e, 
\eqn{f(y) = (y, y^2, \ldots, y^r)'}.  The default is \code{r=4}.  This argument is not used
when \code{F.user} is specified.
}
  \item{d}{
The dimension of the central subspace defined above.  This must be specified by the user
when \code{weight.type="L1"}.
If unspecified by the user this function will use the sequential permutation testing procedure, described in Section 8.2
of Cook, Forzani, and Rothman (2012), to select \code{d}.  
}
  \item{F.user}{
A matrix with \eqn{n} rows and \eqn{r} columns, where the \eqn{i}th row is \eqn{f(y_i)} defined above. 
This argument is optional, and will typically be used when polynomial basis functions are not desired.
}
\item{weight.type}{
The type of weight matrix estimate \eqn{\widehat W} to use.  
Let \eqn{\widehat\Delta} be the observed residual sample covariance matrix for the multivariate
regression of \code{X} on \eqn{f}\code{(Y)} with \eqn{n-r-1} scaling.
There are three options for \eqn{\widehat W}: 
\itemize{
\item \code{weight.type="sample"} uses a Moore-Penrose generalized inverse of \eqn{\widehat\Delta} for \eqn{\widehat W},
when \eqn{p \leq n-r-1} this becomes the inverse of \eqn{\widehat\Delta};
\item \code{weight.type="diag"} uses the inverse of the diagonal matrix with the same diagonal as \eqn{\widehat\Delta}
for \eqn{\widehat W};
\item \code{weight.type="L1"} uses the L1-penalized inverse of \eqn{\widehat\Delta} described in equation (5.4) of Cook, Forzani,
and Rothman (2012).  In this case, \code{lam.vec} and \code{d} must be specified by the user.
The glasso algorithm of Friedman et al. (2008) is used through the R package \code{glasso}.
}
}
\item{lam.vec}{
A vector of candidate tuning parameter values to use when \code{weight.type="L1"}.  If this vector has more than one entry,
then \code{kfold} cross validation will be performed to select the optimal tuning parameter value.  
}
\item{kfold}{
The number of folds to use in cross-validation to select the optimal tuning parameter when  \code{weight.type="L1"}.
Only used if \code{lam.vec} has more than one entry.
}
  \item{silent}{
Logical. When \code{silent=FALSE}, progress updates are printed.  
}
  \item{qrtol}{
The tolerance for calls to \code{qr.solve()}. 
}
  \item{cov.tol}{
The convergence tolerance for the QUIC algorithm used when \code{weight.type="L1"}. 
}
  \item{cov.maxit}{
The maximum number of iterations allowed for the QUIC algorithm used when \code{weight.type="L1"}. 
}
  \item{NPERM}{
The number of permutations to used in the sequential permutation testing procedure to select \eqn{d}. 
Only used when \code{d} is unspecified.
}
\item{level}{
The significance level to use to terminate the sequential permutation testing procedure to select \eqn{d}. 
}
}
\details{
See Cook, Forzani, and Rothman (2012) more information.
}
\value{
A list with
\item{Gamhat}{this is \eqn{\widehat\Gamma} described above.}
\item{bhat}{this is \eqn{\widehat\beta} described above.}
\item{Rmat}{this is \eqn{\widehat W\widehat\Gamma(\widehat\Gamma'\widehat W \widehat\Gamma)^{-1}}.}
\item{What}{this is \eqn{\widehat W} described above.}
\item{d}{this is \eqn{d} described above.}
\item{r}{this is \eqn{r} described above.}
\item{GWG}{this is \eqn{\widehat\Gamma'\widehat W \widehat\Gamma}}
\item{fc}{a matrix with \eqn{n} rows and \eqn{r} columns where the \eqn{i}th row is \eqn{f(y_i) - \bar f}.}
\item{Xc}{a matrix with \eqn{n} rows and \eqn{p} columns where the \eqn{i}th row is \eqn{x_i - \bar x}.}
\item{y}{the vector of \eqn{n} response measurements.}
\item{mx}{this is \eqn{\bar x} described above.}
\item{mf}{this is \eqn{\bar f} described above.}
\item{best.lam}{this is selected tuning parameter value used when \code{weight.type="L1"}, will be \code{NULL} otherwise.}
\item{lam.vec}{this is the vector of candidate tuning parameter values used when 
\code{weight.type="L1"}, will be \code{NULL} otherwise.}
\item{err.vec}{this is the vector of validation errors from cross validation, one error for each entry in \code{lam.vec}.
Will be \code{NULL} unless \code{weight.type="L1"} and \code{lam.vec} has more than one entry.}
\item{test.info}{a dataframe that summarizes the results from the sequential testing procedure.  Will be \code{NULL}
unless \code{d} is unspecified.
}
}
\references{
Cook, R. D., Forzani, L., and Rothman, A. J. (2012).  
Estimating sufficient reductions of the predictors in abundant high-dimensional regressions. 
Annals of Statistics 40(1), 353-384.

Friedman, J.,  Hastie, T., and Tibshirani R. (2008). 
Sparse inverse covariance estimation with the lasso. 
Biostatistics 9(3), 432-441.
}
\author{
Adam J. Rothman
}
\seealso{
\code{\link{pred.response}}
}

\examples{
set.seed(1)
n=20
p=30
d=2
y=sqrt(12)*runif(n)
Gam=matrix(rnorm(p*d), nrow=p, ncol=d)
beta=diag(2)
E=matrix(0.5*rnorm(n*p), nrow=n, ncol=p)
V=matrix(c(1, sqrt(12), sqrt(12), 12.8), nrow=2, ncol=2)
tmp=eigen(V, symmetric=TRUE)
V.msqrt=tcrossprod(tmp$vec*rep(tmp$val^(-0.5), each=2), tmp$vec)
Fyc=cbind(y-sqrt(3),y^2-4)\%*\%V.msqrt
X=0+Fyc\%*\%t(beta)\%*\%t(Gam) + E

fit=fit.pfc(X=X, y=y, r=3, weight.type="sample")
## display hypothesis testing information for selecting d
fit$test.info
##  make a response versus fitted values plot
plot(pred.response(fit), y)
} 
