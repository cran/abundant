\name{pred.response}
\alias{pred.response}
\title{
Predict the response with the fitted high-dimensional principal fitted components model
}
\description{
Let \eqn{x\in R^p} denote the values of the \eqn{p} predictors.
This function computes \eqn{\widehat E(Y|X=x)} using equation (8.1)
of Cook, Forzani, and Rothman (2012).
}
\usage{
pred.response(fit, newx=NULL)
}

\arguments{
  \item{fit}{
The object returned by \code{fit.pfc()}. 
}
  \item{newx}{
A matrix with \eqn{N} rows and \eqn{p} columns where each row is an instance of \eqn{x} described above.
If this argument is unspecified, then the fitted values are returned, i.e, 
\code{newx=X}, where \code{X} was the predictor matrix used in
the call to \code{fit.pfc()}.  
}
 
}
\details{
See Cook, Forzani, and Rothman (2012) for more information.
}
\value{
A vector of response prediction with \code{nrow(newx)} entries. 
}
\references{
Cook, R. D., Forzani, L., and Rothman, A. J. (2012).  
Estimating sufficient reductions of the predictors in abundant high-dimensional regressions. 
Annals of Statistics 40(1), 353-384.
}
\author{
Adam J. Rothman
}
 
\seealso{
\code{\link{fit.pfc}}
}

\examples{
set.seed(1)
n=25
p=50
d=1
true.G = matrix(rnorm(p*d), nrow=p, ncol=d)
y=rnorm(n)
fy = y
E=matrix(rnorm(n*p), nrow=n, ncol=p) 
X=fy\%*\%t(true.G) + E
fit=fit.pfc(X=X, r=4, d=d, y=y, weight.type="diag")
fitted.values=pred.response(fit)
mean((y-fitted.values)^2)
plot(fitted.values, y)

n.new=100
y.new=rnorm(n.new)
fy.new=y.new
E.new=matrix(rnorm(n.new*p), nrow=n.new, ncol=p) 
X.new = fy.new\%*\%t(true.G) + E.new
mean((y.new - pred.response(fit, newx=X.new))^2)  

} 
