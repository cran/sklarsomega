% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sklarsomega.R
\name{residuals.sklarsomega}
\alias{residuals.sklarsomega}
\title{Extract model residuals.}
\usage{
\method{residuals}{sklarsomega}(object, ...)
}
\arguments{
\item{object}{an object of class \code{sklarsomega}, typically the result of a call to \code{\link{sklars.omega}}.}

\item{\dots}{additional arguments.}
}
\value{
A vector of residuals.
}
\description{
Extract model residuals.
}
\details{
Although residuals may not be terribly useful in this context, we provide residuals nonetheless. Said residuals are computed by first applying the probability integral transform, then applying the inverse probability integral transform, then pre-multiplying by the inverse of the square root of the (fitted) copula correlation matrix. For nominal or ordinal scores, the distributional transform approximation is used.
}
\examples{
# Fit a subset of the cartilage data, assuming a Laplace marginal distribution.
# Produce a normal probability plot of the residuals, and overlay the line y = x.

data(cartilage)
data.cart = as.matrix(cartilage)[1:100, ]
colnames(data.cart) = c("c.1.1", "c.2.1")
fit.lap = sklars.omega(data.cart, level = "balance", control = list(dist = "laplace"))
summary(fit.lap)
res = residuals(fit.lap)
qqnorm(res, pch = 20)
abline(0, 1, col = "blue", lwd = 2)
}
\references{
Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
}
\seealso{
\code{\link{sklars.omega}}
}
