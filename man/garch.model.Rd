\name{GARCH.model}
\alias{GARCH.model}
\title{
  Definition of an GARCH model.
}
\description{
  \command{GARCH.model} defines an GARCH model as input for the
  argument \code{user} of \command{\link{tsgui}}
}
\usage{
GARCH.model(a.max=3, b.max=3)	  
}
\arguments{
  \item{a.max,b.max}{non-negative integer.
    Any GARCH(p,q) model can be modelled in the gui
    with \eqn{p \le}{p <=}\code{a.max} and  \eqn{q \le}{q <=}\code{b.max}.
    Note that \code{a.max} + \code{b.max} should not be larger than
    6. Otherwise gui cannot manage the abundance of parameters.
  }
}

\value{
  A list of parameters is returned.
}
%\references{}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{\command{\link{tsgui}}, \command{\link{ARMA.model}}}
  
\examples{  
if (interactive())
  tsgui(user=list(ARCH=GARCH.model(a.max=5, b.max=0)),
        included.models = "GARCH")
}

\keyword{ts}
\keyword{models}
