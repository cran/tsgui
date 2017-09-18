\name{ARMA.model}
\alias{ARMA.model}
\title{
  Definition of an ARMA model.
}
\description{
  \command{ARMA.model} defines an ARMA model as input for the
  argument \code{user} of \command{\link{tsgui}}
}
\usage{
ARMA.model(p.max=3, q.max=3)	  
}
\arguments{
  \item{p.max,q.max}{non-negative integer.
    Any ARMA(p,q) model can be modelled in the gui
    with \eqn{p \le}{p <=}\code{p.max} and  \eqn{q \le}{q <=}\code{q.max}.
    Note that \code{p.max} + \code{q.max} should not be larger than
    6. Otherwise gui cannot manage the abundance of parameters.
  }
}

\value{
  A list of parameters is returned.
}
%\references{}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{\command{\link{tsgui}}, \command{\link{GARCH.model}}}
  
\examples{  
  if (interactive())
    tsgui(user=list(AR=ARMA.model(p.max=5, q.max=0)),
          included.models = "ARMA")
}

\keyword{ts}
\keyword{models}
