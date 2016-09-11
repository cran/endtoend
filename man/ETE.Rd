\name{ETE}
\alias{ETE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical transmissions/receptions for a L-limited End to End model
}
\description{
This function computes the expected value of the number of transmissions/receptions for End to End model with L-limited retransmissions per packet.
}
\usage{
ETE(p1, p2, L, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
Data success probability
}
  \item{p2}{
ACK success probability
}
  \item{L}{
Maximum number of retransmissions.
}
  \item{N}{
Number of Hops
}
}

\details{
When there is no limitation, \code{L} value must be set as \code{L=Inf}.
}
\value{
The ouput is a matrix containing the following values:

\item{1}{Success Probability}
\item{2}{Expected Data Transmissions}
\item{3}{Expected ACK Transmissions}
\item{4}{Expected Total Transmissions}
\item{5}{Expected Data Receptions}
\item{6}{Expected ACK Receptions}
\item{7}{Expected Total Receptions}
}
\references{
Heimlicher, S., Nuggehalli, P., & May, M. (2007). End-to-end vs. hop-by-hop transport. ACM SIGMETRICS Performance Evaluation Review, 35(3), 59.

Heimlicher, S., Karaliopoulos, M., Levy, H., & May, M. (2007). End-to-end vs. Hop-by-hop Transport under Intermittent Connectivity (Invited Paper). Proceedings of the First International Conference on Autonomic Computing and Communication Systems.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{MCETE}}
}
\examples{
#An N=5 End to End system with limited L=7 retransmission per hop
ETE(p1=0.65,p2=0.4,L=7,N=5)

#An unlimited N=5 End to End system
ETE(p1=0.65,p2=0.4,L=Inf,N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{End to End}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
