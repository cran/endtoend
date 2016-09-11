\name{MCETE}
\alias{MCETE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo transmissions/receptions simulations for a L-limited End to End model
}
\description{
This function compute the mean of the number of transmissions/receptions for End to End model with L-limited retransmissions per packet simulating via Monte Carlo.
}
\usage{
MCETE(p1, p2, L, N, M = 5000)
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
Maximum number of retransmissions
}
  \item{N}{
Number of Hops
}
  \item{M}{
Number of Monte Carlo Simulations
}
}
\value{
The ouput is a matrix containing the following values:

\item{1}{MC Success Probability}
\item{2}{MC Mean Data Transmissions}
\item{3}{MC Mean ACK Transmissions}
\item{4}{MC Mean Total Transmissions}
\item{5}{MC Mean Data Receptions}
\item{6}{MC Mean ACK Receptions}
\item{7}{MC Mean Total Receptions}
}
\references{
Heimlicher, S., Nuggehalli, P., & May, M. (2007). End-to-end vs. hop-by-hop transport. ACM SIGMETRICS Performance Evaluation Review, 35(3), 59.

Heimlicher, S., Karaliopoulos, M., Levy, H., & May, M. (2007). End-to-end vs. Hop-by-hop Transport under Intermittent Connectivity (Invited Paper). Proceedings of the First International Conference on Autonomic Computing and Communication Systems.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{ETE}}
}
\examples{
#Monte Carlo simulations for an N=5 End to End system
#with limited L=7 retransmission per hop

MCETE(p1=0.65,p2=0.4,L=7,N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{End to End}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
