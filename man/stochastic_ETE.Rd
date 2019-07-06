\name{stochastic_ETE}
\alias{stochastic_ETE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Probabilities Monte Carlo transmissions/receptions simulations for a L-limited End to End model
}
\description{
This function compute the mean of the number of transmissions/receptions for End to End model with L-limited retransmissions per packet simulating via Monte Carlo.
}
\usage{
stochastic_ETE(dist1,p11,p12,dist2,p21,p22,L,N,M=10^5,printout=TRUE,plotspdf=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dist1}{
For the data success probability: probability density function. Options are "uniform" and "beta".
}
  \item{p11}{
For the data success probability: lower limit of the uniform distribution (\code{dist1 == "uniform"}) or shape1 (alpha) paremeter of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{p12}{
For the data success probability: upper limit of the uniform distribution (\code{dist1 == "uniform"}) or shape2 (beta) paremeter of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{dist2}{
For the ACK success probability: probability density function. Options are "uniform" and "beta".
}
  \item{p21}{
For the ACK success probability: lower limit of the uniform distribution (\code{dist1 == "uniform"}) or shape1 (alpha) paremeter of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{p22}{
For the ACK success probability: upper limit of the uniform distribution (\code{dist1 == "uniform"}) or shape2 (beta) paremeter of a Beta distribution (\code{dist1 == "beta"}).
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
  \item{printout}{
If \code{TRUE} (by default), the function prints some outputs and plots
}
  \item{plotspdf}{
If \code{TRUE} (by default), the function exports all plots in pdf in the working directory
}
}
\value{
The ouput is a matrix containing two elements:

\item{data}{a dataframe containing all Monte Carlo replications}
\item{stats}{descriptive statistics}
for
\item{1}{p1}
\item{2}{p2}
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
\code{\link{ETE}},\code{\link{MCETE}}
}
\examples{
#Monte Carlo simulations for an N=5 End to End system
#with limited L=7 retransmission per hop

#We now consider p1 ~ Uniform(0.2,0.6)
dist1 = "uniform"
p11 = 0.2
p12 = 0.6

#and p2 ~ Beta(3,1)
dist2 = "beta"
p21 = 3
p22 = 1

#no outputs and plots
out = stochastic_ETE(dist1,p11,p12,dist2,p21,p22,L=7,N=5,M=5*10^3,printout=FALSE,plotspdf=FALSE)
out$data  #simulations
out$stats #resume

#uncommnet next line for outputs plots and pdf file
#out = stochastic_ETE(dist1,p11,p12,dist2,p21,p22,L=7,N=5,M=5*10^3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{End to End}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
