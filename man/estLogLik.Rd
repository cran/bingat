\name{estLogLik}
\alias{estLogLik}
\title{Estimate the Log Likelihood Value}

\description{This function estimates log likelihood value for a given graph.}

\usage{estLogLik(data, type, g, tau)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{g}{A single columned data frame to estimate the likelihood for.}
  \item{tau}{A single value used in estimating the likelihood.}
}

\value{The likelihood value of the input graph g.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	brainll <- estLogLik(braingraphs, "adjMatrix", braingstar, braintau)
	brainll
}
