\name{estTau}
\alias{estTau}
\title{Estimate Tau}

\description{This function estimates tau for a given set of graphs.}

\usage{estTau(data, type, gstar)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{gstar}{A single columned data frame to be used as the g-star of the data set.}
}

\value{The tau value for the data based on g.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	braintau
}
