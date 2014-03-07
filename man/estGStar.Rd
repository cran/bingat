\name{estGStar}
\alias{estGStar}
\title{Estimate G-Star}

\description{This function estimates the g-star graph for a given set of graphs.}

\usage{estGStar(data, cutoff = 0.5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\value{A single vector that is the gstar is returned.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)

	braingstar <- estGStar(braingraphs) 
	braingstar[1:25]
}