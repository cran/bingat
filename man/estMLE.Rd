\name{estMLE}
\alias{estMLE}
\title{Estimate the MLE Parameters}

\description{This function estimates the MLE parameters g-star and tau for a given set of graphs.}

\usage{estMLE(data, type, cutoff = .5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\details{
	Essentially this function calls both \code{estGStar} and \code{estTau} and returns the results. \cr
}

\value{A list containing g-star and tau named gstar and tau respectively.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)

	brainmle <- estMLE(braingraphs, "adjMatrix") 
	brainmle
}
