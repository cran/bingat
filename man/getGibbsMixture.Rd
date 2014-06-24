\name{getGibbsMixture}
\alias{getGibbsMixture}
\title{Group Splitter}

\description{This function splits the data into groups based on the gibbs criteria.}

\usage{getGibbsMixture(data, type, desiredGroups, maxIter = 50, digits = 3, cutoff = .5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{desiredGroups}{The number of groups to test for.}
  \item{maxIter}{The maximum number of iterations to run searching for an optimal split.}
  \item{digits}{The number of digits to round internal values to when checking the stop criteria.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\value{
	A list that contains information about the group splits.  The list contains the final weights,
	gstars and taus for every group, a boolean indicating convergence, the number of iterations
	it took, and the group for each graph.
}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)

	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5) 
}
