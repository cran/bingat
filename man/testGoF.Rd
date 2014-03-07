\name{testGoF}
\alias{testGoF}
\title{Test the Goodness of Fit}

\description{This function tests the goodness of fit for given a set of graphs.}

\usage{testGoF(data, type, plot = TRUE, cutoff = .5, main)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{plot}{A boolean to create a plot of the results or not.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
  \item{main}{The title for the plot if applicable.}
}

\value{
	A list containing information about the goodness of fit and potentially a plot.  The list contains
	the Pearson statistics, degrees of freedom, and p-value, the G statistics and p-value, the Chi Squared
	statistics and p-value and finally the table with the observed and expected counts.
}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	gof <- testGoF(braingraphs, "adjMatrix")
}
