\name{getNumNodes}
\alias{getNumNodes}
\title{Get the Number of Nodes in a Graph}

\description{This function will return the number of nodes for a given of graph.}

\usage{getNumNodes(data, type)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
}

\value{The number of individual nodes in the given graph.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainnodes
}
