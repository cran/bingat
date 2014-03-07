\name{glrtReg}
\alias{glrtReg}
\title{General Likelihood Ratio Test for Graphs}

\description{This function performs the general likelihood ratio test.}

\usage{glrtReg(data, type, groups, cutoff = .5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{groups}{A vector of 0/1s used to denote the starting group splits.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\value{The GLRT value for the \code{glmReg}.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	
	reg <- glrtReg(braingraphs, "adjMatrix", grps)
	reg
}
