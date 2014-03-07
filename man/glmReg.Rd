\name{glmReg}
\alias{glmReg}
\title{Fitting Generalized Linear Regression Model to Graphs}

\description{This function fits a generalized linear regression model to a symmetric adjacency matrix.}

\usage{glmReg(data, type, groups, cutoff = .5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{groups}{A vector of 0/1s used to denote the starting group splits.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\value{
	A list containing the mean graph of each group (b0, b0+b1), the group differences (b1), the hamming error, 
	the log likelihood, and tau.
}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	glm <- glmReg(braingraphs, "adjMatrix", grps) 
	glm
}
