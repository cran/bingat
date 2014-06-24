\name{glrtPvalue}
\alias{glrtPvalue}
\title{GLRT Regression Results}

\description{This function returns the p-value of the significance of b1 in the regression model.}

\usage{glrtPvalue(data, type, groups, PBmethod = FALSE, bootstraps = 10, cutoff = .5)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used.}
  \item{groups}{A vector of 0/1s used to denote the starting group splits.}
  \item{PBmethod}{A boolean that determines if a parametric method of calculating the p-value should be used or not.}
  \item{bootstraps}{The number of bootstraps.}
  \item{cutoff}{The percentage of 1's a node must have to be kept in the gstar. Default is 0.50.}
}

\details{
	The non parametric method in \code{PBmethod} allows for bootstrapping to sample the data from each individual group with 
	replacement and calculates the GLRT for each bootstrapped dataset.  Then the number that are greater than original GLRT 
	result are used to calculate the p-value. \cr
	
	The parametric method in \code{PBmethod} estimates a central graph and tau of the entire data ignoring groups 
	and then randomly generates data from the Gibbs distribution. It assigns data to groups based on the given \code{covs} and 
	calculates the GLRT value and the number that are greater than original GLRT result are used to calculate the p-value.     
}

\value{A list containing the results from \code{glm.reg} and \code{GLRT.reg} as well as the p-value.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	numBootStraps <- 1 ### This is set low for speed
	
	btrue <- glrtPvalue(braingraphs, "adjMatrix", grps, TRUE, numBootStraps) 
	btrue
}
