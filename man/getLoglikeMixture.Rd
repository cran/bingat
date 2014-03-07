\name{getLoglikeMixture}
\alias{getLoglikeMixture}
\title{Group Finder}

\description{This function takes group splits and determines the likelihood of those groups.}

\usage{getLoglikeMixture(data, mixture)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{mixture}{The output of the \code{getGibbsMixture} function.}
}

\value{A list containing the BIC criteria and the log likelihood named bic and ll respectively.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5)
	brainlm <- getLoglikeMixture(braingraphs, braingm)
	brainlm
	
	### By running the loglik mixture over several groups you can find which is the optimal
	\dontrun{
		mixtures <- NULL
		for(i in 1:5){
			tempgm <- getGibbsMixture(braingraphs, "adjMatrix", i)
			mixtures[i] <- getLoglikeMixture(braingraphs, tempgm)$bic
		}
		
		bestgroupnum <- which(min(mixtures) == mixtures)
		bestgroupnum
	}
}
