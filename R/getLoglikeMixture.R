getLoglikeMixture <-
function(data, mixture){ 
if(missing(data) || missing(mixture))
stop("data and/or mixture is missing.")

type <- mixture$type
numGroups <- mixture$numgroups
numGraphs <- ncol(data)
nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)

#Calculate distances to the gstars
distances <- matrix(NA, ncol(data), numGroups)
for(j in 1:numGroups)
distances[,j] <- apply(data, 2, function(x){calcDistance(x, mixture$gstars[[j]], type)})

#Calculate the weighted loglikelihood values
LL <- 0
for(i in 1:numGraphs){
tempLL <- 0
for(j in 1:numGroups){
tau <- mixture$taus[j]
weight <- mixture$weights[j]
tempLL <- tempLL + weight * (1+exp(-tau))^(-edges) * exp(-tau*distances[i, j]) 

}
LL <- LL + log(tempLL)
}

BIC <- -2*LL + numGroups * (2+edges) * log(numGraphs)

return(list(ll=LL, bic=BIC))
}
