estLogLik <-
function(data, type, g, tau){
if(missing(data) || missing(type) || missing(g) || missing(tau))
stop("data, type, g, and/or tau is missing.")

nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)
normConst <- (1+exp(-tau))^(-edges)    
normConst <- ifelse(normConst==0, .Machine$double.xmin, normConst) #Adjust normConst if it is 0

if(class(g) == "data.frame" || class(g) == "matrix") #Check g is a single vector
g <- g[,1]

distToGStar <- apply(data, 2, function(x, g, type) {calcDistance(x, g, type)}, g=g, type=type)
LogLik <- ncol(data) * log(normConst) - tau * sum(distToGStar)

return(LogLik)
}
