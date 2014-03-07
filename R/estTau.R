estTau <-
function(data, type, gstar){
if(missing(data) || missing(type) || missing(gstar))
stop("data, type, and/or gstar is missing.")

nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)

if(class(gstar) == "data.frame" || class(gstar) == "matrix") #Check gstar is a single vector
gstar <- gstar[,1]

distToGStar <- apply(data, 2, function(x, g, type) {calcDistance(x, g, type)}, g=gstar, type=type)

sumDist <- sum(distToGStar)
sumDist <- ifelse(sumDist==0, .Machine$double.xmin, sumDist) #Adjust sumDist if it is 0

num <- (ncol(data) * edges)^(-1) * sumDist
den <- 1-num
tau <- -log(num/den)

return(tau)
}
