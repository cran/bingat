graphNetworkPlot <-
function(data, type, main="Network Plot", labels, groupCounts, groupLabels){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

#Set up plot dot colors
if(missing(groupCounts)){
myColors <- "red"
}else{
allColors <- rainbow((length(groupCounts)+1))
myColors <- NULL
for(i in 1:length(groupCounts))
myColors <- c(myColors, rep(allColors[i], groupCounts[i]))
}

#Take only the first column of data if it is multi columned
if(!is.numeric(data))
data <- data[,1]

y <- vec2mat(data, type)
g <- network(as.matrix(y), directed=FALSE)
if(!missing(labels))
network.vertex.names(g) <- as.data.frame(labels)

plot(g, mode="circle", vertex.col=myColors, label=network.vertex.names(g), main=main, edge.col="black")

if(!missing(groupLabels))
legend("topright", inset=.005, legend=groupLabels, fill=allColors, horiz=FALSE)
}
