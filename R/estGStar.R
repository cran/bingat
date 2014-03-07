estGStar <-
function(data, cutoff = 0.5){
if(missing(data))
stop("data is missing.")

mean <- rowSums(data)/ncol(data)
gstar <- 1*(mean > cutoff)

return(gstar)
}
