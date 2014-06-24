glmReg <-
function(data, type, groups, cutoff=.5){
if(dim(table(groups)) != 2){
gstar <- estGStar(data, cutoff)
b0 <- gstar
b0b1 <- NULL
b1 <- NULL
hammingError <- NULL
tau <- estTau(data, type, gstar)
loglik <- estLogLik(data, type, gstar, tau)
}else{
if(max(groups) != 1 || min(groups) != 0)
stop("'groups' must use 0 and 1 to denote groups.")

#Estimate the gstars for each group
b0 <- estGStar(data[,groups==0], cutoff)  
b0b1 <- estGStar(data[,groups==1], cutoff)  
b1 <- xor(b0, b0b1)*1

index <- as.matrix(1:length(groups))
hammingError <- apply(index, 1, function(x, b0, b1, grps, datap, typ){
calcDistance(datap[,x], xor(b0, b1*grps[x]), typ)
}, b0=b0, b1=b1, grps=groups, datap=data, typ=type)

#Get the tau and loglik values
gstar <- matrix(0, nrow(data), ncol(data))
for(i in 1:length(groups))
gstar[,i] <- as.matrix(xor(b0, b1*groups[i])) 

tau <- estTau(data, type, gstar)
loglik <- estLogLik(data, type, gstar, tau)
}

results <- list(b0, b1, b0b1, hammingError, loglik, tau)
names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")

return(results)
}
