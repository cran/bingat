glmReg <-
function(data, type, groups, cutoff=.5){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

#Make sure the data is in the proper format or if possible convert it to that format 
if(class(data) == "list" && (tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt")){
data <- do.call(cbind, data)
}else if(length(dim(data)) == 2 && (tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt" || tolower(type) == "diag")){
data <- data
}else{
stop(sprintf("%s is unknown and/or the dimensions of data are unknown.", as.character(type)))
}

if(missing(groups)){
gstar <- estGStar(data, cutoff)
b0 <- gstar
b0b1 <- NULL
b1 <- NULL
hammingError <- NULL
}else{
if(dim(table(groups)) != 2)
stop("groups must have exactly two groups.")

if(max(groups) != 1 || min(groups) != 0)
stop("groups must use 0 and 1 to denote groups.")

#Estimate the gstars for each group
b0 <- estGStar(data[,groups==0], cutoff)  
b0b1 <- estGStar(data[,groups==1], cutoff)  
b1 <- xor(b0, b0b1)*1

gstar <- b1

index <- as.matrix(1:length(groups))
hammingError <- apply(index, 1, function(x, b0, b1, grps, datap){calcDistance(datap[,x], xor(b0, b1*grps[x]))}, b0=b0, b1=b1, grps=groups, datap=data)
}

#Get the tau and logliklihood for our data
tau <- estTau(data, type, gstar)
loglik <- estLogLik(data, type, gstar, tau)

results <- list(b0, b1, b0b1, hammingError, loglik, tau)
names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")

return(results)
}
