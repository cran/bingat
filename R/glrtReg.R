glrtReg <-
function(data, type, groups, cutoff=.5){
if(missing(data) || missing(groups))
stop("data and/or grps is missing.")

regCoeffHA <- glmReg(data, type, groups, cutoff)
loglikHA <- regCoeffHA$loglik

regCoeffH0 <- glmReg(data, type, cutoff=cutoff)
loglikH0 <- regCoeffH0$loglik

glrt <- 2 * (loglikHA-loglikH0)
return(glrt)
}
