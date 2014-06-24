glrtReg <-
function(data, type, groups, cutoff=.5){
regCoeffHA <- glmReg(data, type, groups, cutoff)
loglikHA <- regCoeffHA$loglik

regCoeffH0 <- glmReg(data, type, rep(0, ncol(data)), cutoff)
loglikH0 <- regCoeffH0$loglik

glrt <- 2 * (loglikHA-loglikH0)
return(glrt)
}
