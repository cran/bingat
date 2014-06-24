glrtPvalue <-
function(data, type, groups, PBmethod=FALSE, bootstraps=10, cutoff=.5){
if(missing(data) || missing(type) || missing(groups))
stop("data, type, and/or groups is missing.")

if(bootstraps <= 0)
stop("bootstraps must be an integer greater than 0.")

#Make sure the data is in the proper format or if possible convert it to that format 
if(class(data) == "list" && (tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt")){
data <- do.call(cbind, data)
}else if(length(dim(data)) == 2 && (tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt" || tolower(type) == "diag")){
data <- data
}else{
stop(sprintf("%s is unknown and/or the dimensions of data are unknown.", as.character(type)))
}

reg <- glmReg(data, type, groups, cutoff)
glrt <- list(glrtReg(data, type, groups, cutoff))
names(glrt) <- "GLRT"
bootsample <- 1:bootstraps

if(PBmethod){ #Parametric Bootstrap
reg.coeff.H0 <- glmReg(data, type, rep(0, ncol(data)), cutoff)
b0.H0 <- reg.coeff.H0$b0.covs0
tau.H0 <- reg.coeff.H0$tau

GLRT.bootstrap <- apply(as.matrix(bootsample), 1, function(x, gstar, tau, N, grps, type, co){
data.MC <- rGibbs(gstar, tau, type, N)
glrt.mc <- glrtReg(data.MC, type, grps, co)
return(glrt.mc)
},gstar=b0.H0, tau=tau.H0, N=ncol(data), grps=groups, type=type, co=cutoff)
}else{ #Nonparametric Bootstrap: Resampling with replacement
x1 <- sum(groups==1)
x2 <- sum(groups==0)

GLRT.bootstrap <- apply(as.matrix(bootsample), 1, function(x, x1, x2, data, grps, type, co){
p1 <- sample(x1+x2, replace=TRUE)
data1b <- data[,p1[c(1:x1)]]
data2b <- data[,p1[c((x1+1):(x1+x2))]]
dataMCnp <- cbind(data1b, data2b)
glrt.mcnp <- glrtReg(dataMCnp, type, grps, co)
return(glrt.mcnp)
}, x1=x1, x2=x2, data=data, grps=groups, type=type, co=cutoff)
}

p.value <- list(sum(GLRT.bootstrap >= glrt)/bootstraps)
names(p.value) <- "pvalue"

#Label our output list
l <- list(reg, glrt, p.value)
keys <- unique(unlist(lapply(l, names)))
results<- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)

return(results)
}
