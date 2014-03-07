estMLE <-
function(data, type, cutoff = .5){
gstar <- estGStar(data, cutoff)
tau <- estTau(data, type, gstar)

return(list(gstar=gstar, tau=tau))
}
