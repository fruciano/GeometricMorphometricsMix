### Ancillary functions
### (not exported)


rep.row=function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col=function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}