DimReductionlist <- setClass("DimReductionlist", slots=c(Components = "data.frame",cum.var = "numeric", sdev="numeric", rotation="matrix",center= "numeric", scale="ANY"))

Set.DimReductionlist <- function(O){
Comp.df <- as.data.frame(O$x)
sd <-O$sdev
names(sd) <- colnames(Comp.df)
prop.var <- (O$sdev)^2/sum((O$sdev)^2)
names(prop.var) <- colnames(Comp.df)
Object <- DimReductionlist(Components =Comp.df, cum.var = prop.var, sdev = sd, rotation = O$rotation, center= O$center, scale = O$scale)
return(Object)
}
