gety<-function(x,z){
  y<-(z-x)^2-x^2
  y<-y^.5
  return(y)
}

x<-seq(0,1,0.01)
y<-gety(x,2)
x<-c(x,rev(x))
y<-c(y,rev(-y))
x<-c(x,-x)
y<-c(y,rev(y))



plot(y,x,type="line",ylab="beta 2",xlab="beta 1",main="Contour plot for beta_1 and beta_2")
