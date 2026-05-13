

sqfunc<-function(x,b,s){
  (1/2)*sqrt(b*(4*s+b*x^2))
}

sat<-function(x,a,b,B2,s=1E-10){
  a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s)
}

#B2 = breakpoint
#start=c(B2=mean(x)+(max(x)-mean(x))/2,a=a.guess,b=round(max(c(slopes,0.0001)),5)),

curve(sat(x,a=2,b=2,B2=4),0,8)

a<-2
b<-2
B2<-4
s=1E-10
b2<-3

curve(a+0*x,-2,8,ylim=c(-5,10),col='blue')
curve((1/2)*b*(B2)+0*x,-2,8,add=T,col='red')
curve(sqfunc(-x,b,s)+0*x,-2,8,add=T,col='purple')
curve(-sqfunc(B2-x,b,s)+0*x,-2,8,add=T,col='green')

curve(-(b2/2)*(x-B2),-2,8,add=T,col='orange')


curve(-sqfunc(B2-x,b2,s)+0*x,-2,8,add=T,col='pink')



satdecay<-function(x,a,b,b2,B2,s=1E-10){
  a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s) - (-b2/2)*(x-B2) - sqfunc(B2-x,-b2,s)
}


curve(sat(x,a=2,b=2,B2=4),0,8)
curve(satdecay(x,a=2,b=2,b2=-1,B2=4),0,8,add=T,col='purple')



