Paired_sample = function(dbar,D,sd,n,type,alpha){
  df = n-1
  s = sd/sqrt(n)
  t = (dbar-D)/s
  p_value = round(pt(t,df),4)
  
  x <- seq(-4,4,length=100)*s + D
  hx <- dnorm(x,D,s)
  
  if (type == 'Two Tail') {
    plot(x, hx, type="n", ylab="",xlab = "D",
         main="Two Tail Paired Sample Z-Test", axes=FALSE)
    if (p_value > 0.5) {
      p_value = 2*(1-p_value)
    } else {p_value = 2*p_value}
    t1 = (paste("P-value:",round(p_value,7)))
    zl = round(D + qt(alpha/2,df)*s,4)
    zu = round(D + qt(1-alpha/2,df)*s,4)
    t2 = (paste("Critical values: [",zl,",",zu,"]"))
    lb = round(dbar + qt(alpha/2,df)*s,4)
    ub = round(dbar + qt(1-alpha/2,df)*s,4)
    t3 =  (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    
    i <- x <= zl 
    j  = x >=zu
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(zl,zu,D-3*s, D-2*s,D-s,D,
                 D+s,D+2*s,D+3*s,dbar), pos=0) 
    
    text(x=zl, y=(dnorm(D,D,s) - dnorm(zl,D,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(D,D,s) - dnorm(zl,D,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,D, s), col='blue', length=0.1, lwd=1)
    
    text(x=zu, y=(dnorm(D,D,s) - dnorm(zu,D,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(D,D,s) - dnorm(zu,D,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,D, s), col='blue', length=0.1, lwd=1)
    
    text(x=dbar, y=(dnorm(D,D,s)-0)/1.6, labels='Sample Mean (dbar)', col='Green')
    arrows(x0=dbar, y0=(dnorm(D,D,s)-0)/1.6, x1=dbar, y1=0, col='Green', length=0.1, lwd=1)
    
  }
  if (type == 'Left Tail') {
    plot(x, hx, type="n", ylab="",xlab = "D",
         main="Left Tail Paired Sample Z-Test", axes=FALSE)
    t1 = (paste("P-value:",p_value))
    zl = round(D + qt(alpha,df)*s,4)
    t2= (paste("Critical values: [",zl,"]"))
    lb = round((dbar) + qt(alpha/2,df)*s,4)
    ub = round((dbar) + qt(1-alpha/2,df)*s,4)
    t3 = (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    i <- x <= zl 
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    axis(1, at=c(D-3*s, D-2*s,D-s,D,D+s,D+2*s,D+3*s,zl,dbar), pos=0)
    text(x=zl, y=(dnorm(D,D,s) - dnorm(zl,D,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(D,D,s) - dnorm(zl,D,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,D, s), col='blue', length=0.1, lwd=1)
    text(x=dbar, y=(dnorm(D,D,s) - 0)/1.6, labels='Sample Mean (dbar)', col='Green')
    arrows(x0=dbar, y0=(dnorm(D,D,s) - 0)/1.6, x1=dbar, y1=0, col='Green', length=0.1, lwd=1)
    
  }
  if (type == 'Right Tail') {
    plot(x, hx, type="n", ylab="",xlab = "D",
         main="Right Tail Paired Sample Z-Test", axes=FALSE)
    t1 = (paste("P-value:",1-p_value))
    zu = round(D + qt(1-alpha,df)*s,4)
    t2= (paste("Critical values: [",zu,"]"))
    lb = round((dbar) + qt(alpha/2,df)*s,4)
    ub = round((dbar) + qt(1-alpha/2,df)*s,4)
    t3 = (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    j  = x >=zu
    lines(x, hx)
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(D-3*s, D-2*s,D-s,D,
                 D+s,D+2*s,D+3*s,zu,dbar), pos=0)
    text(x=zu, y=(dnorm(D,D,s) - dnorm(zu,D,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(D,D,s) - dnorm(zu,D,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,D, s), col='blue', length=0.1, lwd=1)
    text(x=dbar, y=(dnorm(D,D,s) - 0)/1.6, labels='Sample Mean (dbar)', col='Green')
    arrows(x0=dbar, y0=(dnorm(D,D,s) - 0)/1.6, x1=dbar, y1=0, col='Green', length=0.1, lwd=1)
  }
  
  return(list(t1,t2,t3))  
  
}


