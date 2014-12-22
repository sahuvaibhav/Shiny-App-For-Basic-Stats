pooled = function(x1,x2,du,s1,s2,n1,n2,type,alpha){
  df = n1+n2-2
  s = sqrt(((s1^2*(n1-1) + s2^2*(n2-1))/(n1+n2-2))*(1/n1+1/n2))
  t = ((x1-x2) - du)/s
  p_value = round(pt(t,df),4)

  x <- seq(-4,4,length=100)*s + du
  hx <- dnorm(x,du,s)

  if (type == 'Two Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dMu",
         main="Two Tail Independent Sample Pooled Variance Test", axes=FALSE)
    if (p_value > 0.5) {
      p_value = 2*(1-p_value)
    } else {p_value = 2*p_value}
    t1 = (paste("P-value:",round(p_value,7)))
    zl = round(du + qt(alpha/2,df)*s,4)
    zu = round(du + qt(1-alpha/2,df)*s,4)
    t2 = (paste("Critical values: [",zl,",",zu,"]"))
    lb = round((x1-x2) + qt(alpha/2,df)*s,4)
    ub = round((x1-x2) + qt(1-alpha/2,df)*s,4)
    t3 =  (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    
    i <- x <= zl 
    j  = x >=zu
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(zl,zu,du-3*s, du-2*s,du-s,du,
                 du+s,du+2*s,du+3*s,x1-x2), pos=0) 
    
    text(x=zl, y=(dnorm(du,du,s) - dnorm(zl,du,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(du,du,s) - dnorm(zl,du,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,du, s), col='blue', length=0.1, lwd=1)
                      
    text(x=zu, y=(dnorm(du,du,s) - dnorm(zu,du,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(du,du,s) - dnorm(zu,du,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,du, s), col='blue', length=0.1, lwd=1)
                             
    text(x=x1-x2, y=(dnorm(du,du,s)-0)/1.6, labels='Sample Mean (x1-x2)', col='Green')
    arrows(x0=x1-x2, y0=(dnorm(du,du,s)-0)/1.6, x1=x1-x2, y1=0, col='Green', length=0.1, lwd=1)
   
  }
  if (type == 'Left Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dMu",
         main="Left Tail Independent Sample Pooled Variance Test", axes=FALSE)
    t1 = (paste("P-value:",p_value))
    zl = round(du + qt(alpha,df)*s,4)
    t2= (paste("Critical values: [",zl,"]"))
    lb = round((x1-x2) + qt(alpha/2,df)*s,4)
    ub = round((x1-x2) + qt(1-alpha/2,df)*s,4)
    t3 = (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    i <- x <= zl 
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    axis(1, at=c(du-3*s, du-2*s,du-s,du,du+s,du+2*s,du+3*s,zl,x1-x2), pos=0)
    text(x=zl, y=(dnorm(du,du,s) - dnorm(zl,du,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(du,du,s) - dnorm(zl,du,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,du, s), col='blue', length=0.1, lwd=1)
    text(x=x1-x2, y=(dnorm(du,du,s) - 0)/1.6, labels='Sample Mean (x1-x2)', col='Green')
    arrows(x0=x1-x2, y0=(dnorm(du,du,s) - 0)/1.6, x1=x1-x2, y1=0, col='Green', length=0.1, lwd=1)
  
  }
  if (type == 'Right Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dMu",
         main="Right Tail Independent Sample Pooled Variance Test", axes=FALSE)
    t1 = (paste("P-value:",1-p_value))
    zu = round(du + qt(1-alpha,df)*s,4)
    t2= (paste("Critical values: [",zu,"]"))
    lb = round((x1-x2) + qt(alpha/2,df)*s,4)
    ub = round((x1-x2) + qt(1-alpha/2,df)*s,4)
    t3 = (paste("Confidence Interval:[",lb, ",",ub,"]" ))
    j  = x >=zu
    lines(x, hx)
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(du-3*s, du-2*s,du-s,du,
                 du+s,du+2*s,du+3*s,zu,x1-x2), pos=0)
    text(x=zu, y=(dnorm(du,du,s) - dnorm(zu,du,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(du,du,s) - dnorm(zu,du,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,du, s), col='blue', length=0.1, lwd=1)
    text(x=x1-x2, y=(dnorm(du,du,s) - 0)/1.6, labels='Sample Mean (x1-x2)', col='Green')
    arrows(x0=x1-x2, y0=(dnorm(du,du,s) - 0)/1.6, x1=x1-x2, y1=0, col='Green', length=0.1, lwd=1)
  }
    
  return(list(t1,t2,t3))  

}


