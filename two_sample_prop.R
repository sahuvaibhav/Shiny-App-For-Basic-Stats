# Two Independent sample - population standard deviations are known 

two_sample_prop = function(p1,p2,dp,n1,n2,type,alpha){
  pbar = (p1*n1+p2*n2)/(n1+n2)
  qbar = 1-pbar
  s= sqrt(pbar*qbar*(1/n1+1/n2))
  p_value = round(pnorm((p1-p2),dp,s),4)
  x <- seq(-4,4,length=100)*s + dp
  hx <- dnorm(x,dp,s)
  
  if (type == 'Two Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dp",
        main="Two Tail Two Population Proportion Z-Test", axes=FALSE)
    if (p_value > 0.5) {
      p_value = 2*(1-p_value)
    } else {p_value = 2*p_value}
    t1 = paste("P-value:",p_value)
    zl = round(qnorm(alpha/2,dp,s),4)
    zu = round(qnorm(1-alpha/2,dp,s),4)
    t2 = paste("Critical values: [",zl,",",zu,"]")
    lb = round(qnorm(alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    ub = round(qnorm(1-alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    t3 = paste("Confidence Interval:[",lb, ",",ub,"]" )
    i <- x <= zl 
    j  = x >=zu
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(zl,zu,dp-3*s, dp-2*s,dp-s,dp,
                 dp+s,dp+2*s,dp+3*s,p1-p2), pos=0) 
    
    text(x=zl, y=(dnorm(dp,dp,s) - dnorm(zl,dp,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(dp,dp,s) - dnorm(zl,dp,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,dp, s), col='blue', length=0.1, lwd=1)
    
    text(x=zu, y=(dnorm(dp,dp,s) - dnorm(zu,dp,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(dp,dp,s) - dnorm(zu,dp,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,dp, s), col='blue', length=0.1, lwd=1)
    
    text(x=p1-p2, y=(dnorm(dp,dp,s)-0)/1.6, labels='Sample Mean (p1-p2)', col='Green')
    arrows(x0=p1-p2, y0=(dnorm(dp,dp,s)-0)/1.6, x1=p1-p2, y1=0, col='Green', length=0.1, lwd=1)
    
  }
  if (type == 'Left Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dp",
         main="Left Tail Two Population Proportion Z-Test", axes=FALSE)
    t1 = paste("P-value:",p_value)
    zl = round(qnorm(alpha,dp,s),4)
    t2 = paste("Critical values: [",zl,"]")
    lb = round(qnorm(alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    ub = round(qnorm(1-alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    t3 = paste("Confidence Interval:[",lb, ",",ub,"]" )
    i <- x <= zl 
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    axis(1, at=c(dp-3*s, dp-2*s,dp-s,dp,dp+s,dp+2*s,dp+3*s,zl,p1-p2), pos=0)
    text(x=zl, y=(dnorm(dp,dp,s) - dnorm(zl,dp,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(dp,dp,s) - dnorm(zl,dp,s))/1.2, x1=zl-0.3*s, y1=dnorm(zl-0.3*s,dp, s), col='blue', length=0.1, lwd=1)
    text(x=p1-p2, y=(dnorm(dp,dp,s) - 0)/1.6, labels='Sample Mean (p1-p2)', col='Green')
    arrows(x0=p1-p2, y0=(dnorm(dp,dp,s) - 0)/1.6, x1=p1-p2, y1=0, col='Green', length=0.1, lwd=1)
  
  }
  if (type == 'Right Tail') {
    plot(x, hx, type="n", ylab="",xlab = "dp",
         main="Right Tail Two Population Proportion Z-Test", axes=FALSE)
    t1 = paste("P-value:",1-p_value)
    zu = round(qnorm(1-alpha,dp,s),4)
    t2 = paste("Critical values: [",zu,"]")
    lb = round(qnorm(alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    ub = round(qnorm(1-alpha/2, p1-p2,sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)),4)
    t3 = paste("Confidence Interval:[",lb, ",",ub,"]" )
    j  = x >=zu
    lines(x, hx)
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(dp-3*s, dp-2*s,dp-s,dp,
                 dp+s,dp+2*s,dp+3*s,zu,p1-p2), pos=0)
    text(x=zu, y=(dnorm(dp,dp,s) - dnorm(zu,dp,s))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(dp,dp,s) - dnorm(zu,dp,s))/1.2, x1=zu+0.3*s, y1=dnorm(zu+0.3*s,dp, s), col='blue', length=0.1, lwd=1)
    text(x=p1-p2, y=(dnorm(dp,dp,s) - 0)/1.6, labels='Sample Mean (p1-p2)', col='Green')
    arrows(x0=p1-p2, y0=(dnorm(dp,dp,s) - 0)/1.6, x1=p1-p2, y1=0, col='Green', length=0.1, lwd=1)
  }  
  return(list(t1,t2,t3))
  
  
}
##### Problems ####

