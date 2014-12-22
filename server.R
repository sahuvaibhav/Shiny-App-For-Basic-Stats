
source("pooled.R")
source('unpooled.R')
source('two_sample_prop.R')
source('Paired_sample.R')
options("scipen"=100, "digits"=4)
library("pastecs")
library("psych")
library("car")
library("xtable")
normal = function(mean, sd, lb,ub){
  x <- seq(-4,4,length=100)*sd + mean
  hx <- dnorm(x,mean,sd)
  
  plot(x, hx, type="n", ylab="",
       main="Normal Distribution", axes=FALSE)
  
  i <- x >= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red")
  
  area <- pnorm(ub, mean, sd) - pnorm(lb, mean, sd)
  result <- paste("P(",lb,"< Mu <",ub,") =",
                  signif(area, digits=3))
  mtext(result,3)
  axis(1, at=seq(mean-4*sd, mean+4*sd, sd), pos=0) 
}
hypothesis_test = function(mean,SmeanHT,sd,sdchoice,type,alpha,n){
  x <- seq(-4,4,length=100)*sd/sqrt(n) + mean
  hx <- dnorm(x,mean,sd/sqrt(n))

  if (type == '1') {
    plot(x, hx, type="n", ylab="",xlab = "Xbar",main="Two Tail Hypothesis Test", axes=FALSE)
    if (sdchoice == '1'){
      zl = qnorm(alpha/2,mean,sd/sqrt(n))
      zu = qnorm(1-alpha/2,mean,sd/sqrt(n))
    } else {
      zl = mean+ qt(alpha/2,n-1)*sd/sqrt(n)
      zu = mean+qt(1-alpha/2,n-1)*sd/sqrt(n)
    }
    i <- x <= zl 
    j  = x >=zu
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(zl,zu,mean-3*sd/sqrt(n), mean-2*sd/sqrt(n),mean-sd/sqrt(n),mean,
                 mean+sd/sqrt(n),mean+2*sd/sqrt(n),mean+3*sd/sqrt(n),SmeanHT), pos=0) 
    
    text(x=zl, y=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zl,mean,sd/sqrt(n)))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zl,mean,sd/sqrt(n)))/1.2, x1=zl-0.3*sd/sqrt(n), y1=dnorm(zl-0.3*sd/sqrt(n),mean, sd/sqrt(n)), col='blue', length=0.1, lwd=1)
    
    text(x=zu, y=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zu,mean,sd/sqrt(n)))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zu,mean,sd/sqrt(n)))/1.2, x1=zu+0.3*sd/sqrt(n), y1=dnorm(zu+0.3*sd/sqrt(n),mean, sd/sqrt(n)), col='blue', length=0.1, lwd=1)
    
    text(x=SmeanHT, y=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, labels='Sample Mean', col='Green')
    arrows(x0=SmeanHT, y0=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, x1=SmeanHT, y1=0, col='Green', length=0.1, lwd=1)
  }
  if (type == '2') {
    plot(x, hx, type="n", ylab="",xlab = "Xbar",
         main="Left Tail Hypothesis Test", axes=FALSE)
    if (sdchoice == '1'){
      zl = qnorm(alpha,mean,sd/sqrt(n))
    } else {
      zl = mean+ qt(alpha,n-1)*sd/sqrt(n)
    }
    i <- x <= zl 
    lines(x, hx)
    polygon(c(min(x),x[i],zl), c(0,hx[i],0), col="red")
    axis(1, at=c(mean-3*sd/sqrt(n), mean-2*sd/sqrt(n),mean-sd/sqrt(n),mean,
                 mean+sd/sqrt(n),mean+2*sd/sqrt(n),mean+3*sd/sqrt(n),zl,SmeanHT), pos=0)
    text(x=zl, y=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zl,mean,sd/sqrt(n)))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zl, y0=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zl,mean,sd/sqrt(n)))/1.2, x1=zl-0.3*sd/sqrt(n), y1=dnorm(zl-0.3*sd/sqrt(n),mean, sd/sqrt(n)), col='blue', length=0.1, lwd=1)
    text(x=SmeanHT, y=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, labels='Sample Mean', col='Green')
    arrows(x0=SmeanHT, y0=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, x1=SmeanHT, y1=0, col='Green', length=0.1, lwd=1)
  }
  
  if (type == '3') {
    plot(x, hx, type="n", ylab="",xlab = "Xbar",
         main="Right Tail Hypothesis Test", axes=FALSE)
    if (sdchoice == '1'){
      zu = qnorm(1-alpha,mean,sd/sqrt(n))
    } else {
      zu = mean+qt(1-alpha,n-1)*sd/sqrt(n)
    }
    j  = x >=zu
    lines(x, hx)
    polygon(c(zu,x[j],max(x)), c(0,hx[j],0), col="red")
    axis(1, at=c(mean-3*sd/sqrt(n), mean-2*sd/sqrt(n),mean-sd/sqrt(n),mean,
                 mean+sd/sqrt(n),mean+2*sd/sqrt(n),mean+3*sd/sqrt(n),zu,SmeanHT), pos=0)
    text(x=zu, y=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zu,mean,sd/sqrt(n)))/1.2, labels='Rejection Region', col='blue')
    arrows(x0=zu, y0=(dnorm(mean,mean,sd/sqrt(n)) - dnorm(zu,mean,sd/sqrt(n)))/1.2, x1=zu+0.3*sd/sqrt(n), y1=dnorm(zu+0.3*sd/sqrt(n),mean, sd/sqrt(n)), col='blue', length=0.1, lwd=1)
    text(x=SmeanHT, y=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, labels='Sample Mean', col='Green')
    arrows(x0=SmeanHT, y0=(dnorm(mean,mean,sd/sqrt(n)) - 0)/1.6, x1=SmeanHT, y1=0, col='Green', length=0.1, lwd=1)
  }
}

hyp_test_res = function(mean,SmeanHT,sd,sdchoice,type,alpha,n){
  if (type == '1'){
    if (sdchoice == '1'){
      zl = round(qnorm(alpha/2,mean, sd/sqrt(n)),5)
      zu = round(qnorm(1-alpha/2,mean, sd/sqrt(n)),5)
      p_value = round(pnorm(SmeanHT,mean,sd/sqrt(n)),5)
          if (p_value >0.5){
              p_value = 2*(1- p_value)
          } else {
            p_value = 2*p_value }
    } else {
      zl = round(mean + qt(alpha/2,n-1)*sd/sqrt(n),5)
      zu = round(mean + qt(1-alpha/2,n-1)*sd/sqrt(n),5)
      p_value = round(pt((SmeanHT-mean)*sqrt(n)/sd, n-1),5)
      if (p_value >0.5){
        p_value = 2*(1- p_value)
      } else {
        p_value = 2*p_value }
    }
    if (p_value <= alpha)    {
      text = paste("As p-value is less than alpha(for 2 tail test) we can reject null hypothesis")        
    }
    if (p_value > alpha)    {
      text = paste("As p-value is greater than alpha(for 2 tail test) we can NOT reject null hypothesis")        
    }
    text0 = paste("p-value:" , p_value)
    text2 = paste("Critical Values:" , zl ,"," ,zu)

  }
 
  if (type == '2'){
    if (sdchoice == "1"){
      zl = round(qnorm(alpha,mean, sd/sqrt(n)),5)
      p_value = round(pnorm(SmeanHT,mean,sd/sqrt(n)),5)
    } else {
      zl = round(mean + qt(alpha,n-1)*sd/sqrt(n),5)
      p_value = round(pt((SmeanHT-mean)*sqrt(n)/sd, n-1),5)
    } 
      
    if (p_value <= alpha)    {
      text = paste("As p-value is less than alpha we can reject the null hypothesis(H0).")        
    }
    if (p_value > alpha)    {
      text = paste("As p-value is greater than alpha we can NOT reject the null hypothesis(H0)")        
    }
    text0 = paste("p-value:" , p_value)
    text2 = paste("Critical Value:" , zl)
    
  }
  if (type == '3'){
    if(sdchoice == '1') {
      zu = round(qnorm(1-alpha,mean, sd/sqrt(n)),5)
      p_value =round(1- pnorm(SmeanHT,mean,sd/sqrt(n)),5)
    } else {
      zu = round(mean + qt(1-alpha,n-1)*sd/sqrt(n),5)
      p_value =round(1- pt((SmeanHT-mean)*sqrt(n)/sd, n-1),5)
   
    }
    if (p_value <= alpha)    {
      text = paste("As p-value is less than alpha we can reject the null hypothesis(H0).")        
    }
    if (p_value > alpha)    {
      text = paste("As p-value is greater than alpha we can NOT reject the null hypothesis(H0)")        
    }
    text0 = paste("p-value:" , p_value)
    text2 = paste("Critical Value:" ,zu)
  }
  
  return(c(text0,text,text2))
  
}

shinyServer(function(input, output){
   output$mapCI = renderPlot({normal(input$meanCI,input$sdCI,input$lbCI,input$ubCI)})
   output$mapHT = renderPlot({hypothesis_test(input$meanHT,input$SmeanHT,input$sdHT,input$sdchoice,input$type,input$alpha,input$n)})
   alt_hyp = reactive({
     if (input$type == '1'){
         return(paste("Alternate Hypothesis Ha: Mu != ", input$meanHT))
     }
     if (input$type == '2'){
       return(paste("Alternate Hypothesis Ha: Mu < ", input$meanHT))
     }
     if (input$type == '3'){
       return(paste("Alternate Hypothesis Ha: Mu > ", input$meanHT))
     }
   })
  
    output$H0 =  renderText({paste("Null Hypothesis H0: Mu = ", input$meanHT)})
    output$Ha =  renderText({alt_hyp()})
    output$result1 =  renderText({hyp_test_res(input$meanHT,input$SmeanHT,input$sdHT,input$sdchoice,input$type,input$alpha,input$n)[1]})
    output$result2 =  renderText({hyp_test_res(input$meanHT,input$SmeanHT,input$sdHT,input$sdchoice,input$type,input$alpha,input$n)[2]})
    output$result3 =  renderText({hyp_test_res(input$meanHT,input$SmeanHT,input$sdHT,input$sdchoice,input$type,input$alpha,input$n)[3]})
   
## Hypotheis Test for two populations start ##
   
   text_ind = reactive({
     if(input$VarType == "Pooled Variance"){
       return(pooled(input$x1,input$x2,input$du,input$s1,input$s2,input$n1,
                     input$n2,input$itype,input$alphaI))  
     }
     if(input$VarType == "Unpooled Variance"){
       return(unpooled(input$x1,input$x2,input$du,input$s1,input$s2,input$n1,
                       input$n2,input$itype,input$alphaI)) }
   })
   
   plot_ind = reactive({
     if(input$VarType == "Pooled Variance"){
       pooled(input$x1,input$x2,input$du,
              input$s1,input$s2,input$n1,input$n2,input$itype,input$alphaI)  }
     if(input$VarType == "Unpooled Variance"){
       unpooled(input$x1,input$x2,input$du,input$s1,input$s2,input$n1,
                input$n2,input$itype,input$alphaI) }
   })
   output$plotIndependent = renderPlot({plot_ind()})
   output$textIndependent1 = renderText({text_ind()[[1]]})
   output$textIndependent2 = renderText({text_ind()[[2]]})
   output$textIndependent3 = renderText({text_ind()[[3]]})
   
   text_Pr = reactive({
     return(two_sample_prop(input$p1,input$p2,input$dp,
                            input$n1Pr,input$n2Pr,input$prtype,input$alphaPr))  
     
   })
   
   output$plotProportions = renderPlot({two_sample_prop(input$p1,input$p2,input$dp,
                                                        input$n1Pr,input$n2Pr,input$prtype,input$alphaPr)})
   output$textProp1 = renderText({text_Pr()[[1]]})
   output$textProp2 = renderText({text_Pr()[[2]]})
   output$textProp3 = renderText({text_Pr()[[3]]})
   
   text_P = reactive({
     return(Paired_sample(input$dbar,input$D,input$sdbar,input$nP,input$ptype,input$alphaP))  
     
   })
   
   output$plotPaired = renderPlot({Paired_sample(input$dbar,input$D,
                                                 input$sdbar,input$nP,input$ptype,input$alphaP)})
   output$txtPaired1 = renderText({text_P()[[1]]})
   output$txtPaired2 = renderText({text_P()[[2]]})
   output$txtPaired3 = renderText({text_P()[[3]]})
   
## Hypotheis Test for two populations end ##

## ANOVA STARTS##
data2 = reactive({
  data = data()
  var = input$variable
  return(data$var)
})

data = reactive({
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  else 
    data = read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
  return(data)
})
filterCatVar = reactive({
  data = data()
  data.temp = data
  data.temp[colnames(data)] = lapply(data[colnames(data)],as.factor)
  levels.temp = lapply(data.temp[colnames(data.temp)],levels)
  length.levels = lapply(levels.temp,length)
  length.levels = length.levels[length.levels<30]
  return(names(length.levels))
})



output$filename = renderText({
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  paste("File Name :",inFile$name, "         Size:" , inFile$size/1000, "KB",sep = " ")
})
output$contents <- renderTable({
  if (is.null(data()))
    return(NULL)
  #data = read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
  stat.desc(data(),norm = T)
  #describe(read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote))
})
output$variable = renderUI({
  #inFile <- input$file1
  if (is.null(data()))
    return(NULL)
  variables = colnames(data())
  checkboxGroupInput("variable", "Choose Variables", variables,selected = variables[[1]])
})
#     output$t = renderText(
#       paste(filterCatVar(),collapse = ",")
#       #paste(input$variable[1],length(input$variable))
#       )

output$plots <- renderUI({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$variable))
    return(NULL)
  plot_output_list <- lapply(1:length(input$variable), function(i) {
    #plotname <- paste(input$variable[[i]], i, sep="")
    plotname <- paste("plot", i, sep="")
    plotOutput(plotname, height = 350, width = 700)
  })
  
  # Convert the list to a tagList - this is necessary for the list of items
  # to display properly.
  do.call(tagList, plot_output_list)
  
})


observe({ for (i in 1:length(input$variable)) {
  
  local({
    my_i <- i
    plotname <- paste("plot",my_i, sep="")
    data = data() 
    output[[plotname]] <- renderPlot({
      par(oma = c( 2, 2, 2, 2 ) )
      nf<- layout( matrix( c( 1, 2 ), 1, 2, byrow = TRUE),c(2, 2), c( 5, 5 ), TRUE )
      layout.show(nf )
      par(mfrow = c(1,2))
      hx <- dnorm(data[[input$variable[[my_i]]]],mean(data[[input$variable[[my_i]]]]),sd(data[[input$variable[[my_i]]]])/sqrt(length(data[[input$variable[[my_i]]]])))
      boxplot(data[[input$variable[[my_i]]]],col = "red",xlab = NULL)
      hist(data[[input$variable[[my_i]]]], col = "red",main = NULL,prob = T,xlab = NULL)
      curve(dnorm(x,mean(data[[input$variable[[my_i]]]]),sd(data[[input$variable[[my_i]]]])),add = T, col = "green", lwd = 2)
      ##  Create an overall title.
      mtext( paste(input$variable[[my_i]]), outer = TRUE,cex = 1.5,col = "Blue")
      #hist(data[[input$variable[[my_i]]]])
    })
  })
}
})

output$varAnova = renderUI({
  if (is.null(data()))
    return(NULL)
  variables = colnames(data())
  catVar = filterCatVar()
  #checkboxGroupInput("variable", "Choose Variables", variables,selected = variables[[1]])
  switch(input$TypeAnova,
         "anova" = div(class = "row-fluid",
                       div(class="span5",selectInput("DepVar", "Dependent Variable",variables)),
                       div(class="span5",checkboxGroupInput("IndepVar", "Independent Variables", catVar,selected = catVar[1]))
         ),
         #                "ancova" = div(class = "row-fluid",
         #                               div(class="span5",selectInput("DepVar", "Dependent Variable",variables)),
         #                               div(class="span5",checkboxGroupInput("IndepVar", "Independent Variables", catVar,selected = catVar[1]))
         #                ),
         "manova" =  div(class = "row-fluid",
                         div(class="span5",checkboxGroupInput("MDepVar", "Dependent Variables", variables,selected = variables[[1]])),
                         div(class="span5",checkboxGroupInput("MIndepVar", "Independent Variables", catVar,selected = catVar[1]))
         ) )
  
})

output$textAnova1 = renderText({
  paste("Response Variable:",input$DepVar,sep = " ")
})
output$textAnova2 = renderText({
  paste("Independent Variables:",paste(input$IndepVar,collapse = ","),sep = " ")
})

anova.fit = reactive({
  if (is.null(input$IndepVar))
    return(NULL)
  data1 = data()
  data1[input$IndepVar] = lapply(data1[input$IndepVar],as.factor)
  data1[input$DepVar] = lapply(data1[input$DepVar],as.numeric)
  e1 = paste(input$DepVar,paste(input$IndepVar,collapse = "*"),sep = "~")
  fit = aov(formula(e1),data = data1)
  return(fit)
  
}
)
output$SummaryAnova <- renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  
  summary(anova.fit())
})

anova.plot = reactive({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  layout(matrix(c(1,2,3,4),2,2))
  plot(anova.fit())
})

output$PlotAnova = renderPlot({
  anova.plot()
})

output$TukeyPlot <- renderPlot({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  
  plot(TukeyHSD(anova.fit()))
})

output$Tukey <- renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  
  data.frame(TukeyHSD(anova.fit())[1])
})

bartlett = reactive({
  if (is.null(input$IndepVar))
    return(NULL)
  data1 = data()
  data1[input$IndepVar] = lapply(data1[input$IndepVar],as.factor)
  data1[input$DepVar] = lapply(data1[input$DepVar],as.numeric)
  e1 = paste(input$DepVar,paste("interaction(",paste(input$IndepVar,collapse = ","),")",sep = ""),sep = "~")
  b.fit = bartlett.test(formula(e1),data = data1)
  data.bartlett = cbind(b.fit[1],b.fit[2],b.fit[3],b.fit[4],b.fit[5])
  colnames(data.bartlett) = c("Bartlett's K-squared","df","p.value","data.name","method")
  return(data.bartlett)
  
})
output$Bartlett <- renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  
  bartlett()
})

lavene = reactive({
  if (is.null(input$IndepVar))
    return(NULL)
  data1 = data()
  data1[input$IndepVar] = lapply(data1[input$IndepVar],as.factor)
  data1[input$DepVar] = lapply(data1[input$DepVar],as.numeric)
  e1 = paste(input$DepVar,paste(input$IndepVar,collapse = "*"),sep = "~")
  l.fit = leveneTest(formula(e1),data = data1)
  data.lavene = cbind(l.fit[1],l.fit[2],l.fit[3])
  colnames(data.lavene) = c("df","F value","Pr(>F)")
  return(data.lavene)
})

output$Levene <- renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$IndepVar))
    return(NULL)
  lavene()
})

summary_manova = reactive({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  data1 = data()
  data1[input$MIndepVar] = lapply(data1[input$MIndepVar],as.factor)
  data1[input$MDepVar] = lapply(data1[input$MDepVar],as.numeric)
  e1 = paste(paste("cbind(",paste(input$MDepVar,collapse = ","),")",sep = ""),paste(input$MIndepVar,collapse = "*"),sep = "~")
  fit = manova(formula(e1),data = data1)
  summ.manova= summary.aov(fit)
  tables <- list()
  for (i in 1:length(summ.manova)){
    tables[i] = paste(strong("Response Variable: "),input$MDepVar[i],br(),br(),print(xtable(summ.manova[1]),type = "html"),br(),sep = "" )
  }
  
  return(lapply(tables, paste))
  
  
})
output$SummaryManova= renderUI({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  #print(summary.aov(manova.fit()))
  #   s1 = print(xtable(summary.aov(manova.fit())[1]),type = "html")
  #   s3 = print("response Variable:")
  #   s2 = print(xtable(summary.aov(manova.fit())[2]),type = "html")
  out = summary_manova()
  return(div(HTML(out),class="shiny-html-output"))
  #HTML(print(xtable(summary.aov(manova.fit())[1]),type = "html"))
  #HTML(print(xtable(summary.aov(manova.fit())[2]),type = "html"))
  # HTML(print(xtable(manova.fit()),type = "html"))
})

# output$fitManova= renderTable({
#   if (is.null(data()))
#     return(NULL)
#   manova.fit()
#   
# })

manova.fit = reactive({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  data1 = data()
  data1[input$MIndepVar] = lapply(data1[input$MIndepVar],as.factor)
  data1[input$MDepVar] = lapply(data1[input$MDepVar],as.numeric)
  e1 = paste(paste("cbind(",paste(input$MDepVar,collapse = ","),")",sep = ""),paste(input$MIndepVar,collapse = "*"),sep = "~")
  fit = manova(formula(e1),data = data1)
  return(fit)
})
output$PillaiManova= renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  #   out = print(xtable(summary(manova.fit(),test = "Pillai")),type= "html")
  #   return(div(HTML(out),class="shiny-html-output"))
  data.frame(summary(manova.fit(),test = "Pillai")[4])
})
output$WilksManova= renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  data.frame(summary(manova.fit(),test = "Wilks")[4])
})
output$HotLawManova= renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  data.frame(summary(manova.fit(),test = "Hotelling-Lawley")[4])
})
output$RoyManova= renderTable({
  if (is.null(data()))
    return(NULL)
  if (is.null(input$MIndepVar))
    return(NULL)
  if (is.null(input$MDepVar))
    return(NULL)
  data.frame(summary(manova.fit(),test = "Roy")[4])
})


## ANOVA ENDS##

}
  )
