shinyUI(navbarPage("Stats",
  tabPanel("Confidence Interval",
  
  sidebarLayout(
    sidebarPanel(
      p("Provide the mean, Standard Deviation, 
        Upper Bound and Lower Bound to plot confidence Interval on Normal Plot"),
      numericInput("meanCI", label = h4("Mean"), value = 0),
      numericInput("sdCI", label = h4("Standard Deviation(not Variance)"), value = 1),
      numericInput("lbCI", label = h4("Lower Bound"), value = -1),
      numericInput("ubCI", label = h4("Upper Bound"), value = 1)
      
      ),
    mainPanel(plotOutput("mapCI"))
      
    )
    ),
  tabPanel("Hypthesis Testing",
           
           sidebarLayout(
             sidebarPanel(h4("Hypothesis Testing for Population Mean"),
               p("Select the type of your Hypothesis Test"),
               radioButtons("type", label = "",
                            choices = list("Two Tail" = 1, "Left Tail" = 2, "Right Tail" = 3),
                            selected = 1),
               numericInput("meanHT", label = h4("Population Mean(Mu)"), value = 0),
               numericInput("SmeanHT", label = h4("Sample Mean(Xbar)"), value = 0),
               radioButtons("sdchoice", label = h4("Standard Deviation"),
                            choices = list("Population Standard Deviation(Sigma)" = 1, 
                                           "Sample Standard Deviation(s)" = 2),
                            selected = 1),
               numericInput("sdHT", label = "",value = 1),
               numericInput("alpha", label = h4("Alpha"), value = 0.05),
               numericInput("n", label = h4("Sample Size(n)"), value = 100)
               
             ),
             mainPanel(p(strong(span("Z-test:",style = "color:blue")), "We are applying Z-test for hypotheis testing 
                       when population standard deviation is known."),
             p(strong(span("t-test:",style = "color:blue")), "We are applying t-test for hypotheis testing 
                       when population standard deviation is nor known.",
                  span("Strictly assuming population is normally distribured.",
                       style = "color:blue")),
               plotOutput("mapHT"),
                       textOutput("H0"),
                       textOutput("Ha"),
                       h4("Test Results"),
                       textOutput("result1"),
                       textOutput("result2"),
                       textOutput("result3"))
  )),
  tabPanel("Hypothesis Testing Two Populations",
    # Sidebar with a slider input
    #sidebarLayout(
    #sidebarPanel(h4("Select Comparison Method"),
                 tabsetPanel(
                   tabPanel(h6(p("Paired Sample")),
                            sidebarPanel(
                              radioButtons("ptype",h5("Tail Type"),c("Two Tail","Left Tail","Right Tail"),selected = "Two Tail"),
                              numericInput("dbar", label = "dbar", value = 0.0),
                              numericInput("D", label = "D", value = 0.0),
                              numericInput("sdbar", label = "sd", value = 1.0),
                              numericInput("nP", label = "n", value = 100.0),
                              numericInput("alphaP", label = "alpha", value = 0.05)
                              
                            ),
                            mainPanel(plotOutput("plotPaired"),
                                      textOutput("txtPaired1"),
                                      textOutput("txtPaired2"),
                                      textOutput("txtPaired3"))
                   ),
                   tabPanel(h6(p("Independent Sample")),
                            sidebarPanel(
                              div(class = "row-fluid",
                                  div(class="span5",radioButtons("itype",h5("Tail Type"),c("Two Tail","Left Tail","Right Tail"),selected = "Two Tail")),
                                  div(class="span5",radioButtons("VarType",h5("Variance Type"),c("Pooled Variance","Unpooled Variance"),selected = "Pooled Variance"))
                              ),
                              #radioButtons("VarType",h5("Variance Type"),c("Pooled Variance","Unpooled Variance"),selected = "Pooled Variance"),
                              div(class="row-fluid",
                                  div(class="span5",numericInput("x1", label = "x1", value = 2.0)),
                                  div(class="span5",numericInput("x2", label = "x2", value = 1.0))
                              ),
                              div(class="row-fluid",
                                  div(class="span5",numericInput("s1", label = "s1", value = 1.0)),
                                  div(class="span5",numericInput("s2", label = "s2", value = 1.0))
                              ),
                              div(class="row-fluid",
                                  div(class="span5",numericInput("n1", label = "n1", value = 100.0)),
                                  div(class="span5",numericInput("n2", label = "n2", value = 100.0))
                              ),
                              numericInput("du", label = "du", value = 0.0),
                              numericInput("alphaI", label = "alpha", value = 0.05)
                            ),
                            mainPanel(plotOutput("plotIndependent"),
                                      textOutput("textIndependent1"),
                                      textOutput("textIndependent2"),
                                      textOutput("textIndependent3"))
                   ),
                   
                   tabPanel(h6("Proportions"),
                            sidebarPanel(
                              radioButtons("prtype",h5("Tail Type"),c("Two Tail","Left Tail","Right Tail"),selected = "Two Tail"),
                              div(class="row-fluid",
                                  div(class="span5",numericInput("p1", label = "p1", value = 0.2)),
                                  div(class="span5",numericInput("p2", label = "p2", value = 0.3))
                              ),
                              div(class="row-fluid",
                                  div(class="span5",numericInput("n1Pr", label = "n1", value = 100.0)),
                                  div(class="span5",numericInput("n2Pr", label = "n2", value = 100.0))
                              ),
                              numericInput("dp", label = "dp", value = 0.0),
                              numericInput("alphaPr", label = "alpha", value = 0.05)
                              
                            ),
                            mainPanel(plotOutput("plotProportions"),
                                      textOutput("textProp1"),
                                      textOutput("textProp2"),
                                      textOutput("textProp3"))
                   )
                   
                   
    #             ))
    ,
    # Show a plot of the generated distribution
    mainPanel())
    
  ),
  tabPanel("ANOVA",
           tabsetPanel(tabPanel(h6(p("Data Load & Summary")),
                                sidebarPanel("Load CSV file",
                                             fileInput('file1',"",
                                                       accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                             #                                       tags$hr(),
                                             checkboxInput('header', 'Header', TRUE),
                                             div(class = "row-fluid",
                                                 div(class="span5",radioButtons("sep",h5("Separator"),c(Comma=',',Semicolon=';',Tab='\t'),'Comma')),
                                                 div(class="span5",radioButtons("quote",h5("Quote"),c(None='','Double Quote'='"','Single Quote'="'"),'Double Quote'))
                                             ),
                                             uiOutput("variable")
                                ),
                                mainPanel(
                                  tabsetPanel(tabPanel("Summary",
                                                       textOutput('filename'),
                                                       tableOutput('contents'),
                                                       strong(h5(div("Notes:",style = "color:red"))),
                                                       
                                                       p("For more Info on Data Summary refer",a("Link.",href = "https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/R/stat.desc.R?revision=2&root=pastecs&pathrev=4")),
                                                       p(span("- SEMean <- StdDev/sqrt(Nbrval)",style = "color:blue")),
                                                       p(span("- CI.mean.0.95 = confidence interval on the mean at 95%",style = "color:blue")),
                                                       p(span("- CIMean <- qt((0.5+p/2), (Nbrval-1))*SEMean",style = "color:blue")),
                                                       p(span("- Coef.var = Coeffficient of Variation = std.dev/mean",style = "color:blue")),
                                                       p(span("- Skew <- sum((x-mean(x))^3)/(length(x)*sqrt(var(x))^3)        # From e1071 R library",style = "color:blue")),
                                                       p(span("- Kurt <- sum((x-mean(x))^4)/(length(x)*var(x)^2) - 3",style = "color:blue")),
                                                       p(span("- SE <- sqrt(6*Nbrval*(Nbrval-1)/(Nbrval-2)/(Nbrval+1)/(Nbrval+3))",style = "color:blue")),
                                                       p(span("- Skew.2SE <- Skew/(2*SE)    if skew.2SE > 1 then skewness is significantly different than zero",style = "color:blue")),
                                                       p(span("- SE <- sqrt(24*Nbrval*((Nbrval-1)2)/(Nbrval-3)/(Nbrval-2)/(Nbrval+3)/(Nbrval+5))",style = "color:blue")),
                                                       p(span("- Kurt.2SE <- Kurt/(2*SE)   # if Kurt.2SE > 1, then skewness is significantly different than zero",style = "color:blue")),
                                                       p(span("- normtest.W = the statistic of a Shapiro-Wilk test of normality <- shapiro.test(x)$statistic",style = "color:blue")),
                                                       p(span("- normtest.p = associated probability of Shapiro-Wilk test of normality <- shapiro.test(x)$p.value",style = "color:blue")),
                                                       p(span("- if normtest.p > 0.05 ; sample data is from normally distributed population (http://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test)",style = "color:blue"))
                                                       
                                  ),
                                  tabPanel("Plots", 
                                           uiOutput('plots')
                                           #                                     textOutput("t")
                                  )
                                  ))  
           ),
           tabPanel("Anova",
                    sidebarPanel(
                      radioButtons('TypeAnova', 'Type of Anova',
                                   c(ANOVA='anova',
                                     #                                 ANCOVA='ancova',
                                     MANOVA='manova'),
                                   'ANOVA'),
                      strong(h6(span("*** Do not use SAME Dependent and Independent Variables",style = "color:red"))),
                      uiOutput("varAnova")              
                    ),
                    mainPanel(tabsetPanel(tabPanel("ANOVA",
                                                   textOutput('textAnova1'),
                                                   textOutput('textAnova2'),
                                                   tableOutput('SummaryAnova'),
                                                   plotOutput('PlotAnova'),
                                                   strong(h5("Tukey HSD(Honestly Significant Difference)")),
                                                   tableOutput('Tukey'),
                                                   plotOutput('TukeyPlot'),
                                                   strong(h5("Bartlett test of homogeneity of variances")),
                                                   h6("***There must be at least 2 observations in each group"),
                                                   tableOutput('Bartlett'),
                                                   strong(h5("Levene's Test for Homogeneity of Variance")),
                                                   h6("***Levene's test is not appropriate with quantitative explanatory variables"),
                                                   tableOutput('Levene')
                                                   
                                                   
                    ),
                    tabPanel("MANOVA",
                             uiOutput("SummaryManova"),
                             strong(h5("Pillai's Test")),
                             tableOutput("PillaiManova"),
                             strong(h5("Wilks's Test")),
                             tableOutput("WilksManova"),
                             strong(h5("Hotelling-Lawley's Test")),
                             tableOutput("HotLawManova"),
                             strong(h5("Roy's Test")),
                             tableOutput("RoyManova")
                             
                    )
                    )
                    
                    )
           )
           ), mainPanel()
  

    
    )
))
