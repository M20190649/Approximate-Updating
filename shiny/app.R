library(shiny)
library(shinythemes)
library(ggplot2)
source('ADVI.R')
source('MIVB.R')
# Define UI for application
ui <- navbarPage(theme = shinytheme('cerulean'),
                 tabPanel('T0'),
                 tabPanel('Normal Model',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run1', 'Run'),
                              numericInput('M1', 'M:', min = 1, max = 500, value = 50),
                              numericInput('t1', 'T:', min = 50, max = 500, value = 200),
                              sliderInput('sigmasq1', 'Sigma Squared:', min = 0.01, max = 10, value = 1),
                              sliderInput('datasets1', 'Datasets:', min = 1, max = 4, value = 2)
                            ),
                            mainPanel(
                              plotOutput('distPlot1')
                            )
                          )
                 ),
                 tabPanel('AR1 No Mean',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run2', 'Run'),
                              numericInput('M2', 'M:', min = 1, max = 500, value = 50),
                              numericInput('t2', 'T:', min = 50, max = 500, value = 200),
                              sliderInput('sigmasq2', 'Sigma Squared:', min = 0.01, max = 10, value = 1),
                              sliderInput('phi2', 'Phi:', min = -0.99, max = 0.99, value = 0.8),
                              sliderInput('datasets2', 'Datasets:', min = 1, max = 4, value = 2)
                            ),
                            mainPanel(
                              plotOutput('distPlot2')
                            )
                          )
                 ),
                 tabPanel('AR1 With Mean',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run3', 'Run'),
                              numericInput('M3', 'M:', min = 1, max = 500, value = 50),
                              numericInput('t3', 'T:', min = 50, max = 500, value = 200),
                              sliderInput('sigmasq3', 'Sigma Squared:', min = 0.01, max = 10, value = 1),
                              sliderInput('phi3', 'Phi:', min = -0.99, max = 0.99, value = 0.8),
                              sliderInput('gamma3', 'Gamma:', min = -2, max = 2, value = 0.5),
                              sliderInput('datasets3', 'Datasets:', min = 1, max = 4, value = 2)
                            ),
                            mainPanel(
                              plotOutput('distPlot3')
                            )
                          )
                 ),
                 tabPanel('Dynamic Linear Model',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run4', 'Run'),
                              numericInput('M4', 'M:', min = 1, max = 500, value = 50),
                              numericInput('t4', 'T:', min = 50, max = 500, value = 200),
                              sliderInput('sigmasqy4', 'Sigma Squared (Y):', min = 0.01, max = 10, value = 1),
                              sliderInput('sigmasqx4', 'Sigma Squared (X):', min = 0.01, max = 10, value = 1),
                              sliderInput('phi4', 'Phi:', min = -0.99, max = 0.99, value = 0.8),
                              sliderInput('gamma4', 'Gamma:', min = -2, max = 2, value = 0.5),
                              sliderInput('datasets4', 'Datasets:', min = 1, max = 4, value = 2),
                              checkboxInput('meanfield4', 'Diagonal Covariance', value = TRUE),
                              checkboxInput('xderiv4', 'Fit Latent States', value = FALSE),
                              checkboxInput('var4', 'Fit Variance', value = FALSE)
                            ),
                            mainPanel(
                              plotOutput('distPlot4'),
                              plotOutput('distPlot4b')
                            )
                          )
                          
                 ),
                 tabPanel('Stochastic Volatility Model',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run5', 'Run'),
                              numericInput('M5', 'M:', min = 1, max = 500, value = 50),
                              numericInput('t5', 'T:', min = 50, max = 500, value = 200),
                              sliderInput('sigmasq5', 'Sigma Squared:', min = 0.01, max = 10, value = 1),
                              sliderInput('phi5', 'Phi:', min = -0.99, max = 0.99, value = 0.8),
                              sliderInput('gamma5', 'Gamma:', min = -2, max = 2, value = 0.5),
                              sliderInput('datasets5', 'Datasets:', min = 1, max = 4, value = 2),
                              checkboxInput('meanfield5', 'Diagonal Covariance', value = TRUE)
                            ),
                            mainPanel(
                              plotOutput('distPlot5'),
                              plotOutput('distPlot5b')
                            )
                          )
                 ),
                 tabPanel('MIVB',
                          sidebarLayout(
                            sidebarPanel(
                              actionButton('run6', 'Run'),
                              numericInput('M6', 'M:', min = 1, max = 500, value = 50),
                              selectInput('model6', 'Model:', choices=c('DLM', 'SVM')),
                              numericInput('t6', 'T:', min = 50, max = 1000, value = 200),
                              numericInput('s6', 'S:', min = 1, max = 1000, value = 50),
                              sliderInput('sigmasqy6', 'Sigma Squared (Y):', min = 0.01, max = 10, value = 1),
                              sliderInput('sigmasqx6', 'Sigma Squared (X):', min = 0.01, max = 10, value = 1),
                              sliderInput('phi6', 'Phi:', min = -0.99, max = 0.99, value = 0.8),
                              sliderInput('gamma6', 'Gamma:', min = -2, max = 2, value = 0.5),
                              sliderInput('datasets6', 'Datasets:', min = 1, max = 4, value = 2)
                            ),
                            mainPanel(
                              plotOutput('distPlot6'),
                              plotOutput('distPlot6b')
                            )
                          )
                          
                 )
                          
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  VBN <- reactive({
    if(input$run1 == 0){
      return()
    } 
    isolate({
      return(ADVI(input$t1, input$datasets1, input$sigmasq1, 'Normal', input$M1))
    })
  })
  
  VBAR1 <- reactive({
    if (input$run2 == 0){
      return()
    }
    isolate({
      return(ADVI(input$t2, input$datasets2, c(input$sigmasq2, input$phi2), 'AR1-Zero Mean', input$M2))
    })
  })
  
  
  VBAR1M <- reactive({
    if (input$run3 == 0) {
      return()
    }
    isolate({
      return(ADVI(input$t3, input$datasets3, c(input$sigmasq3, input$phi3, input$gamma3), 'AR1 with mean', input$M3))
    })
  })
  
  VBDLM <- reactive({
    if (input$run4 == 0) {
      return()
    }
    isolate({
      return(ADVI(input$t4, input$datasets4, 
                   c(input$sigmasqy4, input$sigmasqx4, input$phi4, input$gamma4), 'DLM', input$M4,
                   meanfield = input$meanfield4, xderiv = input$xderiv4, var = input$var4))
    })
  })
  
  VBSVM <- reactive({
    if (input$run5 == 0) {
      return()
    }
    isolate({
      return(ADVI(input$t5, input$datasets5,
                  c(input$sigmasq5, input$phi5, input$gamma5), 'SVM', input$M5,
                  meanfield = input$meanfield5))
    })
  })
  
  VBMIVB <- reactive({
    if (input$run6 == 0){
      return()
    } 
    isolate({
      return(MIVB(input$t6, input$s6, input$datasets6, 
                  c(input$sigmasqy6, input$sigmasqx6, input$phi6, input$gamma6),
                  input$model6, input$M6))
    })
  })
  
  output$distPlot1 <- renderPlot({
    VBfit = VBN()
    ggplot() + geom_line(data=VBfit$ADVI, aes(Support, Density, colour=Version)) + 
      geom_line(data=VBfit$True, aes(Support, Density)) +
      facet_grid(Variable~Dataset, scales='free') + theme_bw() +
      ggtitle(expression(paste(x[t],' = ',sigma,epsilon[t])))
    
  })
  
  output$distPlot2 <- renderPlot({
    VBfit = VBAR1()
    ggplot() + geom_line(data=VBfit$ADVI, aes(Support, Density, colour=Version)) + 
      geom_density(data=VBfit$True, aes(x=Draws)) +
      facet_grid(Variable~Dataset, scales='free', labeller=label_parsed) + theme_bw() + 
      theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
            axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank()) +
      ggtitle(expression(paste(x[t],' = ', phi, x[t-1], ' + ', sigma,epsilon[t])))
  })
  
  output$distPlot3 <- renderPlot({
    VBfit = VBAR1M()
    ggplot() + geom_line(data=VBfit$ADVI, aes(Support, Density, colour=Version)) + 
        geom_density(data=VBfit$True, aes(x=Draws)) +
        facet_grid(Variable~Dataset, scales='free', labeller=label_parsed) + theme_bw() + 
        theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
            axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank()) +
        ggtitle(expression(paste(x[t],' = ', gamma, ' + ', phi, '(', x[t-1], ' - ', gamma, ') + ', sigma,epsilon[t])))

  })
  
  output$distPlot4 <- renderPlot({
    VBfit = VBDLM()
    ggplot() + geom_line(data=VBfit$ADVI, aes(Support, Density, colour=Version)) + 
        geom_density(data=VBfit$True, aes(x=Draws)) +
        facet_grid(Variable~Dataset, scales='free', labeller=label_parsed) + theme_bw() + 
        theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
              axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank()) +
        ggtitle(expression(paste(y[t], ' = ', x[t], ' + ', sigma[y], nu[t], '       ', 
                               x[t],' = ', gamma, ' + ', phi, '(', x[t-1], ' - ', gamma, ') + ', sigma[x],epsilon[t])))
  })
  
  output$distPlot4b <- renderPlot({
    VBfit = VBDLM()
    ggplot(VBfit$Xt) + geom_line(aes(x=t, y=Actual)) + 
      geom_ribbon(aes(x=t, ymin=L95, ymax=U95, colour=Method), alpha=0) + theme_bw() + facet_grid(Dataset~.)
  })
  
  output$distplot5 <- renderPlot({
    VBfit = VBSVM()
    ggplot() + geom_line(data=VBfit$ADVI, aes(Support, Density, colour=Version)) + 
      geom_density(data=VBfit$True, aes(x=Draws)) +
      facet_grid(Variable~Dataset, scales='free', labeller=label_parsed) + theme_bw() + 
      theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
            axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank()) +
      ggtitle(expression(paste(y[t], ' = ', x[t], ' + log(', nu[t]^2, ')       ', 
                               x[t],' = ', gamma, ' + ', phi, '(', x[t-1], ' - ', gamma, ') + ', sigma[x],epsilon[t])))
  })
  
  output$distPlot5b <- renderPlot({
    VBfit = VBSVM()
    ggplot(VBfit$Xt) + geom_line(aes(x=t, y=Actual)) + 
      geom_ribbon(aes(x=t, ymin=L95, ymax=U95, colour=Method), alpha=0) + theme_bw() + facet_grid(Dataset~.)
  })
  
  output$distPlot6 <- renderPlot({
    VBfit = VBMIVB()
    ggplot() + geom_line(data=VBfit$VB, aes(Support, Density), colour='red') + 
      geom_density(data=VBfit$MCMC, aes(x=Draws)) +
      facet_grid(Variable~Dataset, scales='free', labeller=label_parsed) + theme_bw() + 
      theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
            axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank()) +
      ggtitle(expression(paste(y[t], ' = ', x[t], ' + ', sigma[y], nu[t], '       ', 
                               x[t],' = ', gamma, ' + ', phi, '(', x[t-1], ' - ', gamma, ') + ', sigma[x],epsilon[t])))
  })
  
  output$distPlot6b <- renderPlot({
    VBfit = VBMIVB()
    ggplot(VBfit$Xt) + geom_line(aes(x=t, y=Actual)) + 
      geom_ribbon(aes(x=t, ymin=L95, ymax=U95, colour=Method), alpha=0) + theme_bw() + facet_grid(Dataset~.)
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

