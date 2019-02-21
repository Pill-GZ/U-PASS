# user interface

ui <- fluidPage(
  titlePanel("GWAS power calculator"),
  
  fluidRow(
    
    column(3,
           br(),
           wellPanel("Sample sizes",
                     numericInput("n", "Totoal number of subjects (cases + controls):", 
                                  value = 10000, min = 500, max = 1000000, step = 100),
                     sliderInput("phi", "Fraction of cases:", value = 1/2, min = 0.01, max = 0.99, step = .01)
           ),
           br(),
           wellPanel("False discovery control",
                     selectInput("type_I_error_criteria", "Criteria for false discovery",
                                 list("Type I error", "False discovery rate (FDR)", "Family-wise error rate (FWER)")),
                     conditionalPanel(
                       condition = "input.type_I_error_criteria == 'Type I error'",
                       shinyWidgets::sliderTextInput(inputId = "alpha", 
                                                     label = "Target type I error rate:", 
                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                     selected = 0.05)
                     ),
                     conditionalPanel(
                       condition = "input.type_I_error_criteria == 'False discovery rate (FDR)'",
                       shinyWidgets::sliderTextInput(inputId = "alpha", 
                                                     label = "Target FDR:",
                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                     selected = 0.05)
                     ),
                     conditionalPanel(
                       condition = "input.type_I_error_criteria == 'Family-wise error rate (FWER)'",
                       shinyWidgets::sliderTextInput(inputId = "alpha", 
                                                     label = "Target FWER:", 
                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                     selected = 0.05)
                     )
           )
    ),
    
    column(9,
           # plotOutput("OR.RAF.plot"),
           plotlyOutput("OR.RAF.plotly")
    )
  )
)