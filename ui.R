# user interface

ui <- navbarPage("GWAS power calculator",
                 
                 tabPanel("OR-RAF power diagram", id = "tabset",
                          tags$style(type = 'text/css', '.navbar {font-size: 15px;}', 
                                     '.navbar-default .navbar-brand {font-size: 30px;}'),
                          includeHTML("summary.html"),
                          fluidRow(
                            
                            column(3,
                                   br(),
                                   
                                   #### sample sizes ####
                                   
                                   wellPanel(selectInput("sample_size_specification", "Sample size specification",
                                                         list("Total number of subjects + fraction of controls", 
                                                              "Number of cases + number of controls")),
                                             conditionalPanel(
                                               condition = "input.sample_size_specification == 'Total number of subjects + fraction of controls'",
                                               numericInput("n", "Totoal number of subjects (cases + controls):", 
                                                            value = 20000, min = 500, max = 1000000, step = 100),
                                               sliderInput("phi", "Fraction of cases:", value = 1/2, min = 0.01, max = 0.99, step = .01)
                                             ),
                                             conditionalPanel(
                                               condition = "input.sample_size_specification == 'Number of cases + number of controls'",
                                               numericInput("n1", "Number of cases:", 
                                                            value = 10000, min = 250, max = 1000000, step = 50),
                                               numericInput("n2", "Number of controls:", 
                                                            value = 10000, min = 250, max = 1000000, step = 50)
                                             )
                                   ),
                                   # br(),
                                   
                                   #### false discovery control ####
                                   
                                   wellPanel(selectInput("type_I_error_criteria", "Criteria for false discovery",
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
                                               numericInput("p.FDR", "Effective dimension (number of loci):", 
                                                            value = 100000, min = 250, max = 10000000, step = 50),
                                               shinyWidgets::sliderTextInput(inputId = "alpha.FDR", 
                                                                             label = "Target FDR:",
                                                                             choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                             selected = 0.05)
                                             ),
                                             conditionalPanel(
                                               condition = "input.type_I_error_criteria == 'Family-wise error rate (FWER)'",
                                               numericInput("p.FWER", "Effective dimension (number of loci):", 
                                                            value = 100000, min = 250, max = 10000000, step = 50),
                                               shinyWidgets::sliderTextInput(inputId = "alpha.FWER", 
                                                                             label = "Target FWER:", 
                                                                             choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                             selected = 0.05)
                                             )
                                   ),
                                   # br(),
                                   
                                   #### datasets overlay ####
                                   
                                   wellPanel(checkboxInput("overlay_example_dataset", "Overlay data points", FALSE),
                                             conditionalPanel(
                                               condition = "input.overlay_example_dataset == true",
                                               selectInput("choose_dataset", "Example dataset from EBI",
                                                           list("Breast carcinoma", 
                                                                "Coronary heart disease", 
                                                                "Type II diabetes mellitus"))
                                             )
                                   )
                            ), # end of first column
                            
                            #### OR-RAF diagram display ####
                            
                            column(6, # "fixing height to avoid automatic adjustments"
                                   # textOutput("debug"),
                                   div(style = "height:1200px;", 
                                       plotlyOutput("OR.RAF.plotly", height = "700px"),
                                       br(),
                                       conditionalPanel(
                                         condition = "input.overlay_example_dataset == true",
                                         id = "gene_info",
                                         tags$style(type="text/css", '#gene_info { width:700px; }'),
                                         wellPanel(htmlOutput("selection"))
                                       ) # end of gene_info box
                                   ) # end of "fixing height to avoid automatic adjustments"
                            ) # end of second column
                            
                            #plotlyOutput("OR.RAF.heatmap.plotly")
                            #plotOutput("OR.RAF.plot")
                            
                          ), # end of fluidRow
                          includeHTML("credits.html")
                 ), # end of first panel (OR-RAF)
                 
           tabPanel("Power vs sample sizes",
                    textOutput("power.vec")),
           tabPanel("Design my study")
           
) # end of navbarPage
