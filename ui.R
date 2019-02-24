# user interface

ui <- navbarPage("GWAS power calculator",
                 
                 #### OR-RAF tab ####
                 
                 tabPanel("OR-RAF power diagram", id = "tabset",
                          tags$style(type = 'text/css', '.navbar {font-size: 15px;}', 
                                     '.navbar-default .navbar-brand {font-size: 30px;}'),
                          includeHTML("header.html"),
                          fluidRow(
                            
                            #### parameters and display options ####
                            
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
                                               selectInput("choose_dataset", "Choose a dataset", 
                                                           list("Example dataset from European Bioinformatics Institute (EBI)" = 
                                                                  list("Breast carcinoma", 
                                                                       "Coronary heart disease", 
                                                                       "Type II diabetes mellitus"),
                                                                "User upload" = list("Upload my own data")
                                                                ) # end of drop-down menu options
                                                           ), # end of drop-down menu to select input
                                               conditionalPanel(
                                                 condition = "input.choose_dataset == 'Upload my own data'", 
                                                 # Input: Select a file 
                                                 "Make sure you follow the ",
                                                 tags$a("EBI format", href="https://www.ebi.ac.uk/gwas/docs/fileheaders", target="_blank"),
                                                 ".",
                                                 fileInput("my_data_upload", "Choose TSV File",
                                                           multiple = FALSE,
                                                           accept = c(".tsv" # "text/csv", "text/comma-separated-values,text/plain",
                                                                      ))
                                                 )
                                             ) # end of conditional panel
                                   ) # end of dataset overlay choices
                                   
                            ), # end of first column
                            
                            #### display OR-RAF diagram ####
                            
                            column(6, # "fixing height to avoid automatic adjustments"
                                   # textOutput("debug"),
                                   div(style = "height:1200px;", 
                                       plotlyOutput("OR.RAF.plotly", height = "700px"),
                                       # plotOutput("OR.RAF.plot", height = "700px"), 
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
                            
                          ), # end of fluidRow
                          includeHTML("credits.html")
                 ), # end of first rab panel (OR-RAF)
                 
                 #### Power analysis tab ####
                 
                 tabPanel("Design my study",
                          includeHTML("header.html"),
                          fluidRow(
                            #### Study specifications ####
                            
                            column(3,
                                   br(),
                                   
                                   #### Step 1: target RAF and OR ####
                                   
                                   wellPanel(selectInput("step1_target_OR_RAF", "Step 1: I want to target ...",
                                                         list("Select a target",
                                                              "A specific allele frequency and odds ratio", 
                                                              "A specific signal size (advanced user)")),
                                             conditionalPanel(
                                               condition = "input.step1_target_OR_RAF == 'A specific allele frequency and odds ratio'",
                                               numericInput("target_RAF", "Target risk allele frequency:", 
                                                            value = 0.1, min = 0.001, max = 0.99999, step = 0.001),
                                               numericInput("target_OR", "Target odds ratio:", 
                                                            value = 1.5, min = 1.01, max = 1000000, step = 0.01)
                                             ),
                                             conditionalPanel(
                                               condition = "input.step1_target_OR_RAF == 'A specific signal size (advanced user)'",
                                               shinyWidgets::sliderTextInput(inputId = "target_w2", 
                                                                             label = "Target signal size:",
                                                                             choices = c(as.vector(outer(c(1,2,5), 10^(-6:-2))),0.1),
                                                                             selected = 0.0001)
                                             )
                                   ), # end of step 1: select target 
                                   
                                   
                                   conditionalPanel( 
                                     
                                     #### Step 2: constraint or fixed quantity in the study ####
                                     # only displayed if step 1 is complete
                                     
                                     condition = "input.step1_target_OR_RAF != 'Select a target'",
                                     wellPanel(
                                       selectInput("step2_fixed_quantity", "Step 2: I have a fixed ...",
                                                   list("Select a contraint",
                                                        "Budget, i.e., total number of subjects", 
                                                        "Number of subjects in the Case group",
                                                        "Fraction of Cases among all subjects")),
                                       conditionalPanel(
                                         condition = "input.step2_fixed_quantity == 'Budget, i.e., total number of subjects'",
                                         numericInput("fixed_budget", "Total subjects (Cases + Controls) in study:", 
                                                      value = 20000, min = 500, max = 1000000, step = 100)
                                       ),
                                       conditionalPanel(
                                         condition = "input.step2_fixed_quantity == 'Number of subjects in the Case group'",
                                         numericInput("fixed_cases", "Number of Cases in study:", 
                                                      value = 10000, min = 500, max = 1000000, step = 100)
                                       ),
                                       conditionalPanel(
                                         condition = "input.step2_fixed_quantity == 'Fraction of Cases among all subjects'",
                                         sliderInput("fixed_phi", "Fraction of cases in the study:", value = 1/2, min = 0.01, max = 0.99, step = .01)
                                       )
                                     ), # end of step 2: select a constraint
                                     
                                     conditionalPanel( 
                                       
                                       #### Step 3: choosing false discovery control criteria ####
                                       # only displayed if step 1 and 2 are complete
                                       
                                       condition = "input.step2_fixed_quantity != 'Select a contraint'",
                                       wellPanel(
                                         selectInput("step3_type_I_error_criteria", "Step 3: Criteria for false discovery is ...",
                                                     list("Select a criteria", "Type I error", "False discovery rate (FDR)", "Family-wise error rate (FWER)")),
                                         conditionalPanel(
                                           condition = "input.step3_type_I_error_criteria == 'Type I error'",
                                           shinyWidgets::sliderTextInput(inputId = "design_alpha", 
                                                                         label = "Target type I error rate:", 
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)
                                         ),
                                         conditionalPanel(
                                           condition = "input.step3_type_I_error_criteria == 'False discovery rate (FDR)'",
                                           numericInput("p.FDR", "Effective dimension (number of loci):", 
                                                        value = 100000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "design_alpha.FDR", 
                                                                         label = "Target FDR:",
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)),
                                         conditionalPanel(
                                           condition = "input.step3_type_I_error_criteria == 'Family-wise error rate (FWER)'",
                                           numericInput("p.FWER", "Effective dimension (number of loci):", 
                                                        value = 100000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "design_alpha.FWER", 
                                                                         label = "Target FWER:", 
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05))
                                       ) # end of step 3: select a false discovery control
                                       
                                     ) # end of panels conditional on step 2
                                   ) # end of panels conditional on step 1
                                   
                            ), # end of first column
                            
                            #### display results from power analysis ####
                            
                            column(6, # "fixing height to avoid automatic adjustments"
                                   # textOutput("debug"),
                                   div(style = "height:1200px;", 
                                       textOutput("power.vec")
                                   ) # end of "fixing height to avoid automatic adjustments"
                            ) # end of second column
                          
                          
                          ), # end of fliudRow
                          
                          includeHTML("credits.html")
                 ) # end of second tab panel (study design)
                 
                 #### page ends ####
           
) # end of navbarPage
