# user interface

ui <- shinyUI(tagList(
  introjsUI(),
  navbarPage("GWAS power calculator",
             
             #### OR-RAF tab ####
             
             tabPanel("OR-RAF power diagram", id = "OR-RAF_tab",
                      tags$style(type = 'text/css', '.navbar {font-size: 15px;}', 
                                 '.navbar-default .navbar-brand {font-size: 30px;}'),
                      tags$style(HTML(".introjs-tooltip {max-width: 50%; min-width: 300px;}")),
                      includeHTML("header.html"),
                      actionButton("help_ORRAF", "Take a quick tour of the interface"),
                      
                      fluidRow(
                        
                        #### parameters and display options ####
                        
                        column(3,
                               br(),
                               #### sample sizes ####
                               
                               wellPanel(id = "sample_size",
                                         selectInput("sample_size_specification", "Sample size specification",
                                                     list("Number of subjects + fraction of Cases", 
                                                          "Number of Cases + number of Controls")),
                                         conditionalPanel(
                                           condition = "input.sample_size_specification == 'Number of subjects + fraction of Cases'",
                                           numericInput("n", "Totoal number of subjects:", 
                                                        value = 20000, min = 500, max = 1000000, step = 100),
                                           sliderInput("phi", "Fraction of cases:", value = 1/2, min = 0.01, max = 0.99, step = .01)
                                         ),
                                         conditionalPanel(
                                           condition = "input.sample_size_specification == 'Number of Cases + number of Controls'",
                                           numericInput("n1", "Number of cases:", 
                                                        value = 10000, min = 250, max = 1000000, step = 50),
                                           numericInput("n2", "Number of controls:", 
                                                        value = 10000, min = 250, max = 1000000, step = 50)
                                         )
                               ),
                               # br(),
                               
                               #### false discovery control ####
                               wellPanel(id = "false_discovery_control",
                                         selectInput("type_I_error_criteria", "Criteria for false discovery",
                                                     list("Family-wise error rate (FWER)", 
                                                          #"False discovery rate (FDR)", 
                                                          "Type I error")),
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
                                                        value = 1000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "alpha.FDR", 
                                                                         label = "Target FDR:",
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)
                                         ),
                                         conditionalPanel(
                                           condition = "input.type_I_error_criteria == 'Family-wise error rate (FWER)'",
                                           numericInput("p.FWER", "Effective dimension (number of loci):", 
                                                        value = 1000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "alpha.FWER", 
                                                                         label = "Target FWER:", 
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)
                                         ),
                                         htmlOutput("p_val_cutoff")
                               ),
                               # br(),
                               
                               #### datasets overlay ####
                               wellPanel(id = "overlay_dataset",
                                         checkboxInput("overlay_example_dataset", "Overlay reported findings from NHGRI-EBI GWAS Catalog", FALSE),
                                         conditionalPanel(
                                           condition = "input.overlay_example_dataset == true",
                                           selectInput("choose_dataset", "Choose a dataset", 
                                                       list("Example dataset from the NHGRI-EBI GWAS Catalog" = 
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
                                             tags$a("NHGRI-EBI GWAS Catalog format", href="https://www.ebi.ac.uk/gwas/docs/fileheaders", target="_blank"),
                                             ".",
                                             fileInput("my_data_upload", "Choose TSV File",
                                                       multiple = FALSE,
                                                       accept = c(".tsv" # "text/csv", "text/comma-separated-values,text/plain",
                                                       ))
                                           ), # end of upload my own data
                                           radioButtons("adaptive_sample_size", 
                                                        label = "Adaptive initial sample sizes and false discovery control according to study of selelected loci (experimental)",
                                                        choices = list("Enabled (best guesses from reported text)" = TRUE, 
                                                                       "Disabled (mannually override in boxes above)" = FALSE), 
                                                        selected = TRUE)
                                         ) # end of conditional panel
                               ) # end of dataset overlay choices
                               
                        ), # end of first column
                        
                        #### display OR-RAF diagram ####
                        
                        column(8, # "fixing height to avoid automatic adjustments"
                               #textOutput("debug"),
                               div(style = "height:1200px;", 
                                   div(id = "OR-RAF_diagram",
                                       tags$style(type="text/css", '#OR-RAF_diagram { width:750px; }'),
                                       plotlyOutput("OR.RAF.plotly", height = "700px")),
                                   # plotOutput("OR.RAF.plot", height = "700px"), 
                                   div(id = "gene_info",
                                       tags$style(type="text/css", '#gene_info { width:730px; }'),
                                       conditionalPanel(condition = "input.overlay_example_dataset == true",
                                                        wellPanel(htmlOutput("selection"))
                                       ) # end of gene_info box
                                   )
                               ) # end of "fixing height to avoid automatic adjustments"
                        ) # end of second column
                        
                        #plotlyOutput("OR.RAF.heatmap.plotly")
                        
                      ), # end of fluidRow
                      includeHTML("credits.html")
             ), # end of first rab panel (OR-RAF)
             
             #### Power analysis tab ####
             
             tabPanel("Design my study", id = "design_tab",
                      includeHTML("header.html"),
                      actionButton("help_power_analysis", "Take a quick tour of the interface"),
                      
                      fluidRow(
                        #### Study specifications ####
                        
                        column(3,
                               br(),
                               #### Step 1: target RAF and OR ####
                               
                               introBox(
                                 wellPanel(id = "step1",
                                           selectInput("step1_target_OR_RAF", "Step 1: I want to target a specific ...",
                                                       list("Select a target",
                                                            "Allele frequency and odds ratio", 
                                                            "Signal size per sample (advanced user)")),
                                           conditionalPanel(
                                             condition = "input.step1_target_OR_RAF == 'Allele frequency and odds ratio'",
                                             numericInput("target_RAF", "Target risk allele frequency:", 
                                                          value = 0.001, min = 0.001, max = 0.99999, step = 0.001),
                                             numericInput("target_OR", "Target odds ratio:", 
                                                          value = 5, min = 1.01, max = 1000000, step = 0.01)
                                           ),
                                           conditionalPanel(
                                             condition = "input.step1_target_OR_RAF == 'Signal size per sample (advanced user)'",
                                             shinyWidgets::sliderTextInput(inputId = "target_w2", 
                                                                           label = "Target signal size:",
                                                                           choices = c(as.vector(outer(c(1,2,5), 10^(-6:-2))),0.1),
                                                                           selected = 0.0001)
                                           )
                                 ), # end of step 1: select target 
                                 data.step = 6,
                                 data.intro = "Second tab"
                               ),
                               
                               conditionalPanel( 
                                 
                                 #### Step 2: constraint or fixed quantity in the study ####
                                 # only displayed if step 1 is complete
                                 
                                 condition = "input.step1_target_OR_RAF != 'Select a target'",
                                 wellPanel(id = "step2",
                                   selectInput("step2_fixed_quantity", "Step 2: I have a fixed ...",
                                               list("Select a constraint",
                                                    "Budget / total number of subjects", 
                                                    "Number of Cases",
                                                    "Fraction of Cases")),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Budget / total number of subjects'",
                                     numericInput("fixed_budget", "Total subjects (Cases + Controls):", 
                                                  value = 40000, min = 500, max = 1000000, step = 100)
                                   ),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Number of Cases'",
                                     numericInput("fixed_cases", "Number of Cases in study:", 
                                                  value = 20000, min = 500, max = 1000000, step = 100)
                                   ),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Fraction of Cases'",
                                     sliderInput("fixed_phi", "Fraction of cases in the study:", value = 1/2, min = 0.01, max = 0.99, step = .01)
                                   )
                                 ), # end of step 2: select a constraint
                                 
                                 conditionalPanel( 
                                   
                                   #### Step 3: choosing false discovery control criteria ####
                                   # only displayed if step 1 and 2 are complete
                                   
                                   condition = "input.step2_fixed_quantity != 'Select a constraint'",
                                   wellPanel(id = "step3",
                                     selectInput("step3_type_I_error_criteria", "Step 3: Criteria for false discovery ...",
                                                 list("Select a criterion", "Family-wise error rate (FWER)", 
                                                      #"False discovery rate (FDR)", 
                                                      "Type I error")),
                                     conditionalPanel(
                                       condition = "input.step3_type_I_error_criteria == 'Type I error'",
                                       shinyWidgets::sliderTextInput(inputId = "design_alpha", 
                                                                     label = "Target type I error rate:", 
                                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                     selected = 0.05)
                                     ),
                                     conditionalPanel(
                                       condition = "input.step3_type_I_error_criteria == 'False discovery rate (FDR)'",
                                       numericInput("design_p_FDR", "Effective dimension (number of loci):", 
                                                    value = 100000, min = 250, max = 10000000, step = 50),
                                       shinyWidgets::sliderTextInput(inputId = "design_alpha_FDR", 
                                                                     label = "Target FDR:",
                                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                     selected = 0.05)),
                                     conditionalPanel(
                                       condition = "input.step3_type_I_error_criteria == 'Family-wise error rate (FWER)'",
                                       numericInput("design_p_FWER", "Effective dimension (number of loci):", 
                                                    value = 100000, min = 250, max = 10000000, step = 50),
                                       shinyWidgets::sliderTextInput(inputId = "design_alpha_FWER", 
                                                                     label = "Target FWER:", 
                                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                     selected = 0.05))
                                   ), # end of step 3: select a false discovery control
                                   
                                   conditionalPanel( 
                                     
                                     #### Step 4: choosing non-discovery control target ####
                                     # only displayed if step 1, 2, and 3 are complete
                                     
                                     condition = "input.step3_type_I_error_criteria != 'Select a criterion'",
                                     wellPanel(id = "step4",
                                       selectInput("step4_type_II_error_criteria", "Step 4: Targert for non-discovery ...",
                                                   list("Select a criterion", 
                                                        "Type II error / non-discovery proportion (NDP)", 
                                                        "Family-wise non-discovery rate (FWNDR)")),
                                       conditionalPanel(
                                         condition = "input.step4_type_II_error_criteria == 'Type II error / non-discovery proportion (NDP)'",
                                         shinyWidgets::sliderTextInput(inputId = "design_power", 
                                                                       label = "Target type II error rate / (1-power) / NDP:", 
                                                                       choices = as.vector(outer(c(1,2,5), 10^(-3:-1))),
                                                                       selected = 0.2)
                                       ),
                                       conditionalPanel(
                                         condition = "input.step4_type_II_error_criteria == 'Family-wise non-discovery rate (FWNDR)'",
                                         numericInput("design_sparsity", "Sparsity / number of loci with equal signal or stronger:", 
                                                      value = 100, min = 1, max = 100000, step = 1),
                                         shinyWidgets::sliderTextInput(inputId = "design_FWNDR", 
                                                                       label = "Target FWNDR:",
                                                                       choices = as.vector(outer(c(1,2,5), 10^(-3:-1))),
                                                                       selected = 0.2))
                                     ) # end of step 4: select a non-discovery (type II error) goal
                                     
                                     
                                   ) # end of panels conditional on step 3
                                 ) # end of panels conditional on step 2
                               ) # end of panels conditional on step 1
                               
                        ), # end of first column
                        
                        #### display results from power analysis ####
                        
                        column(9, # "fixing height to avoid automatic adjustments"
                               #textOutput("debug"),
                               div(style = "height:1200px;", 
                                   #textOutput("waiting_for_design2"),
                                   #textOutput("power_vec"),
                                   #textOutput("waiting_for_design"),
                                   #textOutput("debug"),
                                   div(id = "power_analysis_results",
                                       tags$style(type="text/css", '#power_analysis_results { width:700px; }'),
                                       conditionalPanel(condition = "output.waiting_for_design == true",
                                                        br(),
                                                        "Complete your study design on the left"),
                                       conditionalPanel(condition = "output.fixed_n_design == true",
                                                        plotlyOutput("optimal_design_fixed_n", height = "700px")),
                                       conditionalPanel(condition = "output.fixed_n1_design == true",
                                                        plotlyOutput("optimal_design_fixed_n1", height = "700px")),
                                       conditionalPanel(condition = "output.fixed_phi_design == true",
                                                        plotlyOutput("optimal_design_fixed_phi", height = "700px"))
                                       #textOutput("waiting_for_design")
                                   )
                               ) # end of "fixing height to avoid automatic adjustments"
                        ) # end of second column
                        
                        
                      ), # end of fliudRow
                      
                      includeHTML("credits.html")
             ) # end of second tab panel (study design)
             
             #### page ends ####
             ) # end of navbarPage
  ) # end of tagList
) # end of shinyUI
