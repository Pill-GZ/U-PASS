library("rintrojs")
# we strongly recommend using our version of the load spinner
# devtools::install_github("Pill-GZ/shinycssloaders")
library("shinycssloaders")

# user interface

ui <- shinyUI(tagList(
  tags$head(HTML("<title>U-PASS: a unified power analysis of association studies</title>"),
            tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
            tags$style(HTML(".shiny-output-error-validation { color: red; }")),
            tags$style(".rightAlign{float:right;}")),
  introjsUI(),
  navbarPage("U-PASS power calculator", id = "mainNavbarPage", theme = "bootstrap-cosmo-customized.css",
             
             #### OR-RAF tab ############################################################ ####
             tabPanel("OR-RAF power diagram", id = "OR-RAF_tab",
                      tags$style(type = 'text/css', '.navbar {font-size: 20px;}',
                                 '.navbar-default .navbar-brand {font-size: 30px;}'),
                      tags$style(HTML(".introjs-tooltip {max-width: 50%; min-width: 300px;}")),
                      includeHTML("www/header.html"),
                      HTML("<p>We encourage users to take a "), 
                      actionButton("help_ORRAF", "quick tour of the interface"),
                      HTML(", and check out the "),
                      actionButton(inputId="link_to_guide_from_tab1", label="User Guide"),
                      HTML("and detailed documentations in the help pages."),
                      
                      fluidRow(
                        
                        #### parameters and display options ####
                        
                        column(3,
                               br(),
                               #### sample sizes ####
                               
                               wellPanel(id = "sample_size",
                                         selectInput("sample_size_specification", "Sample size specification",
                                                     list("Number of subjects + fraction of cases", 
                                                          "Number of cases + number of controls")),
                                         conditionalPanel(
                                           condition = "input.sample_size_specification == 'Number of subjects + fraction of cases'",
                                           numericInput("n", "Totoal number of subjects:", 
                                                        value = 20000, min = 500, max = 1000000, step = 100),
                                           sliderInput("phi", "Fraction of cases:", value = 1/2, min = 0.05, max = 0.95, step = .05)
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
                                           numericInput("p_FDR", "Effective dimension (number of loci):", 
                                                        value = 1000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "alpha.FDR", 
                                                                         label = "Target FDR:",
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)
                                         ),
                                         conditionalPanel(
                                           condition = "input.type_I_error_criteria == 'Family-wise error rate (FWER)'",
                                           numericInput("p_FWER", "Effective dimension (number of loci):", 
                                                        value = 1000, min = 250, max = 10000000, step = 50),
                                           shinyWidgets::sliderTextInput(inputId = "alpha_FWER", 
                                                                         label = "Target FWER:", 
                                                                         choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                         selected = 0.05)
                                         ),
                                         htmlOutput("p_val_cutoff")
                               ),
                               # br(),
                               
                               #### OR-RAF diagram options ####
                               wellPanel(id = "plot_options", 
                                         tags$b("Plot options"),
                                         checkboxInput("overlay_equipower_curves", "Overlay equi-power curves", TRUE),
                                         checkboxInput("overlay_rare_variant_zones", "Overlay rare-variant zones", TRUE),
                                         conditionalPanel(
                                           condition = "input.overlay_rare_variant_zones == true",
                                           selectInput("rare_variant_zone_specification", "Rare variant is specified as...",
                                                       list("Minimum counts needed to calibrate Fisher's exact test",
                                                            "Absolute variant count in study",
                                                            "Fraction of total subjects in study")),
                                           conditionalPanel(
                                             condition = "input.rare_variant_zone_specification == 'Absolute variant count in study'",
                                             sliderInput(inputId = "rare_variant_threshold_count", 
                                                         label = "Variant counts need to be at least:",
                                                         value = 30, min = 5, max = 50, step = 5)
                                           ),
                                           conditionalPanel(
                                             condition = "input.rare_variant_zone_specification == 'Fraction of total subjects in study'",
                                             shinyWidgets::sliderTextInput(inputId = "rare_variant_threshold_fraction", 
                                                                           label = "Variant fraction need to be at least:",
                                                                           choices = paste0(as.vector(outer(c(1,2,5), 10^(-3:-2))) * 100, "%"),
                                                                           selected = "0.5%")
                                           )
                                         )
                               ),
                               # br(),
                               
                               #### datasets overlay ####
                               wellPanel(id = "data_options",
                                         tags$b("Data options"),
                                         checkboxInput("overlay_example_dataset", "Overlay reported findings from NHGRI-EBI GWAS Catalog", TRUE),
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
                        
                        column(8, 
                               # textOutput("debug"),
                               div(id = "OR-RAF_diagram",
                                   tags$style(type="text/css", '#OR-RAF_diagram { width:750px; }'),
                                   withSpinner(plotlyOutput("OR.RAF.plotly", height = "700px"))),
                               # plotOutput("OR.RAF.plot", height = "700px"), 
                               div(id = "gene_info",
                                   tags$style(type="text/css", '#gene_info { width:730px; }'),
                                   conditionalPanel(condition = "input.overlay_example_dataset == true",
                                                    wellPanel(htmlOutput("selection"))
                                   ) # end of gene_info box
                               )
                        ) # end of second column
                        
                        #plotlyOutput("OR.RAF.heatmap.plotly")
                        
                      ), # end of fluidRow
                      includeHTML("www/credits.html")
             ), # end of first rab panel (OR-RAF)
             
             #### Power analysis tab ####################################################### ####
             
             tabPanel("Design my study", value = "design_tab",
                      includeHTML("www/header.html"),
                      HTML("<p>We encourage users to take a "), 
                      actionButton("help_power_analysis", "quick tour of the interface"),
                      HTML(", and check out the "),
                      actionButton(inputId="link_to_guide_from_tab2", label="User Guide"),
                      HTML("and detailed documentations in the help pages."),
                      
                      fluidRow(
                        #### Study specifications ####
                        
                        column(3,
                               br(),
                               #### Step 1: model specifications ####
                               
                               introBox(
                                 wellPanel(id = "step1",
                                           selectInput("step1_model_specification", "Step 1: Model specifications",
                                                       list("Select a model specification",
                                                            "Allele frequency and odds ratio", 
                                                            "Disease model",
                                                            "Signal size per sample (advanced user)")),
                                           conditionalPanel(
                                             condition = "input.step1_model_specification == 'Allele frequency and odds ratio'",
                                             numericInput("target_RAF", "Target risk allele frequency (in control group):", 
                                                          value = 0.001, min = 0.001, max = 0.99999, step = 0.001),
                                             numericInput("target_OR", "Target odds ratio:", 
                                                          value = 5, min = 1.01, max = 1000000, step = 0.01)
                                           ),
                                           conditionalPanel(
                                             condition = "input.step1_model_specification == 'Disease model'",
                                             selectInput("target_disease_model", "Disease model:",
                                                         list("Multiplicative", "Additive", "Dominant", "Recessive")),
                                             numericInput("target_prevalence", "Disease prevalence in the general population:", 
                                                          value = 0.1, min = 0.01, max = 0.99, step = 0.01),
                                             numericInput("target_RAF_population", "Risk allele frequency in the general population:", 
                                                          value = 0.3, min = 0.01, max = 0.99, step = 0.01),
                                             numericInput("target_GRR", "Genotype relative risk:", 
                                                          value = 1.2, min = 1.01, max = 100, step = 0.01)
                                           ),
                                           conditionalPanel(
                                             condition = "input.step1_model_specification == 'Signal size per sample (advanced user)'",
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
                                 
                                 condition = "input.step1_model_specification != 'Select a model specification'",
                                 wellPanel(id = "step2",
                                   selectInput("step2_fixed_quantity", "Step 2: Sample size constraints",
                                               list("Select a constraint",
                                                    "Budget / total number of subjects", 
                                                    "Number of cases",
                                                    "Fraction of cases")),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Budget / total number of subjects'",
                                     numericInput("fixed_budget", "Total subjects (cases + controls):", 
                                                  value = 20000, min = 500, max = 1000000, step = 100)
                                   ),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Number of cases'",
                                     numericInput("fixed_cases", "Number of cases in study:", 
                                                  value = 10000, min = 500, max = 1000000, step = 100)
                                   ),
                                   conditionalPanel(
                                     condition = "input.step2_fixed_quantity == 'Fraction of cases'",
                                     sliderInput("fixed_phi", "Fraction of cases in the study:", value = 1/2, min = 0.05, max = 0.95, step = .05)
                                   )
                                 ), # end of step 2: select a constraint
                                 
                                 conditionalPanel( 
                                   
                                   #### Step 3: choosing false discovery control criteria ####
                                   # only displayed if step 1 and 2 are complete
                                   
                                   condition = "input.step2_fixed_quantity != 'Select a constraint'",
                                   wellPanel(id = "step3",
                                     selectInput("step3_type_I_error_criteria", "Step 3: Criteria for false discovery",
                                                 list("Select a criterion", "Family-wise error rate (FWER)", 
                                                      #"False discovery rate (FDR)", 
                                                      "Type I error")),
                                     conditionalPanel(
                                       condition = "input.step3_type_I_error_criteria == 'Type I error'",
                                       shinyWidgets::sliderTextInput(inputId = "design_alpha", 
                                                                     label = "Target type I error rate:", 
                                                                     choices = c(as.vector(outer(c(1,5), 10^(-4:-2))),0.1),
                                                                     selected = 0.05)),
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
                                                                     selected = 0.05)),
                                     conditionalPanel(
                                       condition = "input.step3_type_I_error_criteria != 'Select a criterion'",
                                       htmlOutput("design_p_val_cutoff")
                                     )
                                   ), # end of step 3: select a false discovery control
                                   
                                   conditionalPanel( 
                                     
                                     #### Step 4: choosing non-discovery control target ####
                                     # only displayed if step 1, 2, and 3 are complete
                                     
                                     condition = "input.step3_type_I_error_criteria != 'Select a criterion'",
                                     wellPanel(id = "step4",
                                       selectInput("step4_type_II_error_criteria", "Step 4: Target for non-discovery",
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
                        ) # end of second column
                        
                        
                      ), # end of fliudRow
                      
                      includeHTML("www/credits.html")
             ), # end of second tab panel (study design)
             
             
             #### Disease model converter tab ################################################# ####
             
             tabPanel(HTML("Disease model converter"), value = "disease_model_converter",
                      
                      includeHTML("www/header.html"),
                      HTML("<p>We encourage users to take a "), 
                      actionButton("help_model_converter", "quick tour of the interface"),
                      HTML(", and check out the "),
                      actionButton(inputId="link_to_guide_from_tab3", label="User Guide"),
                      HTML("and detailed documentations in the help pages."),
                      
                      fluidRow(
                        
                        br(),
                        #### Disease model specifications ####
                        
                        column(3, offset = 0,
                               wellPanel(id = "disease_model_specification",
                                         selectInput("disease_model", "Disease model:",
                                                     list("Multiplicative", "Additive", "Dominant", "Recessive")),
                                         numericInput("disease_model_prevalence", "Disease prevalence in the population:", 
                                                      value = 0.1, min = 0, max = 0.99, step = 0.01),
                                         numericInput("disease_model_RAF_population", "Risk allele frequency in the population:", 
                                                      value = 0.3, min = 0, max = 0.99, step = 1e-5),
                                         numericInput("disease_model_GRR", "Genotype relative risk:", 
                                                      value = 1.5, min = 1.01, max = 100, step = 0.01),
                                         actionButton(inputId = "use_disease_model_specification",  class = 'rightAlign', 
                                                      label = "Go to power calculator")
                               ) # end of disease model specifications
                        ), # end of left column
                        
                        #### disease model conversions ####
                        
                        column(1, 
                               br(), br(), br(), br(), br(),
                               HTML('<center> <img src="arrow.png" width="100%" /> </center>')),
                        column(3, offset = 0,
                               br(), br(), br(),
                               wellPanel(id = "disease_model_conversion_results",
                                 HTML("<b>Risk allele frequency in controls:</b>"),
                                 panel(
                                   htmlOutput("disease_model_converter_result_f")
                                 ),
                                 HTML("<b>Odds ratio between allele variants:</b>"),
                                 panel(
                                   htmlOutput("disease_model_converter_result_R")
                                 ),
                                 actionButton(inputId = "use_canonical_specification",  class = 'rightAlign', 
                                              label = htmlOutput("disease_model_converter_message"))
                               ) # end of conversion output
                        ) # end of right column
                      ), # enf of disease model converter
                      
                      br(),
                      
                      fluidRow(
                        column(7, offset = 0,
                               includeHTML("www/disease_model_converter.html"),
                               HTML("<p>A few things to note:<ul><li>Disease model parameters may be incompatible. 
                                    Try, e.g., "), 
                               actionButton("disease_model_preset_incompatible", "an incompatible disease model"),
                               HTML(".</li><li>Multiple disease models can map to the same set of canonical parameters. 
                                    Try, e.g., these"), 
                               actionButton("disease_model_preset_multiplicative", "Multiplicative"),
                               actionButton("disease_model_preset_additive", "Additive"),
                               actionButton("disease_model_preset_dominant", "Dominant"),
                               actionButton("disease_model_preset_recessive", "Recessive"),
                               HTML("models.</li></ul>"), 
                               br(), br()
                        ) 
                      ), # end of comments to disease model conversions
                      includeHTML("www/credits.html")
             ), # end of disease model converter
             
             #### Help pages ####
             
             navbarMenu("Help",
                        # Documentation page is inserted below the User Guide
                        tabPanel(HTML('User Guide</a></li>
                                      <li><a href=\"disease_models_revisited.html\" target=\"_blank\">Disease Models Revisited 
                                      <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAQElEQVR42qXKwQkAIAxDUUdxtO6/RBQkQZvSi8I/pL4BoGw/XPkh4XigPmsUgh0626AjRsgxHTkUThsG2T/sIlzdTsp52kSS1wAAAABJRU5ErkJggg=="></a></li>
                                      <li><a href=\"unified_power_analysis.html\" target=\"_blank\">Unified Power Analysis 
                                      <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAQElEQVR42qXKwQkAIAxDUUdxtO6/RBQkQZvSi8I/pL4BoGw/XPkh4XigPmsUgh0626AjRsgxHTkUThsG2T/sIlzdTsp52kSS1wAAAABJRU5ErkJggg==">'), 
                                 value = "user_guide",
                                 fluidRow(
                                   column(6, offset = 3,
                                          withMathJax(includeHTML("www/user_guide.html")),
                                          includeHTML("www/credits.html")
                                   )
                                 )
                        )
             ), # end of help menu bar 
             
             #### About pages ####
             
             navbarMenu("About",
                        tabPanel(HTML("Download and Installation"), value = "download",
                                 fluidRow(
                                   column(6, offset = 3,
                                          includeHTML("www/download_and_installation.html")
                                   )
                                 )
                        ),
                        tabPanel(HTML("Citation and Contact"), value = "contact",
                                 fluidRow(
                                   column(6, offset = 3,
                                          includeHTML("www/contact.html")
                                   )
                                 )
                        )
             ) # end of about menu bar
             
             #### page ends ####
             ) # end of navbarPage
  ) # end of tagList
) # end of shinyUI
