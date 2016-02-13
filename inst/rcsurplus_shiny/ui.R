##########################################################################
# Shiny user interface ###################################################
##########################################################################
shinyUI(navbarPage(
    title = 'RC SURPLUS',
    tabPanel("About",
             navbarPage('About',
                        tabPanel( 'Introduction', fluidRow(helpText( help$Introduction )) ),
                        tabPanel( 'Model', 
                                  fluidRow(
                                      withMathJax(''),
                                      helpText( help$Model1 ),
                                      helpText( help$Model2 ),
                                      helpText( help$Model3 ),
                                      helpText( help$Model4 ),
                                      helpText( help$Model5 ),
                                      helpText( help$Model6 )
                                  )),
                        tabPanel( 'UI', fluidRow(helpText( help$UI )) ),
                        tabPanel( 'Fitting', fluidRow(helpText( help$Fitting )) )
             )
    ),
    tabPanel("Input",
	    navbarPage('Input',
			tabPanel('Data',
                fluidPage(
                    sidebarLayout(
                        sidebarPanel(
                            checkboxInput("do_catch", "Plot catch", TRUE),
                            checkboxInput("do_effort", "Plot effort", FALSE),
                            checkboxInput("do_cpue", "Plot cpue", FALSE),
                            checkboxInput("one_row", "Plot in one row", TRUE)
                        ),
                        mainPanel(
                            plotOutput("plot_data")
                        )
                    )
                )
            )
		)
	),
    tabPanel('Model',
             fluidPage(
                 br(),
                 fluidRow(
                     column(4,
                            h4("Model(s) to fit (multi-choice)"),
                            selectizeInput("spm", "",choices = 
                                               c('Pella-Tomlinson', 'Schaefer', 'Fox', 'Alternative'),
                                           multiple=TRUE),
                            h4("Bounds of uniform priors"),
                            sliderInput("Kprior", "Catchability: \\(K\\) and \\(e^\\rho\\)",   
                                        min = myglobal$priorK[1], 
                                        max = myglobal$priorK[2], 
                                        value = c(myglobal$priorK[1],myglobal$priorK[2]), 
                                        step = 100 ),
                            sliderInput("rprior", "Growth rate: \\(r\\)",   
                                        min = myglobal$priorr[1], 
                                        max = myglobal$priorr[2], 
                                        value = c(myglobal$priorr[1],myglobal$priorr[2]), 
                                        step = 0.01 ),
                            sliderInput("PHIprior", "Elasticity: \\(\\phi\\)",   
                                        min = myglobal$priorPHI[1], 
                                        max = myglobal$priorPHI[2], 
                                        value = c(myglobal$priorPHI[1],myglobal$priorPHI[2]), 
                                        step = 0.01 ),
                            sliderInput("qprior", "Log-catchability: \\(\\log (q)\\) and \\(\\chi\\)",   
                                        min = myglobal$priorq[1], 
                                        max = myglobal$priorq[2], 
                                        value = c(myglobal$priorq[1],myglobal$priorq[2]), 
                                        step = 0.01 ),
                            sliderInput("sprior", "Model log-variance: \\(\\log (\\sigma)\\)",   
                                        min = myglobal$priors[1], 
                                        max = myglobal$priors[2], 
                                        value = c(myglobal$priors[1],myglobal$priors[2]), 
                                        step = 0.01 )
                     ),
                     column(4,
                            h4("MCMC settings"),
                            sliderInput("MCMCn", "Number of iterations for output:",   
                                        min = myglobal$mcmc_n[1], 
                                        max = myglobal$mcmc_n[2], 
                                        value = myglobal$mcmc_n[1], 
                                        step = 100 ),
                            sliderInput("MCMCb", "Burn-in ratio:",   
                                        min = myglobal$mcmc_b[1], 
                                        max = myglobal$mcmc_b[2], 
                                        value = myglobal$mcmc_b[2], 
                                        step = 0.05 ),
                            sliderInput("MCMCt", "Thinning factor:",   
                                        min = myglobal$mcmc_t[1], 
                                        max = myglobal$mcmc_t[2], 
                                        value = myglobal$mcmc_t[1], 
                                        step = 1 ),
                            sliderInput("MCMCc", "Number of chains:",   
                                        min = myglobal$mcmc_c[1], 
                                        max = myglobal$mcmc_c[2], 
                                        value = myglobal$mcmc_c[1], 
                                        step = 1 ),
                            checkboxInput("overrelax", "Over-relax", FALSE),
                            p('Total no. MCMC iterations:'),
                            verbatimTextOutput('model_settings')
                     ),
                     column(4,
                            h4("Convergence diagnostics"),
                            checkboxInput("do_coda", "CODA", FALSE),
                            h4("Model Status"),
                            bsAlert(inputId = "model_inputId"),
                            bsActionButton('fit_model', label = "Fit model(s)", style = "info", size="mini"),
                            verbatimTextOutput('model_status')
                     )
                 )
             )
    ),
    tabPanel('Output',
             navbarPage('Output',
                        tabPanel('Hyperparameters',
                                 fluidPage(
                                     sidebarLayout(
                                         sidebarPanel(
                                             h4('Plot:'),
                                             selectizeInput("spm_hyper", "",choices = c(
                                                 'Pella-Tomlinson',
                                                 'Schaefer',
                                                 'Fox',
                                                 'Alternative'), multiple=TRUE),
                                             checkboxInput("plot_KeRHO", "Carrying capacity: \\(K\\) and/or \\(e^\\rho\\)", FALSE),
                                             checkboxInput("plot_r", "Growth rate: \\(r\\)", FALSE),
                                             checkboxInput("plot_PHI", "Elasticity: \\(\\phi\\)", FALSE),
                                             checkboxInput("plot_qCHI", "Log-catchability: \\(\\log(q)\\) and/or \\(\\chi\\)", FALSE),
                                             checkboxInput("plot_lSIGMA", "Model log-variance: \\(\\log(\\sigma)\\)", FALSE),
                                             checkboxInput("plot_msy", "MSY", FALSE),
                                             checkboxInput("plot_bmsy", "B(MSY)", FALSE),
                                             checkboxInput("plot_fmsy", "F(MSY)", FALSE)
                                         ),
                                         mainPanel(
                                             plotOutput("plot_hyperparameters")
                                         )
                                     )
                                 )
                        ),
                        tabPanel('Fit',
                                 fluidPage(
                                     sidebarLayout(
                                         sidebarPanel(
                                             h4('Observed vs. fitted CPUE:'),
                                             selectizeInput("spm_fit", "",
                                                            choices = c('Pella-Tomlinson',
                                                                        'Schaefer',
                                                                         'Fox',
                                                                         'Alternative'),
                                                            multiple=TRUE),
                                             sliderInput("years_fit", "Years:",   
                                                         min = 1965, max = 1988, 
                                                         value = c(1965,1988), step = 1 ),
                                             checkboxInput("fit_one_row", "Plot in one row", TRUE)
                                         ),
                                         mainPanel(
                                             plotOutput("fit_plot")
                                         )
                                     )
                                 )
                        ),
                        tabPanel('Diagnostics',
                                 fluidPage(
                                     sidebarLayout(
                                         sidebarPanel(
                                             h4('Model:'),
                                             selectizeInput("spm_diag", "",choices = c('Pella-Tomlinson',
                                               'Schaefer',
                                               'Fox',
                                               'Alternative'), multiple=TRUE),
                                             h4('Diagnostic:'),
                                             selectInput("diag", "",choices = c('Summary & Deviance',
                                                                                'CODA - Gelman',
                                                                                'CODA - Geweke',
                                                                                'CODA - Raftery & Lewis'))
                                         ),
                                         mainPanel(
                                             verbatimTextOutput("summary_dev")
                                         )
                                     )
                                 )
                        )
             )
    )
))

