##########################################################################
# Shiny server ###########################################################
##########################################################################
shinyServer(function(input, output, session) {
    
    output$plot_data  <- renderPlot(
        m$plot_data(input)
    )
  
    output$model_settings <- renderPrint({
        cat(floor(input$MCMCt * input$MCMCn / input$MCMCb ))
    })
    
    output$model_status <- renderPrint({
        if(any(input$Kprior>0) | any(input$rprior>0) | any(input$qprior>0) | 
           any(input$sprior>0) |
           input$MCMCn>0 | input$MCMCb>0 | input$MCMCt>0 | input$MCMCc>0 | 
           length(input$spm) == 0 | input$do_coda){
            closeAlert(session = session, "model_alertId")
            createAlert(session, anchorId = "model_inputId", 
                        alertId = "model_alertId", 
                        content = "Model(s) not fitted.", 
                        style = "danger", dismiss = FALSE, append = FALSE )
        }
        if(input$fit_model > m$fit_counter){
            closeAlert(session = session, "model_alertId")
            createAlert(session, anchorId = "model_inputId", 
                        alertId = "model_alertId", 
                        content = "Fitting model(s), please wait.", 
                        style = "warning", dismiss = FALSE, append = FALSE )
            ok <- m$fit_models(input)
            closeAlert(session = session, "model_alertId")
            createAlert(session, anchorId = "model_inputId", 
                        alertId = "model_alertId", 
                        content = "Model(s) fitted.", 
                        style = "success", dismiss = FALSE, append = TRUE)
        }
    })
    
    output$fit_plot  <- renderPlot(
        m$plot_fitted_cpue(input)
    )
    
    output$plot_hyperparameters <- renderPlot({
        m$plot_hyperparameters(input)
    })
    
    output$summary_dev <- renderPrint({
        sep <- c('*** Schaefer ****************************************\n',
                 '*** Alternative *************************************\n')
        if(input$diag == 'Summary & Deviance'){
            if(any(input$spm_diag == 'Pella-Tomlinson')){
                cat('*** Pella-Tomlinson *********************************\n')
                print(m$models[[1]])
                cat('\n')
            }
            if(any(input$spm_diag == 'Schaefer')){
                cat('*** Schaefer ****************************************\n')
                print(m$models[[2]])
                cat('\n')
            }
            if(any(input$spm_diag == 'Fox')){
                cat('*** Fox *********************************************\n')
                print(m$models[[3]])
                cat('\n')
            }
            if(any(input$spm_diag == 'Alternative')){
                cat('*** Alternative *************************************\n')
                print(m$models[[4]])
                cat('\n')
            }
        } else if(input$diag == 'CODA - Gelman' & input$do_coda){
            if(input$spm_diag != 'Alternative') {
                cat(sep[1])
                print(gelman.diag(m$coda_models[[1]]))
                cat('\n')
            }
            if(input$spm_diag != 'Schaefer'){
                cat(sep[2])
                print(gelman.diag(m$coda_models[[2]]))
                cat('\n')
            }
        } else if(input$diag == 'CODA - Geweke' & input$do_coda){
            if(input$spm_diag != 'Alternative') {
                cat(sep[1])
                print(geweke.diag(m$coda_models[[1]], frac1=0.1, frac2=0.5))
                cat('\n')
            }
            if(input$spm_diag != 'Schaefer'){
                cat(sep[2])
                print(geweke.diag(m$coda_models[[2]], frac1=0.1, frac2=0.5))
                cat('\n')
            }
        } else if(input$diag == 'CODA - Raftery & Lewis' & input$do_coda){
            if(input$spm_diag != 'Alternative') {
                cat(sep[1])
                print(raftery.diag(data=m$coda_models[[1]], 
                                   q=0.025, r=0.005, s=0.95, converge.eps=0.001))
                cat('\n')
            }
            if(input$spm_diag != 'Schaefer'){
                cat(sep[2])
                print(raftery.diag(data=m$coda_models[[2]], 
                                   q=0.025, r=0.005, s=0.95, converge.eps=0.001))
                cat('\n')
            }
        }
    })
  
})
