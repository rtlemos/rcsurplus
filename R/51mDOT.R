nereus <- setRefClass(
    Class = 'nereus',
    fields = list(
        fit_counter = 'numeric',
        models = 'list',
        coda_models = 'list',
        model_names = 'character',
        par_titles = 'data.frame',
        IO = '.IAT'
    ),
    methods = list(
        initialize = function(){
            .self$fit_counter <- 0
            .self$models <- vector('list',length=2)
            .self$IO <- .IOT()
            .self$model_names <- c('Pella-Tomlinson', 'Schaefer', 'Fox', 'Alternative')
            .self$par_titles <- data.frame(
                pella_tomlinson=c('K', 'r', 'phi', 'log(q)',
                   'log(sigma)','MSY','B(MSY)','F(MSY)'),
                schaefer=c('K', 'r', 'phi', 'log(q)',
                   'log(sigma)','MSY','B(MSY)','F(MSY)'),
                fox=c('K', 'r', 'phi', 'log(q)',
                    'log(sigma)','MSY','B(MSY)','F(MSY)'),
                alternative = c('exp(rho)', 'r', 'phi', 'chi',
                   'log(sigma)','MSY','B(MSY)','F(MSY)'),
                stringsAsFactors = FALSE)
        },
        fit_models = function(input){
            'Function to fit surplus production models'
            
            if(any(input$spm =='Pella-Tomlinson')){
                .self$get_pella_tomlinson(input)
            } else {
                .self$models[[1]] <- NA
            }
            if(any(input$spm =='Schaefer')){
                .self$get_schaefer(input)
            } else {
                .self$models[[2]] <- NA
            }
            if(any(input$spm =='Fox')){
                .self$get_fox(input)
            } else {
                .self$models[[3]] <- NA
            }
            if(any(input$spm =='Alternative')){
                .self$get_alternative(input)
            } else {
                .self$models[[4]] <- NA
            }
            .self$fit_counter <- .self$fit_counter + 1
        },
        
        get_pella_tomlinson = function(input){
            model_file <- file.path(tempdir(), "pella_tomlinson_model.txt") 
            nobs    <- myglobal$nobs
            Y       <- hake$log_cpue
            C       <- hake$catch
            Ef      <- hake$effort
            lsprior <- input$sprior
            rprior  <- input$rprior
            PHIprior<- input$PHIprior
            Kprior  <- input$Kprior
            lqprior <- input$qprior
            mymodel <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                PHI ~ dunif(PHIprior[1],PHIprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs-1)) {
                    meanP[t+1] <- max( P[t]+(r/PHI)*P[t]*( 1-pow(P[t],PHI) )-C[t]/K, 0.001)
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY  <- r * K * exp( -1*(1/PHI + 1)*log(1+PHI) ) 
                BMSY <- K * exp( -1/PHI *log(1+PHI) ) 
                FMSY <- r / (1+PHI) 
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior","rprior","PHIprior", "Kprior","lqprior")
            params <- c("K", "r", "PHI", "lq", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() { list(  r = mean(input$rprior),
                                          PHI = mean(input$PHIprior),
                                          lq = mean(input$qprior), 
                                          K = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1-input$MCMCb)) )
            .self$models[[1]] <- bugs(data=data, inits=inits,  parameters.to.save=params,  
                                      model.file = model_file, n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE, over.relax = input$overrelax )
            if(input$do_coda){
                outS <- bugs(data=data, inits=inits, 
                             parameters.to.save=c("K", "r", "PHI", "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE, over.relax = input$overrelax )
                .self$coda_models[[1]] <- read.bugs(outS)
            }
        },
        
        get_schaefer = function(input){
            model_file <- file.path(tempdir(), "schaefer_model.txt") 
            nobs    <- myglobal$nobs
            Y       <- hake$log_cpue
            C       <- hake$catch
            Ef      <- hake$effort
            lsprior <- input$sprior
            rprior  <- input$rprior
            Kprior  <- input$Kprior
            lqprior <- input$qprior
            mymodel <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs-1)) {
                    meanP[t+1] <- max( (P[t]+r*P[t]*(1-P[t])-C[t]/K), 0.001) #eps=0.001
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY <- r * K / 4
                BMSY <- K / 2
                FMSY <- r/2
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior","rprior","Kprior","lqprior")
            params <- c("K", "r", "lq", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() { list(  r = mean(input$rprior),  
                                          lq = mean(input$qprior), 
                                          K = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1-input$MCMCb)) )
            .self$models[[2]] <- bugs(data=data, inits=inits,  parameters.to.save=params,  
                                   model.file = model_file, n.chains = input$MCMCc,
                                   n.iter = lchain, n.thin = input$MCMCt, 
                                   n.burnin = floor(lchain * input$MCMCb),
                                   DIC = TRUE, codaPkg = FALSE, over.relax = input$overrelax )
            if(input$do_coda){
                outS <- bugs(data=data, inits=inits, 
                             parameters.to.save=c("K", "r", "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE, over.relax = input$overrelax )
                .self$coda_models[[2]] <- read.bugs(outS)
            }
        },
        
        get_fox = function(input){
            model_file <- file.path(tempdir(), "fox_model.txt") 
            nobs    <- myglobal$nobs
            Y       <- hake$log_cpue
            C       <- hake$catch
            Ef      <- hake$effort
            lsprior <- input$sprior
            rprior  <- input$rprior
            Kprior  <- input$Kprior
            lqprior <- input$qprior
            mymodel <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs-1)) {
                    meanP[t+1] <- max((P[t]+r*P[t]*(1-log(P[t]*K)/log(K))-C[t]/K ),0.001)
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY <- r*K / (exp(1)*log(K))
                BMSY <- K / exp(1)
                FMSY <- r/log(K)
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior","rprior","Kprior","lqprior")
            params <- c("K", "r", "lq", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() { list(  r = mean(input$rprior),  
                                          lq = mean(input$qprior), 
                                          K = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1-input$MCMCb)) )
            .self$models[[3]] <- bugs(data=data, inits=inits,  parameters.to.save=params,  
                                      model.file = model_file, n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE, over.relax = input$overrelax )
            if(input$do_coda){
                outS <- bugs(data=data, inits=inits, 
                             parameters.to.save=c("K", "r", "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE, over.relax = input$overrelax )
                .self$coda_models[[3]] <- read.bugs(outS)
            }
        },
        
        get_alternative = function(input){
            model_file <- file.path(tempdir(), "alternative_model.txt") 
            nobs       <- myglobal$nobs
            Y          <- hake$log_cpue
            C          <- hake$catch
            Ef         <- hake$effort
            lsprior    <- input$sprior
            PHIprior   <- input$PHIprior
            eRHOprior  <- input$Kprior
            CHIprior   <- input$qprior
            mymodel    <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                PHI ~ dunif(PHIprior[1],PHIprior[2])
                CHI ~ dunif(CHIprior[1],CHIprior[2])
                eRHO ~ dunif(eRHOprior[1],eRHOprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                eCHI  <- exp(CHI)
                RHO   <- log(eRHO)
                #######################################
                # Process equation: A[t] = log(B[t]/K)
                #######################################
                A[1] ~ dnorm(0,tau)
                for (t in 1:(nobs-1)) {
                    meanA[t+1] <- (1-PHI) * A[t] - eCHI * Ef[t]
                    A[t+1] ~ dnorm(meanA[t+1], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- CHI + A[t] + RHO
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                BMSY <- exp( (1.0 / PHI) * log(1.0 - PHI) + RHO )
                MSY <-  BMSY * PHI / (1.0 - PHI)
                FMSY <- -log(1.0-PHI)
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior","PHIprior","eRHOprior","CHIprior")
            params <- c("eRHO", "PHI", "CHI", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() { list(  PHI = mean(input$PHIprior),  
                                          CHI = mean(input$qprior), 
                                          eRHO = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1-input$MCMCb)) )
            .self$models[[4]] <- bugs(data=data, inits=inits,  parameters.to.save=params,  
                                      model.file = model_file, n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE, over.relax = input$overrelax )
            if(input$do_coda){
                outA <- bugs(data=data, inits=inits, 
                            parameters.to.save=c("eRHO", "PHI", "CHI", "lsigma"),  
                            model.file = model_file, n.chains = input$MCMCc,
                            n.iter = lchain, n.thin = input$MCMCt, 
                            n.burnin = floor(lchain * input$MCMCb),
                            DIC = FALSE, codaPkg = TRUE, over.relax = input$overrelax )
                .self$coda_models[[4]] <- read.bugs(outA)
            }
        },
        plot_hyperparameters = function(input){
            nmod <- length(input$spm_hyper)
            all_titles <- c('Carrying capacity', 'Growth rate', 'Elasticity', 'Log-catchability',
                            'Log-variance', 'MSY', 'BMSY', 'FMSY')
            pars <- c(input$plot_KeRHO, input$plot_r, input$plot_PHI, input$plot_qCHI, input$plot_lSIGMA,
                       input$plot_msy,input$plot_bmsy,input$plot_fmsy)
            ploc <- list(pella_tomlinson=1:8, schaefer=c(1,2,NA,3:7), fox=c(1,2,NA,3:7), alternative=c(1,NA,2:7))
            plot_titles <- list(pella_tomlinson=all_titles, 
                                schaefer=c(all_titles[1:2],all_titles[4:8]), 
                                fox=c(all_titles[1:2],all_titles[4:8]), 
                                alternative=c(all_titles[1],all_titles[3:8]) )
            if(nmod>0) for(kk in 1:nmod){
                k <- which(input$spm_hyper[kk] == .self$model_names)
                for(i in 1:length(pars)) pars[i] <- (pars[i] & !is.na(ploc[[k]][i]) )
            }
            npars <- sum(pars)
            if(npars==0 | nmod==0) return()
            .self$IO$set_buffer_size(nr=npars, nc=npars)
            i <- 1
            for(ii in 1:length(pars)) if(pars[ii]) {
                j <- 1 #column index
                for(jj in 1:ii) if(pars[jj]){
                    x <- y <- mytx <- myty <- nm <- NULL
                    for(kk in 1:nmod){
                        k <- which(input$spm_hyper[kk] == .self$model_names)
                        ib <- ploc[[k]][ii]
                        jb <- ploc[[k]][jj]
                        if(length(.self$models[[k]]) > 1){
                            aux <- .self$models[[k]]$sims.matrix[,ib]
                            x <- c(x,aux)
                            if(!is.na(ib)) mytx <- plot_titles[[k]][ib]
                            nmaux <- rep(input$spm_hyper[kk], length(aux))
                            nm <- c(nm,nmaux)
                            if(jj < ii){
                                aux <-  .self$models[[k]]$sims.matrix[,jb]
                                y <- c(y,aux)
                                if(!is.na(jb)) myty <- plot_titles[[k]][jb]
                            }
                        }
                    }
                    if(jj < ii){
                        out <- data.frame(x=y, y=x, model=as.factor(nm))
                        mytitle <- c(myty,mytx)
                        .self$IO$get_scatterplot(out, mytitle, xpos=i, ypos=j)
                    } else {
                        out <- data.frame(x=x, model=as.factor(nm))
                        dolegend <- (i == 1 & j == 1 & length(input$spm_hyper)>1)
                        .self$IO$get_density_plot(out, mytitle=mytx, dolegend=dolegend, xpos=i, ypos=i)
                    }
                    j <- j + 1
                }
                i <- i+1
            }
            .self$IO$get_buffer_plot()
        },
        
        plot_data = function(input){
            do_catch  <- input$do_catch
            do_effort <- input$do_effort
            do_cpue   <- input$do_cpue
            one_row   <- input$one_row
            .self$IO$get_data_plot(do_catch, do_effort, do_cpue, one_row)
        },
        
        plot_fitted_cpue = function(input){
            fn <- function( data_vec ){return(exp(quantile(data_vec,c(0.025, 0.5, 0.975))))}
            rn <- (input$years_fit[1]:input$years_fit[2]) - 1964
            .self$IO$set_buffer_size(nr=1, nc=1)
            year <- low95 <- median <- high95 <- mname <- obs <- NULL
            obs_aux <- (hake$catch/hake$effort)[rn]
            nmod <- length(input$spm_fit)
            if(nmod >0) {
                for(kk in 1:nmod){
                    k <- which(input$spm_fit[kk] == .self$model_names)
                    if(length(.self$models[[k]]) > 1){
                        lcpue  <- .self$models[[k]]$sims.list[['mean']]
                        cpue   <- apply(X=lcpue, MARGIN=2, FUN=fn)
                        year   <- c(year, rn+1964)
                        mname  <- c(mname, rep(input$spm_fit[kk], length(rn)))
                        low95  <- c(low95,  cpue[1,rn])
                        median <- c(median, cpue[2,rn])
                        high95 <- c(high95, cpue[3,rn])
                        obs    <- c(obs, obs_aux)
                    }
                }
                out <- data.frame(year=year, low95=low95, median=median, high95=high95, 
                                  obs = obs, model = as.factor(mname))
                .self$IO$get_ts_fit_plot(out, mytitle='Observations vs. median and 95% CI', 
                                         mylabs=c('year','cpue'), xpos=1, ypos=1)
                .self$IO$get_buffer_plot()
            }
        },
        
        plot_fitted_cpue_old = function(input){
            fn <- function( data_vec ){return(exp(quantile(data_vec,c(0.025, 0.5, 0.975))))}
            nr <- if(input$spm_fit == 'Schaefer and Alternative' & !input$fit_one_row) 2 else 1
            nc <- if(input$spm_fit == 'Schaefer and Alternative' &  input$fit_one_row) 2 else 1
            k  <- 1
            rn <- (input$years_fit[1]:input$years_fit[2]) - 1964
            mt <- c('Schaefer', 'Alternative')
            .self$IO$set_buffer_size(nr=nr, nc=nc)
            for(i in 1:nr) for(j in 1:nc){
                if(length(.self$models[[k]]) > 1){
                    lcpue <- .self$models[[k]]$sims.list[['mean']]
                    cpue  <- apply(X=lcpue, MARGIN=2, FUN=fn)
                    out   <- data.frame(year=rn+1964, 
                                        obs=(hake$catch/hake$effort)[rn],
                                        low95=cpue[1,rn], median=cpue[2,rn], high95=cpue[3,rn])
                    .self$IO$get_ts_fit_plot(out, mytitle=mt[k], mylabs=c('year','cpue'), xpos=i, ypos=j)
                }
                k <- k + 1
            }
            .self$IO$get_buffer_plot()
        }
    )
)
