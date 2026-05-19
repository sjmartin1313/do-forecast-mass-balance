#' bayesmetab_mod (modified version of original bayesmetab function in package BASEmetab)
#' This script is adapted from:
#' https://github.com/dgiling/BASEmetab
#'
#' Original license: Creative Commons Attribution 3.0 (CC BY 3.0)
#'
#' Estimates single-station whole-stream metabolic rates from diel dissolved oxygen (DO) curves (see Grace et al. 2015).
#'
#' @param data.dir 	relative or absolute path to the folder containing csv input data files to be read.
#' @param results.dir 	relative or absolute path to the output folder where results (plots and tables) will be written.
#' @param interval 	Integer. The time interval in seconds (e.g. 10 minutes = 600 seconds)
#' @param n.iter 		Integer. Number of MCMC iterations (default = 20000)
#' @param n.burnin 	Integer. Number of iterations of MCMC chains to delete
#' @param update.chains 	Logical. Should the chains automatically update once if not converged? (default = TRUE)
#' @param extra.iter Numeric. Number of extra iterations to run if chains are not converged, as multiple of n.iter (default = 1 times)
#' @param smooth.DO 	Numeric. Proportion of high-frequency fluctuations to filter with fast Fourier transform (default = 0)
#' @param smooth.PAR 	Logical. Should PAR be smoothed with a moving average? (default = FALSE)
#' @param K.init 	Numeric. Initial value of chains for K (/day). Reasonable estimate aids convergence. (default value = 2)
#' @param K.est 		Logical. Should K be estimated with uninformative priors? (default = TRUE)
#' @param K.meas.mean 	Numeric. Mean for informed normal prior distribution when K.est = FALSE
#' @param K.meas.sd 	Numeric. Standard deviation for informed normal prior distribution when K.est = FALSE
#' @param p.est	Logical. Should p be estimated? (default = FALSE)
#' @param theta.est	Logical. Should theta be estimated? (default = FALSE)
#' @param instant 		Logical. Should a table of instantaneous rates be written? (default = FALSE)
#'
#' @return A dataframe and csv file of parameter estimates (mean, SD) and checks of model fit, plots of model fit (see Vignette for details https://github.com/dgiling/BASEmetab/blob/master/vignettes/BASEmetab.pdf).
#'
#' @references Grace et al. (2015) Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation). Limnology and Oceanography: Methods, 13, 103-114.
#'
#' @author Darren Giling, Ralph Mac Nally
#' @modifications by Sonny Martin for the use in "Short-Term Dissolved Oxygen Forecasting in Aquaculture Systems Using a Process-Based Mass-Balance Model"
#' Modifications made:
#' Added R.mean, R.sd, R.median in the output table.
#' Added DO.meas to the list of monitored parameters in the JAGS model.
#' Wrapped results generation in a tryCatch to avoid failure of the whole function when model outputs are invalid or incomplete.
#' Added fallback logic to create empty result rows when model outputs are invalid or incomplete.
#' Implemented helper functions for safe extraction of model outputs and Rhat values.
#' Added checks for valid model output structure before attempting to extract results.
#' Ensured that all required columns are present in the results data frame, adding NA for any missing columns.
#' Replaced direct indexing of model outputs with controlled validation to prevent dimension mismatches.
#' If a result is invalid, a warning is issued and an empty row with the correct structure is added to the output table instead of the function failing.
#' Replaced original code for calculating R2, RMSE, and other metrics with tryCatch blocks.
#' Replaced original instantaneous output logic with mean-based aggregation.
#'
#' @examples
#'
#' ##Link to JAGS
#' library(R2jags)
#' 
#' ##View example data set.
#' #set path to example data.
#' data.dir <- system.file("extdata", package = "BASEmetab")
#' ex.data <- read.csv(file.path(data.dir, "Yallakool_example.csv"))
#' head(ex.data)
#' tail(ex.data)
#'
#' ##Run Example.
#'
#' #set output directory to Output folder in current working directory.
#' results.dir <- file.path(getwd(), "Output")
#' if (dir.exists(results.dir)){} else {
#' dir.create(results.dir)}
#'
#' #run model.
#' results <- bayesmetab_mod(data.dir, results.dir, interval=600)
#'
#' @export
#' @import R2jags

bayesmetab_mod <- function(data.dir, results.dir, interval, n.iter=20000, n.burnin=n.iter*0.5, K.init = 2, 
                      smooth.DO=0, smooth.PAR=FALSE, instant=FALSE, update.chains = TRUE, extra.iter=1,
                      K.est = TRUE, K.meas.mean = 0, K.meas.sd = 4, p.est=FALSE, theta.est=FALSE) 
  {
  start.time<-NULL; start.time<-Sys.time()
  
  # model file
  model.dir <- system.file("tools", package = "BASEmetab")
  
  # define functions
  smooth5 <- function(x) (zoo::rollapply(x, 5, mean, na.rm=T,align="center"))  # moving average of 5 time intervals
  
  # data input and set up output table dataframes
  filenames<-list.files(file.path(data.dir))  
  
  # Set up output tables
  output.table<-NULL
  output.table<-data.frame(File=character(), Date=character(), 
                           GPP.mean=double(), GPP.sd=double(), GPP.median=double(),
                           ER.mean=double(), ER.sd=double(), ER.median=double(), 
                           NEP.mean=double(), NEP.sd=double(), NEP.median=double(), 
                           PR.mean=double(), PR.sd=double(), PR.median=double(),
                           K.mean=double(), K.sd=double(), K.median=double(),
                           theta.mean=double(), theta.sd=double(), theta.median=double(),
                           A.mean=double(), A.sd=double(), A.median=double(), 
                           p.mean=double(), p.sd=double(), p.median=double(),
                           R.mean=double(), R.sd=double(), R.median=double(),
                           R2=double(), PPP=double(), rmse=double(), rmse.relative=double(), mrl.fraction=double(), ER.K.cor=double(), 
                           convergence.check=double(), A.Rhat=double(), K.Rhat=double(), theta.Rhat=double(), p.Rhat=double(), R.Rhat=double(), GPP.Rhat=double(), 
                           DIC=double(), pD=double(),
                           totDailyLight=double(), aveDailyTemp=double(), 
                           interval= double(), smooth.DO=double() , smooth.PAR=logical(), n.iter= double(), n.burnin= double(),
                           stringsAsFactors=FALSE)
  instant.rates<-data.frame(File=character(), Date=character(), interval=integer(), 
                            tempC=double(), I=double(), K.instant=double(), GPP.instant=double(), ER.instant=double(),
                            stringsAsFactors=FALSE)
  
  # Analyse files sequentially
  for (fname in filenames) {
    #fname <- filenames[1]
    
    data<-read.csv(file.path(data.dir,fname), head=T) # read next file
    seconds<-86400
    req.rows <- 86400/interval
    N = nrow(data)
    x = 0:(N-1)

    ## checks
      # check headers
      if (! all(colnames(data) %in% c("Date","Time","I", "tempC", "DO.meas", "atmo.pressure", "salinity"))) {
          stop(paste0("Column headers in input file '", fname, "' do not include 'Date', 'Time', 'I', 'tempC', 'DO.meas', 'atmo.pressure' and 'salinity' (they are case sensitive)")) }
  
      # check dates for "/" and replace with "-"
      data$Date <- gsub("/", "-", data$Date)
        
    ## Smoothing data
    if(smooth.DO > 0) {
      # fast Fourier transform smoothing - low pass filter
      DO.fft = fft(data$DO.meas)
      inx_filter = floor(N/2*(1-smooth.DO))
      filter = rep(1, N)
      filter[inx_filter:(N-inx_filter)] = 0
      DO.fft_filtered = filter * DO.fft
      data$DO.smooth <- Re( fft( DO.fft_filtered, inverse=TRUE) / N )
      noise <- data$DO.meas - data$DO.smooth
    }
    
    if(smooth.PAR == T) {
      data$I.smooth<-c(data$I[1:2],smooth5(data$I),data$I[nrow(data)-1],data$I[nrow(data)]) # moving average over 5 time intervals
    }
    
    # Select dates
    data$Date <- factor(data$Date, levels = unique(data$Date))
    dates <- unique(data$Date)
    n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
    incomp.dates <- dates[n.records != (seconds/interval)]
    if(length(incomp.dates)>0) {
      warning(paste0("Not all dates in file '", fname, "' contain ", req.rows, " rows"))
    }
    dates <- dates[n.records == (seconds/interval)] # select only dates with full days
    
    ## Analyse days sequentially
    for (d in dates) 
    { 
      #d <- dates[2]
      data.sub <- data[data$Date == d,]
      
      # Define data vectors
      num.measurements <- nrow(data.sub)
      tempC <- data.sub$tempC
      salinity <-data.sub$salinity
      atmo.pressure <- data.sub$atmo.pressure
      DO.meas <- if(smooth.DO > 0) data.sub$DO.smooth else data.sub$DO.meas
      PAR     <- if(smooth.PAR == TRUE) data.sub$I.smooth else data.sub$I
      
      # Initial values
      # Set these to something sensible if the model is becoming stuck in a bad parameter space
      # These values here are expressed per timestep, not per day. Divide desired initial K (/day) by the number of timesteps in a day, as shown in default below 
      inits <- function()       {       list(K = K.init / (86400/interval) ) }
      
      # Different random seeds
      kern=as.integer(runif(1000,min=1,max=10000))
      iters=sample(kern,1)
      
      # Set 
      n.chains <- 3
      n.thin <- 10
      p.est.n <- as.numeric(p.est)
      theta.est.n <- as.numeric(theta.est)
      K.est.n <- as.numeric(K.est)
      K.meas.mean.ts <- K.meas.mean / (86400/interval)
      K.meas.sd.ts <- K.meas.sd / (86400/interval)
      data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","salinity","atmo.pressure", "K.init", 
                        "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")  

      # Define monitoring variables
      params=c("A","R","K","K.day","p","theta","tau","ER","GPP","NEP","PR","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled",
               "gppts", "erpts", "kpts", "DO.meas")

      ## Call jags ##
      
      # Set debug = T below to inspect each file for model convergence 
      # (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
      metabfit=NULL
      metabfit <- do.call(R2jags::jags.parallel,
                          list(data=data.list, inits=inits, parameters.to.save=params, model.file = file.path(system.file(package="BASEmetab"), "BASE_metab_model_v2.3.txt"),
                               n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                               n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
                               jags.seed = 123, digits=5))
      
      # print(metabfit, digits=2) # to inspect results of last metabfit
      
      ## diagnostic summaries
      # Rhat (srf) test
      srf<- metabfit$BUGSoutput$summary[,8]
      Rhat.test <- NULL
      Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
      
      # Check for convergence and update once if requested
      if(update.chains == TRUE) {
        if(Rhat.test == "Check convergence") {
          recompile(metabfit)
          metabfit <- update(metabfit, n.iter=n.iter*extra.iter) 
          
          # Rhat (srf) test - second round in case metabfit is updated
          srf<- metabfit$BUGSoutput$summary[,8]
          Rhat.test <- NULL
          Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
        }
      }
      
      # autocorr test
      metabfit.mcmc<-coda::as.mcmc(metabfit)
      ac.lag1 <- coda::autocorr.diag(metabfit.mcmc, lags = 1)
      auto.corr.test <- NULL
      auto.corr.test <- ifelse(any(abs(ac.lag1)>0.1, na.rm=T)==TRUE,"Check ac", "ac OK")
      
      PPP <- metabfit$BUGSoutput$summary["PPfit","mean"] # posterior predictive p-value
      
      DO.mod.means <- metabfit$BUGSoutput$mean$DO.modelled
      DO.mod.sd <- metabfit$BUGSoutput$sd$DO.modelled
      
      R2 = cor(DO.mod.means,DO.meas)^2
      rmse = sqrt(sum((metabfit$BUGSoutput$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
      post.mean.dev <- metabfit$BUGSoutput$mean$deviance
      pD <- metabfit$BUGSoutput$pD
      DIC <- metabfit$BUGSoutput$DIC
      
      DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
      ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
      rmse.relative <- rmse / ptpvar
      
      diff<-metabfit$BUGSoutput$mean$DO.modelled-DO.meas
      mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
      mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # proportion of largest run
      
      ER.K.cor <- cor(metabfit$BUGSoutput$sims.list$ER,metabfit$BUGSoutput$sims.list$K) # plot(metabfit$sims.list$ER ~ metabfit$sims.list$K)
      
      # Safely create results row
      result <- tryCatch({
        # First verify we have valid BUGS output
        if(is.null(metabfit) || is.null(metabfit$BUGSoutput) || is.null(metabfit$BUGSoutput$mean)) {
            stop("Invalid model output structure")
        }

        # Create safe extraction functions
        safe_extract <- function(obj, param, type = "mean") {
            tryCatch({
            if(param %in% names(obj$BUGSoutput[[type]])) {
                val <- obj$BUGSoutput[[type]][[param]]
                if(length(val) > 0) val else NA_real_
            } else {
                NA_real_
            }
            }, error = function(e) NA_real_)
        }
        
        safe_rhat <- function(obj, param) {
            tryCatch({
            if(param %in% rownames(obj$BUGSoutput$summary)) {
                obj$BUGSoutput$summary[param, "Rhat"]
            } else {
                NA_real_
            }
            }, error = function(e) NA_real_)
        }
        
        # Calculate model fit metrics
        DO.mod.means <- safe_extract(metabfit, "DO.modelled")
        R2_val <- tryCatch(cor(DO.mod.means, DO.meas)^2, error = function(e) NA_real_)
        rmse_val <- tryCatch(sqrt(mean((DO.mod.means - DO.meas)^2)), error = function(e) NA_real_)
        
        # Calculate additional metrics
        ptpvar <- tryCatch({
            DO.lag <- DO.meas[2:length(DO.meas)] - DO.meas[1:(length(DO.meas)-1)]
            sqrt(sum((DO.lag)^2)/(length(DO.meas)-1))
        }, error = function(e) NA_real_)
        
        rmse.relative <- if(!is.na(rmse_val) && !is.na(ptpvar)) rmse_val/ptpvar else NA_real_
        
        mrl.fraction <- tryCatch({
            diff <- DO.mod.means - DO.meas
            max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas)
        }, error = function(e) NA_real_)
        
        ER.K.cor <- tryCatch({
            cor(metabfit$BUGSoutput$sims.list$ER, metabfit$BUGSoutput$sims.list$K)
        }, error = function(e) NA_real_)
        
        # Create the results data frame
        res <- data.frame(
            File = as.character(fname),
            Date = as.character(d),
            GPP.mean = safe_extract(metabfit, "GPP"),
            GPP.sd = safe_extract(metabfit, "GPP", "sd"),
            GPP.median = safe_extract(metabfit, "GPP", "median"),
            ER.mean = safe_extract(metabfit, "ER"),
            ER.sd = safe_extract(metabfit, "ER", "sd"),
            ER.median = safe_extract(metabfit, "ER", "median"),
            NEP.mean = safe_extract(metabfit, "NEP"),
            NEP.sd = safe_extract(metabfit, "NEP", "sd"),
            NEP.median = safe_extract(metabfit, "NEP", "median"),
            PR.mean = safe_extract(metabfit, "PR"),
            PR.sd = safe_extract(metabfit, "PR", "sd"),
            PR.median = safe_extract(metabfit, "PR", "median"),
            K.mean = safe_extract(metabfit, "K.day"),
            K.sd = safe_extract(metabfit, "K.day", "sd"),
            K.median = safe_extract(metabfit, "K.day", "median"),
            theta.mean = safe_extract(metabfit, "theta"),
            theta.sd = safe_extract(metabfit, "theta", "sd"),
            theta.median = safe_extract(metabfit, "theta", "median"),
            A.mean = safe_extract(metabfit, "A"),
            A.sd = safe_extract(metabfit, "A", "sd"),
            A.median = safe_extract(metabfit, "A", "median"),
            p.mean = safe_extract(metabfit, "p"),
            p.sd = safe_extract(metabfit, "p", "sd"),
            p.median = safe_extract(metabfit, "p", "median"),
            R.mean = safe_extract(metabfit, "R"),
            R.sd = safe_extract(metabfit, "R", "sd"),
            R.median = safe_extract(metabfit, "R", "median"),
            R2 = R2_val,
            PPP = safe_extract(metabfit, "PPfit"),
            rmse = rmse_val,
            rmse.relative = rmse.relative,
            mrl.fraction = mrl.fraction,
            ER.K.cor = ER.K.cor,
            convergence.check = as.numeric(Rhat.test == "Check convergence"),
            A.Rhat = safe_rhat(metabfit, "A"),
            K.Rhat = safe_rhat(metabfit, "K.day"),
            theta.Rhat = safe_rhat(metabfit, "theta"),
            p.Rhat = safe_rhat(metabfit, "p"),
            R.Rhat = safe_rhat(metabfit, "R"),
            GPP.Rhat = safe_rhat(metabfit, "GPP"),
            DIC = if(!is.null(metabfit$BUGSoutput$DIC)) metabfit$BUGSoutput$DIC else NA_real_,
            pD = if(!is.null(metabfit$BUGSoutput$pD)) metabfit$BUGSoutput$pD else NA_real_,
            totDailyLight = sum(PAR, na.rm = TRUE),
            aveDailyTemp = mean(tempC, na.rm = TRUE),
            interval = interval,
            smooth.DO = smooth.DO,
            smooth.PAR = smooth.PAR,
            n.iter = n.iter,
            n.burnin = n.burnin,
            stringsAsFactors = FALSE
        )
        
        # Verify all required columns are present
        required_cols <- names(output.table)
        missing_cols <- setdiff(required_cols, names(res))
        if(length(missing_cols) > 0) {
            res[missing_cols] <- NA_real_
        }
        
        res
      }, error = function(e) {
        warning("Error creating results for date ", d, ": ", e$message)
        # Create minimal result row
        empty_row <- output.table[1, , drop = FALSE]
        empty_row[1, ] <- NA
        empty_row$File <- fname
        empty_row$Date <- d
        empty_row
      })

      # Add to output table
      if(!is.null(result) && nrow(result) == 1 && ncol(result) == ncol(output.table)) {
        output.table <- rbind(output.table, result)
      } else {
        warning("Invalid result for date ", d)
        # Create empty row with correct structure
        empty_row <- output.table[1, , drop = FALSE]
        empty_row[1, ] <- NA
        empty_row$File <- fname
        empty_row$Date <- d
        output.table <- rbind(output.table, empty_row)
      }

      if(instant == TRUE) {
        param_means <- data.frame(
          File = as.character(rep(fname, seconds/interval)),
          Date = as.character(rep(d, seconds/interval)),
          time = 1:num.measurements,
          K = colMeans(metabfit$BUGSoutput$sims.list$K_instant),
          GPP = colMeans(metabfit$BUGSoutput$sims.list$gppts),
          ER = colMeans(metabfit$BUGSoutput$sims.list$erpts),
          KP = colMeans(metabfit$BUGSoutput$sims.list$kpts),
          DO.modelled = colMeans(metabfit$BUGSoutput$sims.list$DO.modelled),
          DO.meas = colMeans(metabfit$BUGSoutput$sims.list$DO.meas)
        )

        write.csv(param_means, file = file.path(results.dir, paste0("inst_rates_", d, ".csv")), row.names = FALSE)
      }
      
      # diagnostic traceplots and scatterplots
      jpeg(file=file.path(results.dir, paste0(substr(fname, 1,(nchar(fname)-4)),"_", as.character(d), "_1_trace_", gsub(":","",Sys.time()), ".jpg")), width=1200, height=1200, pointsize=30, quality=300)
      
      par(mfrow=c(3,3), mar=c(3,4,2,1), oma=c(0.1,0.1,0.1,0.1))
      traceplot(metabfit, varname=c('A','p','R','K.day','theta'), ask=FALSE, mfrow=c(3,3), new=FALSE)
      
      plot(1:num.measurements,data.sub$DO.meas, type="p",pch=21, col="grey60",cex=0.8, ylim=c(min(DO.mod.means-DO.mod.sd)-0.5,max(DO.mod.means+DO.mod.sd)+0.5), xlab="Timestep", ylab="DO mg/L")
      if(smooth.DO > 0) { points(1:num.measurements,data.sub$DO.smooth,type='l',lwd=2,xlab="Timestep", col="red", cex=0.75)  }
      points(1:num.measurements,DO.mod.means,lwd=1.5, type="l", xlab="Timestep", col="black")  
      points(1:num.measurements,DO.mod.means+DO.mod.sd, type="l", lty=2)
      points(1:num.measurements,DO.mod.means-DO.mod.sd, type="l", lty=2)
      legend(x="topleft", legend=c("DO meas", "DO smooth", "DO fit"), pch=c(1,NA,NA), lty=c(NA,1,1), col=c("grey60","red", "black"), cex=0.75, bty='n')
      
      plot(1:num.measurements,tempC,pch=1,xlab="Timestep" , typ='p', col="grey60")
      legend(x="topleft", legend=c("TempC meas"), pch=c(1), col=c("grey60"), cex=0.75, bty='n')
      
      plot(1:num.measurements,data.sub$I,pch=1,xlab="Timestep" , typ='p', col="grey60", ylab='PAR')
      if(smooth.PAR==TRUE) { points(1:num.measurements,data.sub$I.smooth,type='l',lwd=2,xlab="Timestep", col="red", cex=0.75)  }
      legend(x="topleft", legend=c("PAR meas", "PAR smooth"), pch=c(1,NA), lty=c(NA,1), col=c("grey60","red"), cex=0.75, bty='n')
      
      graphics.off()
      
      # diagnostic density plots
      jpeg(file=file.path(results.dir, paste0(substr(fname, 1,(nchar(fname)-4)),"_", as.character(d), "_2_density_", gsub(":","",Sys.time()), ".jpg")), width=1200, height=1200, pointsize=30, quality=300)

      print(lattice::densityplot(metabfit.mcmc[,c("theta", "K.day", "p", "A", "ER", "GPP")]))

      graphics.off()
      
    }
  }
  
  write.csv(output.table, file=file.path(results.dir, gsub(":","",paste0("BASE_results_",Sys.time(),".csv")))) # write output file
  if(instant == TRUE) {
    write.csv(instant.rates, file=file.path(results.dir, gsub(":","",paste0("instantaneous_rates_",Sys.time(),".csv")))) # write output file
  }
  
  end.time<-NULL; end.time<-Sys.time()
  print(end.time-start.time)
  return(output.table)
  
}