###POWER ANALYSIS FOR G4 CALCULATIONS

###INPUT PARAMETERS
###
###
###THE FOUR MAIN PARAMETERS ARE SENSITIVITY, SPECIFICITY, PREVALENCE, AND BIAS
###IF ONE KNOWS THREE OF THE ABOVE FOUR PARAMETERS, THE FOURTH ONE IS CALCULATED AUTOMATICALLY
###ONLY THREE OF THE FOUR PARAMETERS CAN BE LISTED AT A TIME; ONE PARAMETER MUST REMAIN AS NULL
###BY DEFAULT, SENSITIVITY / SPECIFICITY ARE SET AT 80%, AND PREVALENCE IS SET AT 50%
###N = SAMPLE SIZE FOR THE STUDY (THE TOTAL SAMPLE SIZE FOR THE CONFUSION MATRIX); DEFAULT IS 100
###mu0 = NULL HYPOTHESIS FOR G4 SIGNIFICANCE TESTING; THE POWER ANALYSIS PERFORMS A TWO-TAILED WALD TEST; DEFAULT IS 0.60
###THE POINT ESTIMATE FOR G4 IS BIAS-CORRECTED TO ACCOUNT FOR THE BOOTSTRAP SAMPLING ERROR
###Alpha = TYPE 1 ERROR RATE / STATISTICAL SIGNIFICANCE THRESHOLD; DEFAULT IS 5%
###Repeats = THE NUMBER OF ITERATIONS FOR THE SIMULATION; THE GREATER THE NUMBER, THE GREATER THE PRECISION OF THE POWER ESTIMATE; DEFAULT IS 2500
###AS THE NUMBER OF REPEATS INCREASES, THE TIME TO COMPLETION INCREASES. AT LEAST 1000 REPEATS ARE RECOMMENDED TO MAINTAIN STABLE ESTIMATES
###PARALLEL PROCESSING IS AVAILABLE TO HELP SPEED UP COMPUTATION TIME; INSTRUCTIONS ON THE USAGE OF AVAILABLE CLUSTERS ARE FOUND BELOW
###Boots = FOR EACH REPEAT, THE NUMBER OF REPLICATIONS TO GENERATE A BOOTSTRAP DISTRIBUTION FOR G4 (THIS WILL HELP GENERATE THE 95% CI); DEFAULT IS 250
###Keep.Prevalence = SHOULD THE PREVALENCE RATE BE MAINTAINED DURING THE BOOTSTRAP PROCESS? DEFAULT IS TRUE
###WHEN COMPARING A TEST GROUP TO A GROUND TRUTH / REFERENCE PANEL, IT MAKES SENSE TO MAINTAIN THE PREVALENCE RATE SO THAT EACH BOOTSTRAP RESAMPLE REFLECTS
###THE ORIGINAL SAMPLE'S DISTRIBUTION OF TRUE CLASSIFICATIONS
###IF THE COMPARATOR GROUP IS NOT A GROUND TRUTH / REFERENCE PANEL (FOR EXAMPLE, ANOTHER RATER WHO MAY OR MAY NOT BE CORRECT), THEN "Keep.Prevalence"
###SHOULD BE SET TO FALSE
###
###
###


library(data.table)
library(doSNOW)
library(doRNG)
library(parallel)
library(foreach)
library(splitstackshape)
library(tidyverse)


G4.Power <- function(Sensitivity = 0.80, Specificity = 0.80, Prevalence = 0.50, Bias = NULL, N = 100, mu0 = 0.60, Alpha = 0.050,
			   Repeats = 2500, Boots = 250, Keep.Prevalence = TRUE) {

      if (!is.null(Sensitivity) & !is.null(Specificity) & !is.null(Prevalence) & !is.null(Bias)) {
      print("Error - Only three of the four metrics may be inputted. Please remove one of them to continue")
      }

	pb <- txtProgressBar(max = Repeats, style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)
  
      power.g4 <- foreach(k = 1:Repeats, .packages = c("data.table", "doSNOW", "doRNG", "parallel", "splitstackshape", "tidyverse"),
			.options.snow = opts) %dopar% {
      	
      if (!is.null(Sensitivity) & is.null(Specificity) & !is.null(Prevalence) & !is.null(Bias)) {
      New.Sensitivity = rbeta(1, Sensitivity * 250, (1 - Sensitivity) * 250)
      tp = floor(Prevalence * N * New.Sensitivity)
      fn = floor(Prevalence * N - tp)
      fp = floor(Bias * N - tp)
      tn = N - tp - fn - fp
    } else if (is.null(Sensitivity) & !is.null(Specificity) & !is.null(Prevalence) & !is.null(Bias)) {
      New.Specificity = rbeta(1, Specificity * 250, (1 - Specificity) * 250)
      tn = floor((1 - Prevalence) * N * New.Specificity)
      fp = floor((1 - Prevalence) * N - tn)
      tp = floor(Bias * N - fp)
      fn = N - tn - fp - tp
    } else if (!is.null(Sensitivity) & !is.null(Specificity) & !is.null(Prevalence) & is.null(Bias)) {
      New.Sensitivity = rbeta(1, Sensitivity * 250, (1 - Sensitivity) * 250)
      New.Specificity = rbeta(1, Specificity * 250, (1 - Specificity) * 250)
      tp = floor(Prevalence * N * New.Sensitivity)
      tn = floor((1 - Prevalence) * N * New.Specificity)
      fn = floor(Prevalence * N - tp)
      fp = N - tp - tn - fn
    } else {
      stop("Error - Insufficient input parameters")
    }

      g4 <- sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)

      data <- data.table(Value = c(rep("tp", tp), rep("fn", fn), rep("fp", fp), rep("tn", tn)), Truth = c(rep("positive", tp + fn), rep("negative", tn + fp)))
	
	balance.boot.data = data[rep(seq_len(nrow(data)), Boots), ]
	iboot = sample(1:nrow(balance.boot.data), replace = FALSE)
	balance.boot.data = balance.boot.data[iboot, ]
	  
	matrix.pos = matrix(balance.boot.data[Truth == "positive", Value], nrow = data[Truth == "positive", .N], ncol = Boots)
      matrix.neg = matrix(balance.boot.data[Truth == "negative", Value], nrow = data[Truth == "negative", .N], ncol = Boots)
      
      matrix.all = matrix(balance.boot.data[, Value], nrow = data[, .N], ncol = Boots)
		
      boot.g4 <- foreach(i = 1:Boots) %do% {
      	
      if (Keep.Prevalence == TRUE) {
      	
      bootdata1 <- matrix.pos[, i]
      bootdata2 <- matrix.neg[, i]
      bootdata <- c(bootdata1, bootdata2)

      }
      	
      else if (Keep.Prevalence == FALSE) {
      	
	bootdata <- matrix.all[, i]
		
	}

	tp.boot = length(bootdata[bootdata == "tp"])
	fn.boot = length(bootdata[bootdata == "fn"])
	tn.boot = length(bootdata[bootdata == "tn"])
	fp.boot = length(bootdata[bootdata == "fp"])

      sqrt((tp.boot * tn.boot) / ((tp.boot + fn.boot) * (tn.boot + fp.boot))) * ((tp.boot + fn.boot) * (tn.boot + fp.boot) / ((tp.boot + fp.boot) * (tn.boot +
      fn.boot))) ^ (1/4)
      }

      results <- data.table(
        G4.Metric = round(g4, 5),
        CI.Lower = round(2 * g4 - mean(unlist(boot.g4)) - qnorm(1 - Alpha/2) * sd(unlist(boot.g4)), 5),
        CI.Upper = round(2 * g4 - mean(unlist(boot.g4)) + qnorm(1 - Alpha/2) * sd(unlist(boot.g4)), 5)
      )

      results[, CI.Lower] > mu0
    }

    close(pb)
    print(mean(unlist(power.g4)))
  
}
 			
###SIMULATION IS BASED ON 250 BOOTSTRAP ITERATIONS FOR THE 95% CI, WITH 2500 REPEATS TO ASSESS THE POWER
###WITH THE ABOVE PARAMETERS, SIMULATION TIME IS APPROXIMATELY 30 SECONDS WITH 7 CORES

###SET THE TOTAL NUMBER OF CLUSTERS (I.E., "WORKERS") FOR PARALLEL PROCESSING

detectCores()

cluster.total = detectCores() - 1

###PERFORM PARALLEL PROCESSING TO SPEED UP COMPUTATIONS

cl <- makePSOCKcluster(cluster.total)
registerDoSNOW(cl)

start.time <- Sys.time()
start.time

set.seed(12345, kind = "L'Ecuyer-CMRG")

G4.Power(Sensitivity = 0.90, Specificity = 0.81, Prevalence = 0.07, Bias = NULL, N = 663, mu0 = 0.60, Alpha = 0.050, Repeats = 2500,
Boots = 250, Keep.Prevalence = TRUE)
#G4.Power(Sensitivity = 0.92, Specificity = 0.88, Prevalence = 0.07, Bias = NULL, N = 663, mu0 = 0.60, Alpha = 0.050, Repeats = 2500,
###Boots = 250, Keep.Prevalence = TRUE)
#G4.Power(Sensitivity = 37/46, Specificity = 30/34, Prevalence = 0.60, Bias = NULL, N = 220, mu0 = 0.70, Alpha = 0.050, Repeats = 2500,
###Boots = 250, Keep.Prevalence = TRUE)
#G4.Power(Sensitivity = 37/46, Specificity = 30/34, Prevalence = 0.60, Bias = NULL, N = 963, mu0 = 0.80, Alpha = 0.050, Repeats = 2500,
###Boots = 250, Keep.Prevalence = TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
stopCluster(cl)









