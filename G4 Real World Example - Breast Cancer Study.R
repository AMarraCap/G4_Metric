###EXAMPLE SHOWCASING THE BENEFICIAL EFFECTS OF G4 / P4 / MCC OVER AUROC
###FOR MORE INFORMATION, SEE "https://www.nature.com/articles/s41467-021-26023-2"

library(data.table)
library(tidyverse)
library(splitstackshape)
library(caret)
library(irr)
library(doParallel)
library(doRNG)
library(parallel)
library(foreach)


###EXAMS BY RADIOLOGISTS WITHOUT AI ASSISTANCE

tp = 44
fn = 5
tn = 497
fp = 117
N = tp + fp + tn + fn
prev = (tp + fn) / N
bias = (tp + fp) / N
sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv
kappa = (2 * (tp * tn - fp * fn)) / ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn))
p4 = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
g4 = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
matthews = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)

unaided = data.table(Value = factor(c(rep("tp", tp), rep("fn", fn), rep("tn", tn), rep("fp", fp))),
		  Truth = c(rep("positive", tp + fn), rep("negative", tn + fp)))

###SET THE TOTAL NUMBER OF CLUSTERS (I.E., "WORKERS") FOR PARALLEL PROCESSING

detectCores()

cluster.total = detectCores() - 1

###PERFORM PARALLEL PROCESSING TO SPEED UP COMPUTATIONS

cl <- makePSOCKcluster(cluster.total)
registerDoParallel(cl)

start.time <- Sys.time()
start.time

set.seed(71123, kind = "L'Ecuyer-CMRG")

multiResultClass <- function(p4 = NULL, g4 = NULL, matthews = NULL)
{
  me <- list(
    p4 = p4,
    g4 = g4,
    matthews = matthews
  )

  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}

output <- foreach(k = 1:100000, .packages = c("data.table", "doParallel", "parallel", "tidyverse", "splitstackshape")) %dopar% {
	
result <- multiResultClass()
	
bootdata <- unaided %>% stratified(., group = c("Truth"), size = 0.99999, replace = TRUE)

#iboot = sample(1:nrow(unaided), replace = TRUE)
#bootdata = unaided[iboot, ]
tp = nrow(bootdata[Value == "tp", ])
fn = nrow(bootdata[Value == "fn", ])
tn = nrow(bootdata[Value == "tn", ])
fp = nrow(bootdata[Value == "fp", ])

sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv

result$p4 = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
result$g4 = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
result$matthews = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)
return(result)

}

p4 = NULL
g4 = NULL
matthews = NULL

for (i in 1:100000) {
	
pp4 = output[[i]]$p4
p4 = c(p4, pp4)	

gg4 = output[[i]]$g4
g4 = c(g4, gg4)

mcc = output[[i]]$matthews
matthews = c(matthews, mcc)

}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
stopCluster(cl)
	


hist(p4)
mean(p4)
quantile(p4, c(0.025, 0.975))

hist(g4)
mean(g4)
quantile(g4, c(0.025, 0.975))

hist(matthews)
mean(matthews)
quantile(matthews, c(0.025, 0.975))
hist((matthews + 1) / 2)
mean((matthews + 1) / 2)
quantile((matthews + 1) / 2, c(0.025, 0.975))



###EXAMS BY RADIOLOGISTS WITH AI ASSISTANCE

tp = 45
fn = 4
tn = 542
fp = 72
N = tp + fp + tn + fn
prev = (tp + fn) / N
bias = (tp + fp) / N
sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv
kappa = (2 * (tp * tn - fp * fn)) / ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn))
p4 = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
g4 = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
matthews = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)

aided = data.table(Value = factor(c(rep("tp", tp), rep("fn", fn), rep("tn", tn), rep("fp", fp))),
		Truth = c(rep("positive", tp + fn), rep("negative", tn + fp)))
		
###SET THE TOTAL NUMBER OF CLUSTERS (I.E., "WORKERS") FOR PARALLEL PROCESSING

detectCores()

cluster.total = detectCores() - 1

###PERFORM PARALLEL PROCESSING TO SPEED UP COMPUTATIONS

cl <- makePSOCKcluster(cluster.total)
registerDoParallel(cl)

start.time <- Sys.time()
start.time

set.seed(71123, kind = "L'Ecuyer-CMRG")

multiResultClass <- function(p4 = NULL, g4 = NULL, matthews = NULL)
{
  me <- list(
    p4 = p4,
    g4 = g4,
    matthews = matthews
  )

  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}

output <- foreach(k = 1:100000, .packages = c("data.table", "doParallel", "parallel", "tidyverse", "splitstackshape")) %dopar% {
	
result <- multiResultClass()
	
bootdata <- aided %>% stratified(., group = c("Truth"), size = 0.99999, replace = TRUE)

#iboot = sample(1:nrow(aided), replace = TRUE)
#bootdata = aided[iboot, ]
tp = nrow(bootdata[Value == "tp", ])
fn = nrow(bootdata[Value == "fn", ])
tn = nrow(bootdata[Value == "tn", ])
fp = nrow(bootdata[Value == "fp", ])

sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv

result$p4 = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
result$g4 = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
result$matthews = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)
return(result)

}

p4 = NULL
g4 = NULL
matthews = NULL

for (i in 1:100000) {
	
pp4 = output[[i]]$p4
p4 = c(p4, pp4)	

gg4 = output[[i]]$g4
g4 = c(g4, gg4)

mcc = output[[i]]$matthews
matthews = c(matthews, mcc)

}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
stopCluster(cl)
	


hist(p4)
mean(p4)
quantile(p4, c(0.025, 0.975))

hist(g4)
mean(g4)
quantile(g4, c(0.025, 0.975))

hist(matthews)
mean(matthews)
quantile(matthews, c(0.025, 0.975))
hist((matthews + 1) / 2)
mean((matthews + 1) / 2)
quantile((matthews + 1) / 2, c(0.025, 0.975))	
		


###DIFFERENCE BETWEEN AIDED AND UNAIDED GROUPS

data = data.table(Value = c(aided[, Value], unaided[, Value]), Modality = c(rep("Aided", nrow(aided)), rep("Unaided", nrow(unaided))))
data[, Truth := "positive"]
data[Value == "tp", Truth := "positive"][Value == "fn", Truth := "positive"][Value == "tn", Truth := "negative"][Value == "fp", Truth := "negative"]
data$Modality <- factor(data$Modality)
data$Truth <- factor(data$Truth)
str(data)

###SET THE TOTAL NUMBER OF CLUSTERS (I.E., "WORKERS") FOR PARALLEL PROCESSING

detectCores()

cluster.total = detectCores() - 1

###PERFORM PARALLEL PROCESSING TO SPEED UP COMPUTATIONS

cl <- makePSOCKcluster(cluster.total)
registerDoParallel(cl)

start.time <- Sys.time()
start.time

set.seed(71123, kind = "L'Ecuyer-CMRG")

multiResultClass <- function(diff.p4 = NULL, diff.g4 = NULL, diff.mcc = NULL, diff.mcc.transformed = NULL)
{
  me <- list(
    diff.p4 = diff.p4,
    diff.g4 = diff.g4,
    diff.mcc = diff.mcc,
    diff.mcc.transformed = diff.mcc.transformed
  )

  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}

output <- foreach(k = 1:100000, .packages = c("data.table", "doParallel", "parallel", "tidyverse", "splitstackshape")) %dopar% {
	
result <- multiResultClass()
	
iboot <- data %>% stratified(., group = c("Modality", "Truth"), size = 0.99999, replace = TRUE)

boot.aid = iboot[Modality == "Aided", ]
tp = nrow(boot.aid[Value == "tp", ])
fn = nrow(boot.aid[Value == "fn", ])
tn = nrow(boot.aid[Value == "tn", ])
fp = nrow(boot.aid[Value == "fp", ])

sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv

p4.aid = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
g4.aid = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
mcc.aid = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)
mcc.aid.transformed = (mcc.aid + 1) / 2

boot.unaid = iboot[Modality == "Unaided", ]
tp = nrow(boot.unaid[Value == "tp", ])
fn = nrow(boot.unaid[Value == "fn", ])
tn = nrow(boot.unaid[Value == "tn", ])
fp = nrow(boot.unaid[Value == "fp", ])

sens = tp / (tp + fn)
specif = tn / (tn + fp)
ppv = tp / (tp + fp)
npv = tn / (tn + fn)
fnr = 1 - sens
fpr = 1 - specif
fdr = 1 - ppv
fomitr = 1 - npv

p4.unaid = 4 / (1/sens + 1/specif + 1/ppv + 1/npv)
g4.unaid = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
mcc.unaid = sqrt(sens * specif * ppv * npv) - sqrt(fnr * fpr * fdr * fomitr)
mcc.unaid.transformed = (mcc.unaid + 1) / 2

result$diff.p4 = p4.aid - p4.unaid
result$diff.g4 = g4.aid - g4.unaid
result$diff.mcc = mcc.aid - mcc.unaid
result$diff.mcc.transformed = mcc.aid.transformed - mcc.unaid.transformed
return(result)

}

diff.p4 = NULL
diff.g4 = NULL
diff.mcc = NULL
diff.mcc.transformed = NULL

for (i in 1:100000) {
	
a = output[[i]]$diff.p4
diff.p4 = c(diff.p4, a)	

b = output[[i]]$diff.g4
diff.g4 = c(diff.g4, b)

c = output[[i]]$diff.mcc
diff.mcc = c(diff.mcc, c)

d = output[[i]]$diff.mcc.transformed
diff.mcc.transformed = c(diff.mcc.transformed, d)

}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
stopCluster(cl)



hist(diff.p4)
mean(diff.p4)
quantile(diff.p4, c(0.025, 0.975))

hist(diff.g4)
mean(diff.g4)
quantile(diff.g4, c(0.025, 0.975))

hist(diff.mcc)
mean(diff.mcc)
quantile(diff.mcc, c(0.025, 0.975))

hist(diff.mcc.transformed)
mean(diff.mcc.transformed)
quantile(diff.mcc.transformed, c(0.025, 0.975))



###WITH A PREVALENCE RATE OF 7%, 663 EXAMS WERE CONDUCTED TO DETERMINE IF BREAST CANCER WAS PRESENT IN PATIENTS
###WITHOUT AI ASSISTANCE, THE AUROC WAS 0.924 +- 0.02
###WITH AI ASSISTANCE, THE AUROC WAS 0.962

###FOR MORE INFORMATION, SEE "https://stats.stackexchange.com/questions/483185/determine-how-good-an-auc-is-area-under-the-curve-of-roc"

###AT FIRST IMPRESSION, A READER WOULD BELIEVE THAT THE IMPROVEMENT IS MARGINAL, AND BOTH APPROACHES ARE EXCELLENT (> 90%)
###HOWEVER, THE AI-ASSISTED EXAMS SIGNIFICANTLY INCREASED THE PPV, NOT JUST SPECIFICITY
###AS A RESULT, THE 0.962 AUROC IS UNDER-LEVELED SINCE IT ONLY ACCOUNTS FOR THE SPECIFICITY IMPROVEMENT AND IGNORES PPV

###MEAN P4 WITHOUT AI ASSISTANCE

#P4 = 0.570 (0.521, 0.619)

###MEAN P4 WITH AI ASSISTANCE

#P4 = 0.687 (0.634, 0.740)

###MEAN G4 WITHOUT AI ASSISTANCE

#G4 = 0.666 (0.626, 0.704)

###MEAN G4 WITH AI ASSISTANCE

#G4 = 0.746 (0.706, 0.787)

###MEAN MCC WITHOUT AI ASSISTANCE

#RAW MCC = 0.432 (0.370, 0.491)
#TRANSFORMED MCC = 0.716 (0.685, 0.745)

###MEAN MCC WITH AI ASSISTANCE

#RAW MCC = 0.551 (0.485, 0.615)
#TRANSFORMED MCC = 0.775 (0.742, 0.807)


###DIFFERENCE BETWEEN AIDED AND UNAIDED

#AUROC DIFF = 0.038 (0.028, 0.052)

##P4 DIFF = 0.116 (0.045, 0.188)
##G4 DIFF = 0.080 (0.024, 0.136)
#RAW MCC DIFF = 0.119 (0.030, 0.208)
#TRANSFORMED MCC DIFF = 0.059 (0.015, 0.104)


###FOR THE DIFFERENCE SCORE, A CLUSTERED, STRATIFIED, NON-PARAMETRIC BOOTSTRAP APPROACH IS RECOMMENDED
###HOWEVER, DUE TO LIMITED INFORMATION FROM THE MANUSCRIPT, ONLY CLUSTERING AT THE READER LEVEL CAN BE PERFORMED (MODALITY WILL BE NESTED INSTEAD OF CROSSED)
###DESPITE THE CONSERVATIVE APPROACH, THE DIFFERENCES IN METRICS ARE STILL SIGNIFICANT IN FAVOR OF THE AI









