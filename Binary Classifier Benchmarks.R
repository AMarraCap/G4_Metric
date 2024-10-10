###COMPARE AUROC TO THE BALANCED METRIC FAMILY, ACCURACY, BALANCED ACCURACY, POPANA, KAPPA, ETC. ACROSS DIFFERENT PREVALENCE RATES

library(data.table)
library(performance)
library(ggplot2)
library(cowplot)
library(reshape2)
library(caret)
library(pROC)
library(MASS)
library(parallel)
library(doParallel)
library(splines)

###SIMULATE SYNTHETIC DATA GIVEN AN AUROC SCORE
###FOR MORE INFORMATION, SEE "https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc"

auroc.sim <- seq(0.50, 0.99, 0.01)
#auroc.sim <- c(0.80, 0.90)

prevalence.sim = 0.50
bias.sim = 0.50
Repeats = 250

rock = NULL
prev.sample = NULL
bias.sample = NULL
accuracy = NULL
sens = NULL
specif = NULL
ppv = NULL
npv = NULL
balanced.accuracy = NULL
f1.pos = NULL
f1.neg = NULL
p4 = NULL
g4 = NULL
matthews = NULL
matthews.scaled = NULL
fowlkes.mallows = NULL
youden = NULL
jaccard = NULL
kappa = NULL
kappa.scaled = NULL
gwet = NULL
gwet.scaled = NULL
brennan.prediger = NULL
brennan.prediger.scaled = NULL

rock.sim = NULL
prev.sample.sim = NULL
bias.sample.sim = NULL
accuracy.sim = NULL
sens.sim = NULL
specif.sim = NULL
ppv.sim = NULL
npv.sim = NULL
balanced.accuracy.sim = NULL
f1.pos.sim = NULL
f1.neg.sim = NULL
p4.sim = NULL
g4.sim = NULL
matthews.sim = NULL
matthews.scaled.sim = NULL
fowlkes.mallows.sim = NULL
youden.sim = NULL
jaccard.sim = NULL
kappa.sim = NULL
kappa.scaled.sim = NULL
gwet.sim = NULL
gwet.scaled.sim = NULL
brennan.prediger.sim = NULL
brennan.prediger.scaled.sim = NULL

start.time <- Sys.time()
start.time

###TIME INTERVAL FOR CHECKING PROGRESS (IN SECONDS)

time_interval <- 10

set.seed(11723)

for (j in 1:length(auroc.sim)) {

###CHECK IF THE TIME INTERVAL HAS PASSED

current_time <- Sys.time()
time_elapsed <- as.numeric(difftime(current_time, start.time, units = "secs"))

###IF THE TIME INTERVAL HAS PASSED, PRINT THE PROGRESS

if (time_elapsed > time_interval) {

cat("Current iteration:", j, "\n")
time_interval <- time_interval + 10  # INCREMENT TIME INTERVAL FOR THE NEXT CHECK

}

for (i in 1:Repeats) {
	
auroc = auroc.sim[j]

t <- sqrt(log(1/(1 - auroc)^2))
z <- t - ((2.515517 + 0.802853*t + 0.0103328*t^2) / 
         (1 + 1.432788*t + 0.189269*t^2 + 0.001308*t^3))
d <- z*sqrt(2)

n <- 10000

prediction <- c(rnorm(round(n*(1 - bias.sim), 0), mean = 0), rnorm(round(n*bias.sim, 0), mean = d))

###ADJUST THE BALANCE OF THE DATA BELOW
###DEFAULT IS A PERFECTLY BALANCED DATASET (I.E., PREVALENCE = BIAS = 50%)
###IMBALANCED DATA IS DEFINED AS HAVING A RATIO OF 70%/30% OR GREATER FOR THE TWO CLASSES

truth <- c(rep(0, round(n*(1 - prevalence.sim), 0)), rep(1, round(n*prevalence.sim, 0)))
if (length(truth) != 10000) {
truth = c(0, truth)
}
truth <- factor(truth)

rock[i] = roc(truth, prediction, direction = "<", levels = c("0", "1"))$auc

normalized <- function(x) { (x - min(x)) / (max(x) - min(x)) }
prediction <- normalized(prediction)

prediction <- factor(ifelse(prediction >= 0.500, 1, 0))
matrix = confusionMatrix(prediction, truth, positive = "1")$table

tp = matrix[4]
tn = matrix[1]
fp = matrix[2]
fn = matrix[3]

N = tp + tn + fp + fn
prev.sample[i] = (tp + fn) / N
bias.sample[i] = (tp + fp) / N

smooth = 1e-06

accuracy[i] = (tp + tn) / (tp + fp + tn + fn)
sens[i] = tp / (tp + fn)
specif[i] = tn / (tn + fp)
ppv[i] = tp / (tp + fp)
npv[i] = tn / (tn + fn)
fnr = 1 - sens[i]
fpr = 1 - specif[i]
fdr = 1 - ppv[i]
fomitr = 1 - npv[i]
balanced.accuracy[i] = (sens[i] + specif[i]) / 2

recall = sens[i]
precision = ppv[i]
#f1 = (2 * ppv[i] * sens[i]) / (ppv[i] + sens[i])
f1.pos[i] = (2 * tp) / ((2 * tp) + fn + fp + smooth)
dsc = f1.pos
f1.neg[i] = (2 * tn) / ((2 * tn) + fn + fp + smooth)

#p4[i] = 4 / (1/precision + 1/recall + 1/specif[i] + 1/npv[i])
p4[i] = (4 * tp * tn) / ((4 * tp * tn) + (tp * fp) + (tp * fn) + (fp * tn) + (fn * tn))
#g4[i] = sqrt(sens[i] * specif[i]) * ((prev.sample[i] * (1 - prev.sample[i])) / (bias.sample[i] * (1 - bias.sample[i]))) ^ (1/4)
g4[i] = sqrt((tp * tn) / ((tp + fn) * (tn + fp))) * ((tp + fn) * (tn + fp) / ((tp + fp) * (tn + fn))) ^ (1/4)
#g4[i] = (precision * recall * specif[i] * npv[i]) ^ (1/4)
matthews[i] = sqrt(sens[i] * specif[i] * ppv[i] * npv[i]) - sqrt(fnr * fpr * fdr * fomitr)
matthews.scaled[i] = (matthews[i] + 1) / 2

fowlkes.mallows[i] = sqrt(ppv[i] * sens[i])
youden[i] = sens[i] + specif[i] - 1
jaccard[i] = tp / (tp + fn + fp + smooth)
iou = jaccard[i]

kappa[i] = (2 * (tp * tn - fp * fn)) / ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn))
kappa.scaled[i] = (kappa[i] + 1) / 2
agreement = accuracy[i]
chance = 2 * ((prev.sample[i] + bias.sample[i]) / 2) * (1 - ((prev.sample[i] + bias.sample[i]) / 2))
gwet[i] = (agreement - chance) / (1 - chance)
gwet.scaled[i] = (gwet[i] + 1) / 2
k = 2
chance = 1 / k
brennan.prediger[i] = (agreement - chance) / (1 - chance)
brennan.prediger.scaled[i] = (brennan.prediger[i] + 1) / 2

}

rock.sim[j] = mean(rock)
prev.sample.sim[j] = mean(prev.sample)
bias.sample.sim[j] = mean(bias.sample)
accuracy.sim[j] = mean(accuracy)
sens.sim[j] = mean(sens)
specif.sim[j] = mean(specif)
ppv.sim[j] = mean(ppv)
npv.sim[j] = mean(npv)
balanced.accuracy.sim[j] = mean(balanced.accuracy)
f1.pos.sim[j] = mean(f1.pos)
f1.neg.sim[j] = mean(f1.neg)
p4.sim[j] = mean(p4)
g4.sim[j] = mean(g4)
matthews.sim[j] = mean(matthews)
matthews.scaled.sim[j] = mean(matthews.scaled)
fowlkes.mallows.sim[j] = mean(fowlkes.mallows)
youden.sim[j] = mean(youden)
jaccard.sim[j] = mean(jaccard)
kappa.sim[j] = mean(kappa)
kappa.scaled.sim[j] = mean(kappa.scaled)
gwet.sim[j] = mean(gwet)
gwet.scaled.sim[j] = mean(gwet.scaled)
brennan.prediger.sim[j] = mean(brennan.prediger)
brennan.prediger.scaled.sim[j] = mean(brennan.prediger.scaled)

}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

###AVERAGE SIMULATION TIME FOR EACH PREVALENCE RATE IS 250 SECONDS (APPROXIMATELY 4 MINUTES)

data = data.table(AUROC = rock.sim, Prevalence = prev.sample.sim, Bias = bias.sample.sim,
	 Sensitivity = sens.sim, Specificity = specif.sim, PPV = ppv.sim, NPV = npv.sim,
	 Accuracy = accuracy.sim, Balanced.Accuracy = balanced.accuracy.sim, F1.Positive = f1.pos.sim,
	 F1.Negative = f1.neg.sim, P4 = p4.sim, G4 = g4.sim, MCC = matthews.sim, MCC.Scaled = matthews.scaled.sim,
	 Fowlkes.Mallows = fowlkes.mallows.sim, Youden = youden.sim, Jaccard = jaccard.sim, Kappa = kappa.sim,
	 Kappa.Scaled = kappa.scaled.sim, Gwet = gwet.sim, Gwet.Scaled = gwet.scaled.sim,
	 Brennan.Prediger = brennan.prediger.sim, Brennan.Prediger.Scaled = brennan.prediger.scaled.sim)

data
data[AUROC >= 0.8, ]









