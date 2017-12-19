#Data Analysis for PERAs - tastants applied to labellar sensilla.

x <- c("sciplot", "lme4", "car", "effects", "bbmle", "ggplot2", "gridExtra") 

lapply(x, FUN = function(X) {
  do.call("library", list(X)) 
})

# Megan uses this
setwd("~/Dropbox/Megan/Drosophila_work/Data/pera")

# Ian uses this
#setwd("~/Dropbox/Dworkin_lab/Megan/Drosophila_work/Data/pera")

#Loading data sets
NaCl_data <- read.table("PERA_NaCl_edited.txt", header = T)
NaCl_data <- data.frame(NaCl_data)

NaCl_barom <- read.table("NaCl_barom_readings.txt", header = T)
NaCl_barom$delta_pres <- (NaCl_barom$pressure_6am - NaCl_barom$pressure_9am)

KCl_data <- read.table("PERA_KCl.txt", header = T)
KCl_data <- data.frame(KCl_data)

KCl_barom <- read.table("KCl_barom_readings.txt", header = T)
KCl_barom$delta_pres <- (KCl_barom$pressure_6am - KCl_barom$pressure_9am)

#Checking for duplicate rows in each dataset
sum(duplicated(NaCl_data))
sum(duplicated(KCl_data))

#merging barometric pressure into data sets.
NaCl_data_pres <- merge(NaCl_data, NaCl_barom, by="date_test")
KCl_data_pres <- merge(KCl_data, KCl_barom, by="date_test")

#Getting average response per fly to each stimulus
#NaCl
NaCl_data_pres[24:25, "avg_h2o"] <- NA
NaCl_data_pres[25:26, "avg_sugar"] <- NA
NaCl_data_pres[26:27, "avg_salt"] <- NA

NaCl_data_pres <- transform(NaCl_data_pres, avg_h2o = rowMeans(NaCl_data[, c(7, 10, 13)], na.rm = TRUE))
NaCl_data_pres <- transform(NaCl_data_pres, avg_sugar = rowMeans(NaCl_data[, c(8,11)], na.rm = TRUE))
NaCl_data_pres <- transform(NaCl_data_pres, avg_salt = rowMeans(NaCl_data[, c(9,12)], na.rm = TRUE))

#KCl
KCl_data_pres[24:25, "avg_h2o"] <- NA
KCl_data_pres[25:26, "avg_sugar"] <- NA
KCl_data_pres[26:27, "avg_salt"] <- NA

KCl_data_pres <- transform(KCl_data_pres, avg_h2o = rowMeans(KCl_data[, c(7, 10, 13)], na.rm = TRUE))
KCl_data_pres <- transform(KCl_data_pres, avg_sugar = rowMeans(KCl_data[, c(8,11)], na.rm = TRUE))
KCl_data_pres <- transform(KCl_data_pres, avg_salt = rowMeans(KCl_data[, c(9,12)], na.rm = TRUE))


###First looking at avg responses by individuals tested using labellar receptors
probos_N <- subset(NaCl_data_pres, receptor == "prob")  #should look at tarsi too, but in later analysis.
probos_K <- subset(KCl_data_pres, receptor == "prob")

##getting some preliminary plots of average responses.
probos_N$date_test <- as.numeric(probos_N$date_test)#change this back to factor for analysis
probos_K$date_test <- as.numeric(probos_K$date_test)#change this back to factor for analysis

#plotting avg responses by date.
plot(jitter(probos_N$avg_h2o) ~ probos_N$date_test) 
plot(jitter(probos_K$avg_h2o) ~ probos_K$date_test)

plot(jitter(probos_N$avg_sugar) ~ probos_N$date_test) 
plot(jitter(probos_K$avg_sugar) ~ probos_K$date_test)

plot(jitter(probos_N$avg_salt) ~ probos_N$date_test) 
plot(jitter(probos_K$avg_salt) ~ probos_K$date_test)

#Does fly order matter?
plot(jitter(probos_N$avg_h2o) ~ jitter(probos_N$Fly)) 
plot(jitter(probos_K$avg_h2o) ~ jitter(probos_K$Fly))

plot(jitter(probos_N$avg_sugar) ~ jitter(probos_N$Fly)) 
plot(jitter(probos_K$avg_sugar) ~ jitter(probos_K$Fly))

plot(jitter(probos_N$avg_salt) ~ jitter(probos_N$Fly)) 
plot(jitter(probos_K$avg_salt) ~ jitter(probos_K$Fly))

#Does the order of tastant presentation matter?

#### this function just looks at the proportion in each class. For the sugar before salt..

PropCount <- function(x) {
  crap <- table(x)
  return(c(crap/sum(crap), sample_size  = sum(crap) ))
}

#h2o
with(probos_N, tapply(avg_h2o,Sugar_before_salt, PropCount))
with(probos_K, tapply(avg_h2o,Sugar_before_salt, PropCount))

#sugar
with(probos_N, tapply(avg_sugar,Sugar_before_salt, PropCount))
with(probos_K, tapply(avg_sugar,Sugar_before_salt, PropCount))

#salt
with(probos_N, tapply(avg_salt,Sugar_before_salt, PropCount))
with(probos_K, tapply(avg_salt,Sugar_before_salt, PropCount))

#####plotting mean responses per background per allele with  95% CIs

#water response

tiff("mean_h2o_bybackground&geno2.tif", unit = "px", height = 800, width = 600)
par(mfrow=c(2,1))
lineplot.CI(Allele, avg_h2o, group = Background, data = probos_N, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_h2o", ylim = c(0, 0.3), cex.lab = 1.5, 
                     col = c("blue",  "red"), 
                     pch = c(16,16,16,16,16), main="NaCl",
                     ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x),(1.96*se(x))))

lineplot.CI(Allele, avg_h2o, group = Background, data = probos_K, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_h2o", ylim = c(0, 0.3), cex.lab = 1.5, 
                     col = c("blue", "red"), pch = c(16,16,16,16), main = "KCl",
                     ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

dev.off()

#sugar response
tiff("mean_sugar_bybackground&geno2.tif", unit = "px", height = 800, width = 600)
par(mfrow = c(2,1))
lineplot.CI(Allele, avg_sugar, group = Background, data = probos_N, cex = 1.5, xlab = "Allele", 
            ylab = "Probability of Sugar Response", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue",  "red"), 
            pch = c(16,16,16,16,16), main = "NaCl",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x)))) 

lineplot.CI(Allele, avg_sugar, group = Background, data = probos_K, cex = 1.5, xlab = "Allele", 
            ylab = "Probability of Sugar Response", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue", "red"), 
            pch = c(16,16,16,16), main = "KCl",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

dev.off()


#salt response
tiff("mean_salt_bybackground&geno.tif", unit = "px", height = 800, width = 600)
par(mfrow = c(2,1))
lineplot.CI(Background, avg_salt, group = Allele, data = probos_N, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 0.6), cex.lab = 1.5, 
            col = c("blue",  "red"), 
            pch = c(16,16,16,16,16), main = "NaCl",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_salt, group = Allele, data = probos_K, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 0.6), cex.lab = 1.5, 
            col = c("blue", "red"), 
            pch = c(16,16,16,16), main = "KCl",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))
dev.off()


####Now subsetting by salt concentration, then plotting means (+/-95% CIs) for each subset.

NaCl100 <- subset(probos_N, salt_conc == 100)
NaCl500 <- subset(probos_N, salt_conc == 500)

KCl100 <- subset(probos_K, salt_conc == 100)
KCl500 <- subset(probos_K, salt_conc == 500)

#plotting h2o
tiff("mean_h2o_bybackground&geno&saltconc.tif", unit = "px", height = 800, width = 800)
plot.new()
par(mfrow = c(2,2))
lineplot.CI(Background, avg_h2o, group = Allele, data = NaCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_h2o", ylim = c(0, 0.4), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_h2o, group = Allele, data = NaCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_h2o", ylim = c(0, 0.4), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_h2o, group = Allele, data = KCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_h2o", ylim = c(0, 0.4), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_h2o, group = Allele, data = KCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_h2o", ylim = c(0, 0.4), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))
dev.off()


#plotting sugar responses
tiff("mean_sugar_bybackground&geno&saltconc.tif", unit = "px", height = 800, width = 800)
plot.new()
par(mfrow = c(2,2))
lineplot.CI(Background, avg_sugar, group = Allele, data = NaCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_sugar", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_sugar, group = Allele, data = NaCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_sugar", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_sugar, group = Allele, data = KCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_sugar", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_sugar, group = Allele, data = KCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_sugar", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))
dev.off()

#plotting salt responses
tiff("mean_salt_bybackground&geno&saltconc.tif", unit = "px", height = 800, width = 800)
plot.new()
par(mfrow = c(2,2))
lineplot.CI(Background, avg_salt, group = Allele, data = NaCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_salt, group = Allele, data = NaCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("orange", "blue",  "red", "green", "black"), pch = c(16,16,16,16,16),
            main = "NaCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_salt, group = Allele, data = KCl100, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_100mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))

lineplot.CI(Background, avg_salt, group = Allele, data = KCl500, cex = 1.5, xlab = "Background", 
            ylab = "avg_salt", ylim = c(0, 1), cex.lab = 1.5, 
            col = c("blue", "red", "green", "black"), pch = c(16,16,16,16),
            main = "KCl_500mM",
            ci.fun= function(x) c(mean(x)-(1.96*(se(x))), mean(x)+(1.96*se(x))))
dev.off()

#Subsetting just a little further to get the observed probabilities of response
NaCl100_ore <- subset(NaCl100, Background == "ORE")
NaCl100_sam <- subset(NaCl100, Background == "SAM")

NaCl500_ore <- subset(NaCl500, Background == "ORE")
NaCl500_sam <- subset(NaCl500, Background == "SAM")

KCl100_ore <- subset(KCl100, Background == "ORE")
KCl100_sam <- subset(KCl100, Background == "SAM")

KCl500_ore <- subset(KCl500, Background == "ORE")
KCl500_sam <- subset(KCl500, Background == "SAM")

#Getting observed probability of response and CIs

#NaCl100_ore
mean.prob.switch_NaCl100_ore <- tapply(NaCl100_ore$avg_sugar, NaCl100_ore$Allele, mean)

boot.fn.lower <- function(x=NaCl100_ore$avg_sugar, N=5000) {
  lower.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  lower.CI <- quantile(lower.1, probs=c(0.025))
  lower.CI
}

LCI.obs_NaCl100_ore <- tapply(NaCl100_ore$avg_sugar, NaCl100_ore$Allele, boot.fn.lower)

boot.fn.upper <- function(x=NaCl100_ore$avg_sugar, N=5000) {
  upper.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  upper.CI <- quantile(upper.1, probs=c(0.975))
  upper.CI
}

UCI.obs_NaCl100_ore <- tapply(NaCl100_ore$avg_sugar, NaCl100_ore$Allele, boot.fn.upper)

mean.prob.switch_NaCl100_ore
LCI.obs_NaCl100_ore
UCI.obs_NaCl100_ore

#NaCl100_sam
#Getting observed probability of response and CIs
mean.prob.switch_NaCl100_sam <- tapply(NaCl100_sam$avg_sugar, NaCl100_sam$Allele, mean)

boot.fn.lower <- function(x=NaCl100_sam$avg_sugar, N=5000) {
  lower.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  lower.CI <- quantile(lower.1, probs=c(0.025))
  lower.CI
}

LCI.obs_NaCl100_sam <- tapply(NaCl100_sam$avg_sugar, NaCl100_sam$Allele, boot.fn.lower)

boot.fn.upper <- function(x=NaCl100_sam$avg_sugar, N=5000) {
  upper.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  upper.CI <- quantile(upper.1, probs=c(0.975))
  upper.CI
}

UCI.obs_NaCl100_sam <- tapply(NaCl100_sam$avg_sugar, NaCl100_sam$Allele, boot.fn.upper)

mean.prob.switch_NaCl100_sam
LCI.obs_NaCl100_sam
UCI.obs_NaCl100_sam


###Now for some modeling......

probos_N$subject <- with(probos_N, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))
probos_N$date_test <- as.factor(probos_N$date_test)

#reshaping big data sets into long format
NaCl_reshaped <- reshape(probos_N, 
                               varying = list(7:13), 
                               direction="long", 
                               idvar="subject", 
                               v.names="PER")

NaCl_reshaped <- NaCl_reshaped[order(NaCl_reshaped$subject),]
probos_N <- probos_N[order(probos_N$subject),]

head(NaCl_reshaped) #checking to be sure that reshape worked.
head(probos_N)

tail(NaCl_reshaped)
tail(probos_N)


probos_K$subject <- with(probos_K, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))
probos_K$date_test <- as.factor(probos_K$date_test)

KCl_reshaped <- reshape(probos_K, 
                         varying = list(7:13), 
                         direction="long", 
                         idvar="subject", 
                         v.names="PER")

KCl_reshaped <- KCl_reshaped[order(KCl_reshaped$subject),]
probos_K <- probos_K[order(probos_K$subject),]

head(KCl_reshaped) #checking to be sure that reshape worked.
head(probos_K)

tail(KCl_reshaped)
tail(probos_K)

#now I need to add in columns for treatment
NaCl_reshaped[23:24, "Treatment"] <- NA

NaCl_reshaped$Treatment <- NaCl_reshaped$time
NaCl_reshaped$Treatment <- ifelse(NaCl_reshaped$Treatment == 2, "sugar", ifelse(NaCl_reshaped$Treatment == 5, "sugar", 
                          ifelse(NaCl_reshaped$Treatment == 3, "salt", ifelse(NaCl_reshaped$Treatment == 6, "salt", "h2o"))))
str(NaCl_reshaped)

NaCl_reshaped$Treatment <- as.factor(NaCl_reshaped$Treatment)

#modifying time to reflect true order of presentation
NaCl_sugfirst <- subset(NaCl_reshaped, Sugar_before_salt == "y")
NaCl_saltfirst <- subset(NaCl_reshaped, Sugar_before_salt == "n")

NaCl_saltfirst$time[NaCl_saltfirst$time == 2] <- 3.1
NaCl_saltfirst$time[NaCl_saltfirst$time == 3] <- 2
NaCl_saltfirst$time[NaCl_saltfirst$time == 3.1] <- 3

NaCl_saltfirst$time[NaCl_saltfirst$time == 5] <- 6.1
NaCl_saltfirst$time[NaCl_saltfirst$time == 6] <- 5
NaCl_saltfirst$time[NaCl_saltfirst$time == 6.1] <- 6

NaCl_reshaped <- rbind(NaCl_sugfirst, NaCl_saltfirst)


KCl_reshaped[23:24, "Treatment"] <- NA

KCl_reshaped$Treatment <- KCl_reshaped$time
KCl_reshaped$Treatment <- ifelse(KCl_reshaped$Treatment == 2, "sugar", ifelse(KCl_reshaped$Treatment == 5, "sugar", 
                                                                                ifelse(KCl_reshaped$Treatment == 3, "salt", ifelse(KCl_reshaped$Treatment == 6, "salt", "h2o"))))
str(KCl_reshaped)

KCl_reshaped$Treatment <- as.factor(KCl_reshaped$Treatment)

KCl_sugfirst <- subset(KCl_reshaped, Sugar_before_salt == "y")
KCl_saltfirst <- subset(KCl_reshaped, Sugar_before_salt == "n")

KCl_saltfirst$time[KCl_saltfirst$time == 2] <- 3.1
KCl_saltfirst$time[KCl_saltfirst$time == 3] <- 2
KCl_saltfirst$time[KCl_saltfirst$time == 3.1] <- 3

KCl_saltfirst$time[KCl_saltfirst$time == 5] <- 6.1
KCl_saltfirst$time[KCl_saltfirst$time == 6] <- 5
KCl_saltfirst$time[KCl_saltfirst$time == 6.1] <- 6

KCl_reshaped <- rbind(KCl_sugfirst, KCl_saltfirst)


#modeling NaCl data first.

NaCl_reshaped$Allele <- relevel(NaCl_reshaped$Allele, ref = "wt")

model_trial_1 <- glm(PER ~ Background*Allele*Sex*Treatment + salt_conc + date_test + 
                       Sugar_before_salt + time + delta_pres + avg_salt + avg_h2o, 
                     data=NaCl_reshaped, family=binomial, x=T)
summary(model_trial_1)

###problems start here - model 2 fails to converge
#Warning messages:
#1: Some predictor variables are on very different scales: consider rescaling 
#2: In (function (fn, par, lower = rep.int(-Inf, n), upper = rep.int(Inf,  :failure to converge in 10000 evaluations
#3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :unable to evaluate scaled gradient
#4: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

NaCl_reshaped$salt_conc <- as.factor(NaCl_reshaped$salt_conc) #solved first convergence issue by making salt_conc as factor

model_trial_2 <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         Sugar_before_salt + time + delta_pres + (1 + time + Treatment|subject) + 
                         (1|date_test), data=NaCl_reshaped, family=binomial)

#Tried to get rid of non-convergence warnings by scaling time & rerunning model2. Did not work.
NaCl_reshaped$time_sc <- scale(NaCl_reshaped$time)
NaCl_reshaped$delta_pres_sc <- scale(NaCl_reshaped$delta_pres)
model_trial_2_scaled <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         Sugar_before_salt + time_sc + delta_pres_sc + (1 + time_sc + Treatment|subject) + 
                         (1|date_test), data=NaCl_reshaped, family=binomial)

#Changed the optimizer and number of iterations.
model_trial_2_scaled <- update(model_trial_2_scaled,control=glmerControl(optimizer="bobyqa"))

summary(model_trial_2_scaled)
car::Anova(model_trial_2_scaled)

#Maybe time shouldn't be considered a random effect per this:
#https://stats.stackexchange.com/questions/110004/how-scared-should-we-be-about-convergence-warnings-in-lme4




model_trial_3 <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                              Sugar_before_salt + time + delta_pres + (1 + Treatment|subject) + 
                              (1|date_test), data=NaCl_reshaped, 
                            family=binomial)

summary(model_trial_3_scaled)#this is model1 in my manuscript, also gives convergence warnings 


#getting correlations between random effects
VarCorr(model_trial_3)

#outputting model summary to a file.
#sink(file="model_trial_3.txt")
#summary(model_trial_3)
#sink(NULL)

car::Anova(model_trial_3)


model_trial_4 <- glmer(PER ~ (Background + Allele + Treatment)^2 + Sex + salt_conc + Sugar_before_salt + time + delta_pres + 
                         (1 + time + Treatment|subject) + (1|date_test), data=NaCl_reshaped, family=binomial)

summary(model_trial_4)#model 2 in manuscript
car::Anova(model_trial_4)

model_trial_5 <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + Sugar_before_salt + time + delta_pres + 
                         (1 + Treatment|subject) + (1|date_test), data=NaCl_reshaped, family=binomial)

summary(model_trial_5)#model 3 in manuscript
car::Anova(model_trial_5)

model_trial_6 <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + Sugar_before_salt + time +  
                         (1 + Treatment|subject) + (1|date_test), data=NaCl_reshaped, family=binomial)

summary(model_trial_6)#model 4 in manuscript
car::Anova(model_trial_6)

model_trial_7 <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + time +  
                         (1 + Treatment|subject) + (1|date_test), data=NaCl_reshaped, family=binomial)

summary(model_trial_7)#model 5 in manuscript
car::Anova(model_trial_7)

model_trial_8 <- glmer(PER ~ Background + Allele + Treatment + Sex + time +  
                         (1 + Treatment|subject) + (1|date_test), data=NaCl_reshaped, family=binomial)

summary(model_trial_8)#model 6 in manuscript
car::Anova(model_trial_8)

AIC(model_trial_2,model_trial_3,model_trial_4, model_trial_5, model_trial_6, model_trial_7, model_trial_8)
BIC(model_trial_2,model_trial_3,model_trial_4, model_trial_5, model_trial_6, model_trial_7, model_trial_8)


#reduced for checking signifcance of delta_pres
model_trial_9 <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         Sugar_before_salt + time + (1 + Treatment|subject) + 
                         (1|date_test), data=NaCl_reshaped, 
                       family=binomial)
summary(model_trial_9)
car::Anova(model_trial_9)

#reduced for checking the significance of tastant order
model_trial_10 <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         time + delta_pres + (1 + Treatment|subject) + 
                         (1|date_test), data=NaCl_reshaped, 
                       family=binomial)
summary(model_trial_10)
car::Anova(model_trial_10)

#reduced for checking significance of salt_conc
model_trial_11 <- glmer(PER ~ (Background*Allele*Treatment) + Sex + 
                          Sugar_before_salt + time + delta_pres + (1 + Treatment|subject) + 
                          (1|date_test), data=NaCl_reshaped, 
                        family=binomial)
summary(model_trial_11)
car::Anova(model_trial_11)


#now for the pboot function
pboot <- function(simpleModel,complexModel) {
  s <- simulate(simpleModel)
  L0 <- logLik(refit(simpleModel,s))
  L1 <- logLik(refit(complexModel,s))
  2*(L1-L0)
}


#looking at significance of 3rd order interaction first.
obsdev <- c(2*(logLik(model_trial_3) - logLik(model_trial_4)))

set.seed(101)
third_order<- replicate(1000, 
    pboot(complexModel = model_trial_3, simpleModel=model_trial_4))

mean(third_order>obsdev)

###p-value is 0.03 @seed = 1001, for 100 replicate simulations - need to do 1000.  

#now looking at significance of pressure.
obsdev <- c(2*(logLik(model_trial_3) - logLik(model_trial_9)))

set.seed(1001)
salt <- replicate(1000, 
                        pboot(complexModel = model_trial_3, simpleModel=model_trial_9))

mean(salt>obsdev)

#now looking at significance of tastant order.
obsdev <- c(2*(logLik(model_trial_3) - logLik(model_trial_10)))

set.seed(1001)
TasteOrder <- replicate(1000, 
                  pboot(complexModel = model_trial_3, simpleModel=model_trial_10))

mean(TasteOrder>obsdev)

#now looking at significance of salt_conc.
obsdev <- c(2*(logLik(model_trial_3) - logLik(model_trial_11)))

set.seed(1001)
SaltConc <- replicate(1000, 
                        pboot(complexModel = model_trial_3, simpleModel=model_trial_11))

mean(SaltConc>obsdev)

#conflicting AIC and BIC outcomes, so am going to go with model_trial_3.

effects_model_N1 <- Effect(c("Treatment", "Background", "Allele", "salt_conc"), model_trial_3)  
#plot(effects_model_N1)

prob_sub_NaCl_fit <- as.data.frame(effects_model_N1)

#plotting effects for NaCl
NaCl500 <- subset(prob_sub_NaCl_fit, salt_conc == 500)
NaCl100 <- subset(prob_sub_NaCl_fit, salt_conc == 100)

Treat_h2o_500 <- subset(NaCl500, Treatment == "h2o")
Treat_sug_500 <- subset(NaCl500, Treatment == "sugar")
Treat_salt_500 <- subset(NaCl500, Treatment == "salt")

Treat_h2o_100 <- subset(NaCl100, Treatment == "h2o")
Treat_sug_100 <- subset(NaCl100, Treatment == "sugar")
Treat_salt_100 <- subset(NaCl100, Treatment == "salt")

Treat_h2o_500$Allele <- as.numeric(ordered(Treat_h2o_500$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug_500$Allele <- as.numeric(ordered(Treat_sug_500$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt_500$Allele <- as.numeric(ordered(Treat_salt_500$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))

Treat_h2o_100$Allele <- as.numeric(ordered(Treat_h2o_100$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug_100$Allele <- as.numeric(ordered(Treat_sug_100$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt_100$Allele <- as.numeric(ordered(Treat_salt_100$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))


png(filename = "BxGxT_byconcN.png", height = 1000, width = 1000, units = "px")

plot1 <- ggplot(data = Treat_h2o_100, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3)) + scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot2 <- ggplot(data = Treat_h2o_500, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


plot3 <- ggplot(data = Treat_sug_100, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot4 <- ggplot(data = Treat_sug_500, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


plot5 <- ggplot(data = Treat_salt_100, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of NaCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot6 <- ggplot(data = Treat_salt_500, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of NaCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2)

dev.off()


effects_model_N2 <- Effect(c("Treatment","Background","Allele"), model_trial_3)  
#plot(effects_model_N2)

prob_sub_NaCl_fit2 <- as.data.frame(effects_model_N2)

Treat_h2o2 <- subset(prob_sub_NaCl_fit2, Treatment == "h2o")
Treat_sug2 <- subset(prob_sub_NaCl_fit2, Treatment == "sugar")
Treat_salt2 <- subset(prob_sub_NaCl_fit2, Treatment == "salt")

Treat_h2o2$Allele <- as.numeric(ordered(Treat_h2o2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug2$Allele <- as.numeric(ordered(Treat_sug2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt2$Allele <- as.numeric(ordered(Treat_salt2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))

#####FIGURE1
png(filename = "BxGxT_N.png", height = 1000, width = 700, units = "px")

plot1 <- ggplot(data = Treat_h2o2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) + 
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3)) + scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)

plot2 <- ggplot(data = Treat_sug2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)

plot3 <- ggplot(data = Treat_salt2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of NaCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)


grid.arrange(plot1, plot2, plot3, ncol=1)

dev.off()



effects_model_N3 <- Effect(c("Treatment","Background", "Allele", "Sex"), model_trial_3)  
plot(effects_model_N3)

effects_model_N4 <- Effect(c( "Treatment","Background", "Allele", "Sugar_before_salt"),model_trial_3)
plot(effects_model_N4)

effects_model_N5 <- Effect(c("Background", "Allele", "time"),model_trial_3)
plot(effects_model_N5)

effects_model_N6 <- Effect(c("Treatment", "salt_conc"), model_trial_3)
plot(effects_model_N6)

#using model reduction and pboot function to look at significance of model terms for NaCl responses.



#now modeling KCl data
KCl_reshaped$Allele <- relevel(KCl_reshaped$Allele, ref = "wt")
model_trial_1K <- glm(PER ~ Background*Allele*Sex*Treatment + salt_conc + date_test + 
                       Sugar_before_salt + time + delta_pres + avg_salt + avg_h2o, 
                     data=KCl_reshaped, family=binomial, x=T)
summary(model_trial_1K)

model_trial_2K <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         Sugar_before_salt + time + delta_pres + (1 + time + Treatment|subject) + 
                         (1|date_test), data=KCl_reshaped, 
                       family=binomial)
summary(model_trial_2K)
car::Anova(model_trial_2K)

model_trial_3K <- glmer(PER ~ (Background*Allele*Treatment) + Sex + salt_conc + 
                         Sugar_before_salt + time + delta_pres + (1 + Treatment|subject) + 
                         (1|date_test), data=KCl_reshaped, 
                       family=binomial)
summary(model_trial_3K)

sink(file="model_trial_3K.txt")
summary(model_trial_3K)
sink(NULL)

car::Anova(model_trial_3K)


model_trial_4K <- glmer(PER ~ (Background + Allele + Treatment)^2 + Sex + salt_conc + Sugar_before_salt + time + delta_pres + 
                         (1 + Treatment|subject) + (1|date_test), data=KCl_reshaped, family=binomial)

summary(model_trial_4K)
car::Anova(model_trial_4K)

model_trial_5K <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + Sugar_before_salt + time + delta_pres + (1 + Treatment|subject) + (1|date_test), 
                               data=KCl_reshaped, family=binomial)

summary(model_trial_5K)#this is model 3 in manuscript
car::Anova(model_trial_5K)

model_trial_6K <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + Sugar_before_salt + time +  
                          (1 + Treatment|subject) + (1|date_test), data=KCl_reshaped, family=binomial)
                                              
summary(model_trial_6K)#this is model 4 in manuscript
car::Anova(model_trial_6K)
                        
model_trial_7K <- glmer(PER ~ Background + Allele + Treatment + Sex + salt_conc + time +  
                          (1 + Treatment|subject) + (1|date_test), data=KCl_reshaped, family=binomial)
                                                
summary(model_trial_7K)#this is model 5 in manuscript
car::Anova(model_trial_7K)
                        
model_trial_8K <- glmer(PER ~ Background + Allele + Treatment + Sex + time +  
                          (1 + Treatment|subject) + (1|date_test), data=KCl_reshaped, family=binomial)
                                                
summary(model_trial_8K)#this is model 6 in manuscript
car::Anova(model_trial_8K)
                                                
AIC(model_trial_2K,model_trial_3K,model_trial_4K, model_trial_5K, model_trial_6K, model_trial_7K, model_trial_8K)
BIC(model_trial_2K,model_trial_3K,model_trial_4K, model_trial_5K, model_trial_6K, model_trial_7K, model_trial_8K)

#model_trial_2K is also best for the KCL dataset.

effects_model_K1 <- Effect(c("Treatment", "Background", "Allele", "salt_conc"), model_trial_3K)  
#plot(effects_model_K1)

prob_sub_KCl_fit <- as.data.frame(effects_model_K1)

#plotting effects for KCl
KCl500 <- subset(prob_sub_KCl_fit, salt_conc == 500)
KCl100 <- subset(prob_sub_KCl_fit, salt_conc == 100)

Treat_h2o_500K <- subset(KCl500, Treatment == "h2o")
Treat_sug_500K <- subset(KCl500, Treatment == "sugar")
Treat_salt_500K <- subset(KCl500, Treatment == "salt")

Treat_h2o_100K <- subset(KCl100, Treatment == "h2o")
Treat_sug_100K <- subset(KCl100, Treatment == "sugar")
Treat_salt_100K <- subset(KCl100, Treatment == "salt")

Treat_h2o_500K$Allele <- as.numeric(ordered(Treat_h2o_500K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug_500K$Allele <- as.numeric(ordered(Treat_sug_500K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt_500K$Allele <- as.numeric(ordered(Treat_salt_500K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))

Treat_h2o_100K$Allele <- as.numeric(ordered(Treat_h2o_100K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug_100K$Allele <- as.numeric(ordered(Treat_sug_100K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt_100K$Allele <- as.numeric(ordered(Treat_salt_100K$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))


png(filename = "BxGxT_byconcK.png", height = 1000, width = 1000, units = "px")

plot1 <- ggplot(data = Treat_h2o_100K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3)) + scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot2 <- ggplot(data = Treat_h2o_500K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


plot3 <- ggplot(data = Treat_sug_100K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot4 <- ggplot(data = Treat_sug_500K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


plot5 <- ggplot(data = Treat_salt_100K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of KCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

plot6 <- ggplot(data = Treat_salt_500K, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of KCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)


grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2)

dev.off()


effects_model_K2 <- Effect(c("Treatment","Background", "Allele"), model_trial_3K)  
#plot(effects_model_K2)

prob_sub_KCl_fit2 <- as.data.frame(effects_model_K2)

Treat_h2o2 <- subset(prob_sub_KCl_fit2, Treatment == "h2o")
Treat_sug2 <- subset(prob_sub_KCl_fit2, Treatment == "sugar")
Treat_salt2 <- subset(prob_sub_KCl_fit2, Treatment == "salt")

Treat_h2o2$Allele <- as.numeric(ordered(Treat_h2o2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_sug2$Allele <- as.numeric(ordered(Treat_sug2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))
Treat_salt2$Allele <- as.numeric(ordered(Treat_salt2$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))

#####FIGURE2
png(filename = "BxGxT_K.png", height = 1000, width = 700, units = "px")

plot1 <- ggplot(data = Treat_h2o2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) + 
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Water Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3)) + scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)

plot2 <- ggplot(data = Treat_sug2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of Sugar Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)

plot3 <- ggplot(data = Treat_salt2, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of KCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(),axis.title.y = element_text(vjust = 0.3))+ scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[Etx4]","sd[E3]", "sd[58d]")) + ylim(0,1)


grid.arrange(plot1, plot2, plot3, ncol=1)

dev.off()



#Tarsal responses - only done with NaCl tests
tars_N <- subset(NaCl_data_pres, receptor == "tars")

salt_tars_N <- subset(tars_N, salt1 == 1)
summary(salt_tars_N)

#only done for NaCl tests, plus only 4 of the 1300 individuals tested
#responded to salt.  Therefore I am just going to leave tarsal responses
#out of paper.  



#####modeling responses to salt separately.

#just looking at salt response when stimulus was provided to proboscis
NaCl_sub_salt <- probos_N[, c(1,2,3,4,5,9,12,14,17,20,23,24,25,26)]

# Create a variable for each unique fly as an identifier
NaCl_sub_salt$subject <- with(NaCl_sub_salt, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))

NaCl_sub_salt <- NaCl_sub_salt[order(NaCl_sub_salt$subject),]
head(NaCl_sub_salt)


#here subsetting out salt responses during KCl tests
KCl_sub_salt <- probos_K[, c(1,2,3,4,5,9,12,14,17,20,23,24,25,26)]

# Create a variable for each unique fly as an identifier
KCl_sub_salt$subject <- with(KCl_sub_salt, interaction(Allele, Background, Sex, date_test, Sugar_before_salt, Fly, salt_conc, drop=T, sep="_"))

KCl_sub_salt <- KCl_sub_salt[order(KCl_sub_salt$subject),]
head(KCl_sub_salt)


# reshaping both datasets into long format
NaCl_salt_reshaped <- reshape(NaCl_sub_salt, 
                               varying = list(6:7), 
                               direction="long", 
                               idvar="subject", 
                               v.names="PER")

NaCl_salt_reshaped <- NaCl_salt_reshaped[order(NaCl_salt_reshaped$subject),]

head(NaCl_salt_reshaped) #checking to be sure that reshape worked.
head(NaCl_sub_salt)

tail(NaCl_salt_reshaped)
tail(NaCl_sub_salt)

KCl_salt_reshaped <- reshape(KCl_sub_salt, 
                              varying = list(6:7), 
                              direction="long", 
                              idvar="subject", 
                              v.names="PER")

KCl_salt_reshaped <- KCl_salt_reshaped[order(KCl_salt_reshaped$subject),]

head(KCl_salt_reshaped) #checking to be sure that reshape worked.
head(KCl_sub_salt)


# fit the basic glm just with fixed effects
#set up for modeling
NaCl_salt_reshaped$date_test <- as.factor(NaCl_salt_reshaped$date_test)
NaCl_salt_reshaped$Allele <- relevel(NaCl_salt_reshaped$Allele, "wt")

#originally considered time as numeric, but since we only have first and second assay, and 4 way interaction model failed to converge
#when time was treated as numeric, am considering it a factor.  

NaCl_salt_reshaped$time <- as.factor(NaCl_salt_reshaped$time) 

model_trial_1 <- glm(PER ~ Background*Allele*Sex + salt_conc + date_test + 
                       Sugar_before_salt + time + delta_pres + avg_salt + avg_h2o, 
                     data=NaCl_salt_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", singular.ok = T, type=3)

# Clear signal, but the date effect is large. try lmer

model_trial_2 <- glmer(PER ~ (Background*Allele) + Sex + salt_conc + 
                              Sugar_before_salt + time + (1|subject) + 
                              (1|date_test), data=NaCl_salt_reshaped, 
                            family=binomial)

summary(model_trial_2)

effects_model_N1 <- Effect(c("Background", "Allele"), model_trial_2)  
plot(effects_model_N1)

prob_sub_NaCl_fit <- as.data.frame(effects_model_N1)
prob_sub_NaCl_fit

prob_sub_NaCl_fit$Allele <- as.numeric(ordered(prob_sub_NaCl_fit$Allele, levels = c("wt", "sd1", "etx4","sde3","58d")))

#making plot
png(filename = "BxG_NaCl_full.png", height = 400, width = 600, units = "px")

plot1 <- ggplot(data = prob_sub_NaCl_fit, aes(x = Allele, y = fit, ymin = lower, ymax = upper, colour = Background)) +
  geom_point(position = position_dodge(width = 0.4), cex = 4) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
  scale_colour_manual(values = c("black", "red")) + ylab("Probability of NaCl Response") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.key = element_blank(), axis.title.y = element_text(vjust = 0.3)) + scale_x_discrete(limit=c(1,2,3,4,5),labels = c("wt","sd[1]","sd[ext4]","sd[e3]", "sd[58d]")) + ylim(0,1)

print(plot1)

dev.off()



#Fit is much poorer if we exclude time as a random effect, and differences between the significance levels of fixed effects 
#dramatically changes the significance of fixed effects.  Leaning toward including time, but doing a check with parboot to
#be sure.  

pboot <- function(simpleModel,complexModel) {
  s <- simulate(simpleModel)
  L0 <- logLik(refit(simpleModel,s))
  L1 <- logLik(refit(complexModel,s))
  2*(L1-L0)
}

#obsdev <- c(2*(logLik(model_trial_2_4int) - logLik(model_trial_3)))

#set.seed(1001)
#NaCl_sugar_PB <- replicate(100, 
#    pboot(complexModel = model_trial_2_4int, simpleModel=model_trial_3))

#mean(NaCl_sugar_PB>obsdev)

###Yes - value of 0 suggests we should allow subject response to vary by time; retain as a random effect in model.

#also had some false convergence warnings - need to get rid of time as Ian did on May 8, 2013? If so, how do we account for
#pseudoreplication?

model_trial_4 <- glmer(PER ~ (Background + Allele + Sex + salt_conc)^2 + Sugar_before_salt  + time + (1 + time|subject), 
                       data=NaCl_sugar_reshaped, family=binomial)   

summary(model_trial_4)

#Should we keep date test in model?

car::Anova(model_trial_2_4int, test.statistic="chisq")
car::Anova(model_trial_4, test.statistic="chisq")

#checking with parboot

#obsdev.2 <- c(2*(logLik(model_trial_2_4int) - logLik(model_trial_4)))
#obsdev.2

#set.seed(1001)
#NaCl_sugar2_PB <- replicate(100, 
#    pboot(simpleModel = model_trial_4, complexModel = model_trial_2_4int))

#mean(NaCl_sugar2_PB>obsdev.2)

#value = 0.01, but again with warnings.

#Also checking AIC and BIC values to look at model fits
AIC(model_trial_2_4int,model_trial_3,model_trial_4) 
BIC(model_trial_2_4int,model_trial_3,model_trial_4)

# Model fitting metrics(AIC, BIC) point to using the most complex model - double check all with parboot.

#Does tastant order matter?
model_trial_5 <- glmer(PER ~ (Background + Allele + Sex + salt_conc)^2 + time +
                         (time|subject) + (1|date_test), data=NaCl_sugar_reshaped, family=binomial)
model_trial_5

AIC(model_trial_2_4int, model_trial_5)
BIC(model_trial_2_4int, model_trial_5)

# AICc
ICtab(model_trial_2_4wayint,model_trial_3,model_trial_4,model_trial_5, type="AICc", nobs= model_trial_2_4wayint@dims["n"])

car::Anova(model_trial_2_4int, test.statistic="chisq")
car::Anova(model_trial_5, test.statistic="chisq")  # results for effects of interest are stable both as estimates and in test for deviance. Compute parametric bootstrap. 

#checking model 5 with parboot
#obsdev.3 <- c(2*(logLik(model_trial_2_4int) - logLik(model_trial_5)))
#obsdev.3

#set.seed(1001)
#NaCl_sugar3_PB <- replicate(100, 
#    pboot(simpleModel = model_trial_5, complexModel = model_trial_2_4int))

#mean(NaCl_sugar3_PB>obsdev.3)

#All metrics suggest that the full model and model 5 are roughly equal. Decide which model to estimate parameters from based upon
#results from salt models.



# Plots the interaction between background and allele on probability scale (really convenient)
# (i.e converts back using inverse logit)
effects_model_N1 <- Effect(c("Background", "salt_conc"), model_trial_2_4int)  
plot(effects_model_N1)

effects_model_N2 <- Effect(c("Allele", "salt_conc"), model_trial_2_4int)  
plot(effects_model_N2)

effects_model_N3 <- Effect(c("Background", "Allele"), model_trial_2_4int)  
plot(effects_model_N3)

effects_model_N4 <- Effect(c("Background", "Allele", "Sex"), model_trial_2_4int)  
plot(effects_model_N4)

effects_model_N5 <- Effect(c( "Allele", "salt_conc","Background"),model_trial_2_4int)
plot(effects_model_N5)

effects_model_N6 <- Effect(c("Background","Allele", "avg_sugar"), model_trial_2_4int)  
plot(effects_model_N6)

####Why is the variance so huge for these?  This is not the case at all if you look at the observed probabilities of response
####in the EDA plots?

effects_model_N7 <- Effect(c("Background", "Allele", "Sex", "time"),model_trial_2_4int)
plot(effects_model_N7)

####This is also very interesting.  Time 2 responses seem to be incredibly variable compared to 
####time1 responses.  

effects_model_N7 <- Effect(c("Background","Sugar_before_salt", "time" ),model_trial_2_4int)
plot(effects_model_N7)


# To do this by hand you need to look at 
effects_model_MF$fit
effects_model_MF$lower
effects_model_MF$upper

# But these need to be converted back to probability scale
InvLog <- function(x) {exp(x)/(1+exp(x))}
InvLog(effects_model_MF$fit) # compare this to a call of effects_model

# also there is a convenience function (author of package really thought ahead) to grab all of this
prob_sub_NaCl_fit <- as.data.frame(effects_model_N4)

effects_model2_MF <- allEffects(model_trial_6)
plot(effects_model2_MF)


#####now doing some model fitting for sugar responses during KCl tests.

#not going to run the parametric boostraps for now - am assuming that we are going to be using the 
#full model including all random effects.  

# fit the basic glm just with fixed effects
#set up for modeling
KCl_salt_reshaped$date_test <- as.factor(KCl_salt_reshaped$date_test)
KCl_salt_reshaped$Allele <- relevel(KCl_salt_reshaped$Allele, "wt")

#considering time as a factor - see reasons above.  

KCl_salt_reshaped$time <- as.factor(KCl_salt_reshaped$time) 

model_trial_1 <- glm(PER ~ Background*Allele*Sex + salt_conc + date_test + Sugar_before_salt + time, 
                     data=KCl_sugar_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3) #this will not run and gives the following error:Error in Anova.III.LR.glm(mod, singular.ok = singular.ok)

car::Anova(model_trial_1, singular.ok=T, test.statistic="LR", type = 3)

# Again, clear signal, but the date effect is large. Using lmer

model_trial_2 <- glmer(PER ~ (Background + Allele + Sex + salt_conc + avg_sugar)^2 + Sugar_before_salt + time + (1 + time|subject) + 
                              (1|date_test), data=KCl_salt_reshaped, family=binomial)
summary(model_trial_2)

#just looking at the random effect of time for KCl test sugar response.
model_trial_3 <- glmer(PER ~ (Background + Allele + Sex + salt_conc)^2 + Sugar_before_salt + time + (1|subject) + (1|date_test), 
                       data=KCl_sugar_reshaped, family=binomial) 

summary(model_trial_3)

AIC(model_trial_2, model_trial_3)
BIC(model_trial_2, model_trial_3)

#fits are about the same for these two models, but still including time as a random (and fixed) effect is a good idea to 
#account for pseudoreplication.

# quick check for the fixed effects (confirm with ParBoot or MCMC if necessary)
car::Anova(model_trial_2, test.statistic="chisq")
car::Anova(model_trial_3, test.statistic="chisq")


#Should we keep date test in model?
model_trial_4 <- glmer(PER ~ (Background + Allele + Sex + salt_conc)^2 + Sugar_before_salt  + time + (1 + time|subject), 
                       data=KCl_sugar_reshaped, family=binomial)   

summary(model_trial_4)

AIC(model_trial_2, model_trial_4)
BIC(model_trial_2, model_trial_4)

car::Anova(model_trial_2, test.statistic="chisq")
car::Anova(model_trial_4, test.statistic="chisq")

#model fits are about the same for these two models but fixed effects aren't stable - so important to retain date_test.

#Does tastant order matter?
model_trial_5 <- glmer(PER ~ (Background + Allele + Sex + salt_conc)^2 + time +
                         (time|subject) + (1|date_test), data=KCl_sugar_reshaped, family=binomial)
model_trial_5

AIC(model_trial_2, model_trial_5)
BIC(model_trial_2, model_trial_5)

# AICc
ICtab(model_trial_2,model_trial_3,model_trial_4,model_trial_5, type="AICc", nobs= model_trial_2@dims["n"])

car::Anova(model_trial_2, test.statistic="chisq")
car::Anova(model_trial_5, test.statistic="chisq")  # results for effects of interest are largely stable both as estimates and in test for deviance. Compute parametric bootstrap. 

#What do the effects look like for the full model of KCl test sugar responses?
effects_model_k1 <- Effect(c("Background", "salt_conc"), model_trial_2)
plot(effects_model_k1)

effects_model_k2 <- Effect(c("Allele","salt_conc"), model_trial_2)  
plot(effects_model_k2)

effects_model_k3 <- Effect(c("Background", "Allele", "avg_sugar"), model_trial_2)  
plot(effects_model_k3)

effects_model_k4 <- Effect(c("Background", "Allele"), model_trial_2)  
plot(effects_model_k4)

effects_model_k5 <- Effect(c("Background", "Allele", "Sex", "salt_conc"), model_trial_2)  
plot(effects_model_k5)

effects_model_k6 <- Effect(c("Background", "Allele", "Sex", "time"), model_trial_2)  
plot(effects_model_k6)

effects_model_k6 <- Effect(c("Background", "Allele"), model_trial_2)  
plot(effects_model_k6)


prob_sug_KCl_fit <- as.data.frame(effects_model_k5)

####Now looking at salt responses 

# fit the basic glm just with fixed effects
#set up for modeling
NaCl500_sugar_reshaped$date_test <- as.factor(NaCl500_sugar_reshaped$date_test)
NaCl500_sugar_reshaped$Allele <- relevel(NaCl500_sugar_reshaped$Allele, "wt")

# should time be considered a factor or not here, since we only have first and second assay.
NaCl500_sugar_reshaped$time <- as.factor(NaCl500_sugar_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=NaCl500_sugar_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3)


# Clear signal, but the date effect is large. try lmer

model_trial_2 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt + (1 + time|subject) + 
                         (1|date_test), data=NaCl500_sugar_reshaped, family=binomial)

model_trial_2    

model_trial_3 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt + (1|subject) + (1|date_test), 
                       data=NaCl500_sugar_reshaped, family=binomial) 


model_trial_3

# Do we need to allow subject's behavioral response to vary by time?

pboot <- function(m0,m1) {
  s <- simulate(m0)
  L0 <- logLik(refit(m0,s))
  L1 <- logLik(refit(m1,s))
  2*(L1-L0)
}

obsdev <- c(2*(logLik(model_trial_3) - logLik(model_trial_2)))

set.seed(1)
NaCl_sugar_PB <- replicate(100, pboot(model_trial_2, model_trial_3))

mean(NaCl_sugar_PB>obsdev)

###failing to converge, will leave in for now, but AIC and BIC show little support for allowing 
###time to vary by subject


model_trial_4 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt + (1|subject), 
                       data=NaCl500_sugar_reshaped, family=binomial)   

#Should we keep date test in model?

obsdev <- c(2*(logLik(model_trial_4) - logLik(model_trial_3)))

set.seed(1)
NaCl_sugar2_PB <- replicate(100, pboot(model_trial_3, model_trial_4))

mean(NaCl_sugar2_PB>obsdev)

#value = 0.21, don't keep date_test

#Also checking AIC and BIC values to look at model fits
AIC(model_trial_2,model_trial_3,model_trial_4 ) 
BIC(model_trial_2,model_trial_3,model_trial_4 )


# Taking date_test out of the model for now.....

model_trial_5 <- glmer(PER ~ (Background + Allele + Sex)^2 + Sugar_before_salt + 
                         (1|subject) , data=NaCl500_sugar_reshaped, family=binomial)
model_trial_5


model_trial_6 <- glmer(PER ~ (Background + Allele + Sex)^2 + (1|subject), data=NaCl500_sugar_reshaped, family=binomial) 

model_trial_6

AIC(model_trial_3,model_trial_5,model_trial_6) 
BIC(model_trial_3,model_trial_5,model_trial_6)

#parametric bootstrap for 3rd order interaction

#obsdev <- c(2*(logLik(model_trial_6) - logLik(model_trial_3)))

#set.seed(1001)
#NaCl_sugar3_PB <- replicate(100, pboot(model_trial_3, model_trial_6))

#mean(NaCl_sugar3_PB>obsdev)

#value = 0.19; no support for 3rd order interactions.

#What about for order of tastant presentation?

obsdev <- c(2*(logLik(model_trial_7) - logLik(model_trial_6)))

set.seed(1001)
NaCl_sugar4_PB <- replicate(100, pboot(model_trial_6, model_trial_7))

mean(NaCl_sugar4_PB>obsdev)

#value = 0.5, so tastant order doesn't look important for sugar response

model_trial_7
car::Anova(model_trial_7, test.statistic="chisq")   # results for effects of interest are stable both as estimates and in test for deviance. Compute parametric bootstrap. 

# Plots the interaction between background and allele on probability scale (really convenient)
# (i.e converts back using inverse logit)
effects_model_MF <- Effect(c("Background", "Allele"), model_trial_7)  
plot(effects_model_MF)
# To do this by hand you need to look at 
effects_model_MF$fit
effects_model_MF$lower
effects_model_MF$upper

# But these need to be converted back to probability scale
InvLog <- function(x) {exp(x)/(1+exp(x))}
InvLog(effects_model_MF$fit) # compare this to a call of effects_model

# also there is a convenience function (author of package really thought ahead) to grab all of this
as.data.frame(effects_model_MF)


effects_model2_MF <- allEffects(model_trial_7)
plot(effects_model2_MF)


######Now fitting a model for salt

NaCl100_salt <- NaCl100[, c(1,2,3,4,5,9,12,14)]


# reshaping into long format
NaCl100_salt_reshaped <- reshape(NaCl100_salt, 
                                  varying = list(6:7), 
                                  direction="long", 
                                  idvar="subject", 
                                  v.names="PER")

head(NaCl100_salt_reshaped)
tail(NaCl100_salt_reshaped)

# sorting data by subject    
NaCl100_salt_reshaped<- NaCl100_salt_reshaped[order(NaCl100_salt_reshaped$subject),]

# check that reshaped produced appropriate data set
head(NaCl100_salt_reshaped)    
head(NaCl100)

tail(NaCl100_salt_reshaped)
tail(NaCl100)

#seting up for modeling
NaCl100_salt_reshaped$date_test <- as.factor(NaCl100_salt_reshaped$date_test)
NaCl100_salt_reshaped$Allele <- relevel(NaCl100_salt_reshaped$Allele, "wt")

# should time be considered a factor or not here, since we only have first and second assay.
NaCl100_salt_reshaped$time <- as.factor(NaCl100_salt_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=NaCl100_salt_reshaped, family=binomial, x=T)

model_trial_1

#now fitting random effects
model_trial_2 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt + (1|date_test) 
                       + (time|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)

model_trial_2

#Do we need date_test in the model?
model_trial_3 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt
                       + (time|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)

model_trial_3 

#What about the time by subject random effect?

model_trial_4 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt
                       + (1|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)
model_trial_4

#or the fixed effect of time?
model_trial_5 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt
                       + (1|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)
model_trial_5

AIC(model_trial_1, model_trial_2, model_trial_3, model_trial_4, model_trial_5)
BIC(model_trial_1, model_trial_2, model_trial_3, model_trial_4, model_trial_5)

##These definitely provide evidence for using the more complex models
##Will probably go with model_trial_2, but this should be confirmed with LRboot

#Do we need 3rd order interactions?

model_trial_6 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + Sugar_before_salt + (1|date_test) 
                                        + (time|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)

model_trial_6

model_trial_7 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + Sugar_before_salt 
                       + (time|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)
model_trial_7


AIC(model_trial_2, model_trial_3, model_trial_6, model_trial_7)
BIC(model_trial_2, model_trial_3, model_trial_6, model_trial_7)

###All of these are very similar suggesting that date_test might be unnecessary and 3rd order
###interactions might be too.  Should check this with LRboot

##What about the order of tastant presentation?  
##Noted that it was highly significant in early models

model_trial_8 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + (1|date_test) 
                       + (time|subject), data=NaCl100_salt_reshaped, family=binomial, x=T)

model_trial_8


#LR parboot
obsdev <- c(2*(logLik(model_trial_6) - logLik(model_trial_8)))

set.seed(1001)
NaCl100_salt_PB <- replicate(10, 
                pboot(complexModel = model_trial_6, simpleModel=model_trial_8))

mean(NaCl100_salt_PB>obsdev)

# value is 0.3, but will need to increase replication.  Is this possible?

AIC(model_trial_6, model_trial_8)
BIC(model_trial_6, model_trial_8)

anova(model_trial_6, model_trial_8, test = "Chisq")

Anova(model_trial_6, test = "chisq")

Anova(model_trial_8, test = "chisq")

##The likelihood ratio test says that there is a significant difference (p = 0.001) between these two models, 
##suggesting that I should retain the tastant presentation order.


##in both models 6 and 8, effect of allele is significant.

NaCl500_salt <- NaCl500[, c(1,2,3,4,5,9,12,14)]


# reshaping into long format
NaCl500_salt_reshaped <- reshape(NaCl500_salt, 
                                 varying = list(6:7), 
                                 direction="long", 
                                 idvar="subject", 
                                 v.names="PER")

head(NaCl500_salt_reshaped)
tail(NaCl500_salt_reshaped)

# sorting data by subject    
NaCl500_salt_reshaped<- NaCl500_salt_reshaped[order(NaCl500_salt_reshaped$subject),]

# check that reshaped produced appropriate data set
head(NaCl500_salt_reshaped)    
head(NaCl500)

tail(NaCl500_salt_reshaped)
tail(NaCl500)

#seting up for modeling
NaCl500_salt_reshaped$date_test <- as.factor(NaCl500_salt_reshaped$date_test)
NaCl500_salt_reshaped$Allele <- relevel(NaCl500_salt_reshaped$Allele, "wt")


NaCl500_salt_reshaped$time <- as.factor(NaCl500_salt_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=NaCl500_salt_reshaped, family=binomial, x=T)

model_trial_1

#now fitting random effects
model_trial_2 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt + (1|date_test) 
                       + (time|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)

model_trial_2

#Do we need date_test in the model?
model_trial_3 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt
                       + (time|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)

model_trial_3 

#What about the time by subject random effect?

model_trial_4 <- glmer(PER ~ Background*Allele*Sex + time + Sugar_before_salt
                       + (1|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)
model_trial_4

#or the fixed effect of time?
model_trial_5 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt
                       + (1|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)
model_trial_5

AIC(model_trial_1, model_trial_2, model_trial_3, model_trial_4, model_trial_5)
BIC(model_trial_1, model_trial_2, model_trial_3, model_trial_4, model_trial_5)

##These definitely provide evidence for using the more complex models
##Will probably go with model_trial_2, but this should be confirmed with LRboot

#Do we need 3rd order interactions?

model_trial_6 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + Sugar_before_salt + (1|date_test) 
                       + (time|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)

model_trial_6

model_trial_7 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + Sugar_before_salt 
                       + (time|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)
model_trial_7


AIC(model_trial_2, model_trial_3, model_trial_6, model_trial_7)
BIC(model_trial_2, model_trial_3, model_trial_6, model_trial_7)

###All of these are very similar suggesting that date_test might be unnecessary and 3rd order
###interactions might be too.  Should check this with LRboot

##What about the order of tastant presentation?  
##Noted that it was highly significant in early models

model_trial_8 <- glmer(PER ~ (Background+Allele+Sex)^2 + time + (1|date_test) 
                       + (time|subject), data=NaCl500_salt_reshaped, family=binomial, x=T)

model_trial_8


#LR parboot
obsdev <- c(2*(logLik(model_trial_6) - logLik(model_trial_8)))

set.seed(1001)
NaCl500_salt_PB <- replicate(10, 
                             pboot(complexModel = model_trial_6, simpleModel=model_trial_8))

mean(NaCl500_salt_PB>obsdev)

AIC(model_trial_6, model_trial_8)
BIC(model_trial_6, model_trial_8)

anova(model_trial_6, model_trial_8, test = "Chisq")

Anova(model_trial_6, test = "chisq")
Anova(model_trial_8, test = "chisq")

##In the case of 500 mM NaCl, the difference between the 2 models is small according to AIC and BIC
##The LR test gives p-value of 0.47, suggesting that order of tastant doesn't matter when salt 
##concentration is high. This will need to be confirmed with LR parboot.

####This makes me wonder if we should account for some sort of salt concentration by tastant order effect
####in the models for the combined 100 and 500 mM dataset.

####It turns out that there is pretty substantial overlap in terms of date test for NaCl 100 and 500 tests.
####No going to run a model on the combined dataset.

#looking at sugar response first

probos_sugar <- probos[, c(1,2,3,4,5,8,11,14,17)]

probos_sugar_reshaped <- reshape(probos_sugar, 
                                 varying = list(6:7), 
                                 direction="long", 
                                 idvar="subject", 
                                 v.names="PER")

head(probos_sugar_reshaped)
tail(probos_sugar_reshaped)

# sorting data by subject    
probos_sugar_reshaped<- probos_sugar_reshaped[order(probos_sugar_reshaped$subject),]

# check that reshaped produced appropriate data set
head(probos_sugar_reshaped)    
head(probos)

tail(probos_sugar_reshaped)
tail(probos)

#set up for modeling
probos_sugar_reshaped$date_test <- as.factor(probos_sugar_reshaped$date_test)
probos_sugar_reshaped$Allele <- relevel(probos_sugar_reshaped$Allele, "wt")

###Fitting glm model with all fitted effects 
probos_sugar_reshaped$time <- as.factor(probos_sugar_reshaped$time) 

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + salt_conc + time, 
                     data=probos_sugar_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
Anova(model_trial_1, test.statistic="LR", type=2)
Anova(model_trial_1, test.statistic="LR", type=3)

# Clear signal, but the date effect is large. try lmer

#checking for date_test effect

model_trial_2 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt + salt_conc + time + (1 + time|subject) + 
                         (1|date_test), data=probos_sugar_reshaped, family=binomial)

model_trial_2    

model_trial_3 <- glmer(PER ~ Background*Allele*Sex +  Sugar_before_salt + salt_conc + time + (time|subject), 
                       data=probos_sugar_reshaped, family=binomial) 


model_trial_3

AIC(model_trial_1, model_trial_2, model_trial_3)
BIC(model_trial_1, model_trial_2, model_trial_3)

###According to AIC and BIC, not much differenc between 2 models - 
###model 2, including date)_test looks slightly better, will probably include date_test 

# quick check for the fixed effects (confirm with ParBoot or MCMC)
Anova(model_trial_2, test.statistic="chisq")

Anova(model_trial_3, test.statistic="chisq")

#### no difference in which fixed effects are significant (confirm with parboot)

pboot <- function(simpleModel,complexModel) {
  s <- simulate(simpleModel)
  L0 <- logLik(refit(simpleModel,s))
  L1 <- logLik(refit(complexModel,s))
  2*(L1-L0)
}


# Needs to be:
# obsdev <- c(2*(logLik(complex_model) - logLik(simpler_model)))

obsdev <- c(2*(logLik(model_trial_2) - logLik(model_trial_3)))

set.seed(1001)
probos_sugar_PB <- replicate(100, 
                              pboot(complexModel = model_trial_2, simpleModel=model_trial_3))

mean(probos_sugar_PB>obsdev)

#Should time by subject random effect be in model?
model_trial_4 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt  + salt_conc + time + 
                      (1|subject) + (1|date_test),data=probos_sugar_reshaped, family=binomial)   

AIC(model_trial_2, model_trial_3, model_trial_4)
BIC(model_trial_2, model_trial_3, model_trial_4)

#yes!

##now looking at fixed effects - 3rd order interactions

model_trial_5 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt+ salt_conc + time + 
                         (time|subject) + (1|date_test),data=probos_sugar_reshaped, family=binomial)   
model_trial_5




#Do we need tastant order?
model_trial_6 <- glmer(PER ~ (Background+Allele+Sex)^2 + salt_conc + time + 
                         (time|subject) + (1|date_test),data=probos_sugar_reshaped, family=binomial)   
model_trial_6

AIC(model_trial_2, model_trial_5, model_trial_6)
BIC(model_trial_2, model_trial_5, model_trial_6)

#fit looks slightly better without 3rd order interactions or tastant order.  Confirm with LRparboot.


Anova(model_trial_6, test.statistic="chisq")

####So getting significant effects of allele and background on sugar response.

# Plots the interaction between background and allele on probability scale (really convenient)
# (i.e converts back using inverse logit)
effects_model_MF <- Effect(c("Background", "Allele"), model_trial_6)  
plot(effects_model_MF)
# To do this by hand you need to look at 
effects_model_MF$fit
effects_model_MF$lower
effects_model_MF$upper

# But these need to be converted back to probability scale
InvLog <- function(x) {exp(x)/(1+exp(x))}
InvLog(effects_model_MF$fit) # compare this to a call of effects_model

# also there is a convenience function (author of package really thought ahead) to grab all of this
as.data.frame(effects_model_MF)


effects_model2_MF <- allEffects(model_trial_6)
plot(effects_model2_MF)

####now what about combined NaCl dataset looking at response to salt
probos_salt <- probos[, c(1,2,3,4,5,9,12,14,17)]

probos_salt_reshaped <- reshape(probos_salt, 
                                 varying = list(6:7), 
                                 direction="long", 
                                 idvar="subject", 
                                 v.names="PER")

head(probos_salt_reshaped)
tail(probos_salt_reshaped)

# sorting data by subject    
probos_salt_reshaped<- probos_salt_reshaped[order(probos_salt_reshaped$subject),]

# check that reshaped produced appropriate data set
head(probos_salt_reshaped)    
head(probos)

tail(probos_salt_reshaped)
tail(probos)

#set up for modeling
probos_salt_reshaped$date_test <- as.factor(probos_salt_reshaped$date_test)
probos_salt_reshaped$Allele <- relevel(probos_salt_reshaped$Allele, "wt")

###Fitting glm model with all fitted effects 
probos_salt_reshaped$time <- as.factor(probos_salt_reshaped$time) 

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + salt_conc + time, 
                     data=probos_salt_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
Anova(model_trial_1, test.statistic="LR", type=2)
Anova(model_trial_1, test.statistic="LR", type=3)

# Clear signal, but the date effect is large. try lmer

#checking for date_test effect

model_trial_2 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt + salt_conc + time + (1 + time|subject) + 
                         (1|date_test), data=probos_salt_reshaped, family=binomial)

model_trial_2    

model_trial_3 <- glmer(PER ~ Background*Allele*Sex +  Sugar_before_salt + salt_conc + time + (time|subject), 
                       data=probos_salt_reshaped, family=binomial) 


model_trial_3

AIC(model_trial_1, model_trial_2, model_trial_3)
BIC(model_trial_1, model_trial_2, model_trial_3)

###Again, date_test doesn't look like it needs to be included in the model, but will go ahead with 
###model 2, including date_test.  Glad to see that responses don't seem to be as variable by date
###as we thought they might!

# quick check for the fixed effects (confirm with ParBoot or MCMC)
Anova(model_trial_2, test.statistic="chisq")

Anova(model_trial_3, test.statistic="chisq")

#interesting that neither of these two models suggest that there is an interaction 
#between background and allele.

#Should time by subject random effect be in model?
model_trial_4 <- glmer(PER ~ Background*Allele*Sex + Sugar_before_salt  + salt_conc + time + 
                         (1|subject) + (1|date_test),data=probos_salt_reshaped, family=binomial)   

AIC(model_trial_2, model_trial_3, model_trial_4)
BIC(model_trial_2, model_trial_3, model_trial_4)

#Again - yes!

##now looking at fixed effects - Do we need 3rd order interactions

model_trial_5 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt+ salt_conc + time + 
                         (time|subject) + (1|date_test),data=probos_salt_reshaped, family=binomial)   
model_trial_5

#or tastant order?
model_trial_6 <- glmer(PER ~ (Background+Allele+Sex)^2 + salt_conc + time + 
                         (time|subject) + (1|date_test),data=probos_salt_reshaped, family=binomial)   
model_trial_6

AIC(model_trial_2, model_trial_5, model_trial_6)
BIC(model_trial_2, model_trial_5, model_trial_6)

#fit looks slightly better with 3rd order interactions and tastant order.  Confirm with LRparboot.

###seems like responses for individuals with different alleles change given the salt concentration
###based on interaction plots.  Fitting a background by salt_conc interaction term in the model.

model_trial_7 <- glmer(PER ~ (Background+Allele+Sex)^2 + Background:salt_conc + salt_conc + time + Sugar_before_salt +
                         (time|subject) + (1|date_test),data=probos_salt_reshaped, family=binomial)   
model_trial_7

AIC(model_trial_5, model_trial_7)
BIC(model_trial_5, model_trial_7)

###much better fit when I allow for interactions between salt_conc and Background
###but this should be confirmed with LRparboot
Anova(model_trial_6, test.statistic="chisq")

Anova(model_trial_7, test.statistic="chisq")

####Allele, sex, salt_conc, and tastant order all significant.  Plus a significant interaction
####

# Plots the interaction between background and allele on probability scale (really convenient)
# (i.e converts back using inverse logit)

effects_model_MF <- allEffects(model_trial_6)
plot(effects_model_MF)

effects_model2 <- allEffects(model_trial_7)
plot(effects_model2)

##Now what about for KCl?

##First looking at KCl test response to sugar
KCl100_sugar <- KCl100[, c(1,2,3,4,5,8,11,14)]


# reshaping into long format
KCl100_sugar_reshaped <- reshape(KCl100_sugar, 
                                  varying = list(6:7), 
                                  direction="long", 
                                  idvar="subject", 
                                  v.names="PER")

head(KCl100_sugar_reshaped)
tail(KCl100_sugar_reshaped)

# sorting data by subject    
KCl100_sugar_reshaped<- KCl100_sugar_reshaped[order(KCl100_sugar_reshaped$subject),]

# check that reshaped produced appropriate data set
head(KCl100_sugar_reshaped)    
head(KCl100)

tail(KCl100_sugar_reshaped)
tail(KCl100)

#set up for modeling
KCl100_sugar_reshaped$date_test <- as.factor(KCl100_sugar_reshaped$date_test)
KCl100_sugar_reshaped$Allele <- relevel(KCl100_sugar_reshaped$Allele, "wt")

###Will likely use the full model as in NaCl dataset, but just want to check to see if there are any differences between
###NaCl and KCl 
KCl100_sugar_reshaped$time <- as.factor(KCl100_sugar_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=KCl100_sugar_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3)

#Based upon what we're using for NaCl data:

model_trial_2 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl100_sugar_reshaped, family=binomial)

model_trial_2   

Anova(model_trial_2, test = "chisq")

##Now for salt response during KCl 100 mM tests

##First looking at KCl test response to sugar
KCl100_salt <- KCl100[, c(1,2,3,4,5,9,12,14)]


# reshaping into long format
KCl100_salt_reshaped <- reshape(KCl100_salt, 
                                 varying = list(6:7), 
                                 direction="long", 
                                 idvar="subject", 
                                 v.names="PER")

head(KCl100_salt_reshaped)
tail(KCl100_salt_reshaped)


#set up for modeling
KCl100_salt_reshaped$date_test <- as.factor(KCl100_salt_reshaped$date_test)
KCl100_salt_reshaped$Allele <- relevel(KCl100_salt_reshaped$Allele, "wt")

###Will likely use the full model as in NaCl dataset, but just want to check to see if there are any differences between
###NaCl and KCl 
KCl100_salt_reshaped$time <- as.factor(KCl100_salt_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=KCl100_salt_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3)

#Based upon what we're using for NaCl data:

model_trial_2 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl100_salt_reshaped, family=binomial)

model_trial_2   

Anova(model_trial_2, test = "chisq")

####so significant interaction between background and allele (p = 0.004), confirm with LRparboot


###now looking at KCL 500mM tests

KCl500_sugar <- KCl500[, c(1,2,3,4,5,8,11,14)]


# reshaping into long format
KCl500_sugar_reshaped <- reshape(KCl500_sugar, 
                                 varying = list(6:7), 
                                 direction="long", 
                                 idvar="subject", 
                                 v.names="PER")

head(KCl500_sugar_reshaped)
tail(KCl500_sugar_reshaped)

# sorting data by subject    
KCl500_sugar_reshaped<- KCl500_sugar_reshaped[order(KCl500_sugar_reshaped$subject),]

# check that reshaped produced appropriate data set
head(KCl500_sugar_reshaped)    
head(KCl500)

tail(KCl500_sugar_reshaped)
tail(KCl500)

#set up for modeling
KCl500_sugar_reshaped$date_test <- as.factor(KCl500_sugar_reshaped$date_test)
KCl500_sugar_reshaped$Allele <- relevel(KCl500_sugar_reshaped$Allele, "wt")

###Will likely use the full model as in NaCl dataset, but just want to check to see if there are any differences between
###NaCl and KCl 
KCl500_sugar_reshaped$time <- as.factor(KCl500_sugar_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=KCl500_sugar_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
Anova(model_trial_1, test.statistic="LR", type=2)
Anova(model_trial_1, test.statistic="LR", type=3)

#Based upon what we're using for NaCl data:

model_trial_2 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl500_sugar_reshaped, family=binomial)

model_trial_2   

Anova(model_trial_2, test = "chisq")


##Now for salt response during KCl 100 mM tests

KCl500_salt <- KCl500[, c(1,2,3,4,5,9,12,14)]


# reshaping into long format
KCl500_salt_reshaped <- reshape(KCl500_salt, 
                                varying = list(6:7), 
                                direction="long", 
                                idvar="subject", 
                                v.names="PER")

head(KCl500_salt_reshaped)
tail(KCl500_salt_reshaped)


#set up for modeling
KCl500_salt_reshaped$date_test <- as.factor(KCl500_salt_reshaped$date_test)
KCl500_salt_reshaped$Allele <- relevel(KCl500_salt_reshaped$Allele, "wt")

###Will likely use the full model as in NaCl dataset, but just want to check to see if there are any differences between
###NaCl and KCl 
KCl500_salt_reshaped$time <- as.factor(KCl500_salt_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + Sugar_before_salt + time, 
                     data=KCl500_salt_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3)

#Based upon what we're using for NaCl data:

model_trial_2 <- glmer(PER ~ (Background+Allele+Sex)^2 + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl500_salt_reshaped, family=binomial)

model_trial_2   

Anova(model_trial_2, test = "chisq")

##Now looking at combined dataset

KCl_salt <- KCl_data[, c(1,2,3,4,5,9,12,14,17)]


# reshaping into long format
KCl_salt_reshaped <- reshape(KCl_salt, 
                                varying = list(6:7), 
                                direction="long", 
                                idvar="subject", 
                                v.names="PER")

head(KCl_salt_reshaped)
tail(KCl_salt_reshaped)


#set up for modeling
KCl_salt_reshaped$date_test <- as.factor(KCl_salt_reshaped$date_test)
KCl_salt_reshaped$Allele <- relevel(KCl_salt_reshaped$Allele, "wt")

###Will likely use the full model as in NaCl dataset, but just want to check to see if there are any differences between
###NaCl and KCl 
KCl_salt_reshaped$time <- as.factor(KCl_salt_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Background*Allele*Sex + date_test + salt_conc + Sugar_before_salt + time, 
                     data=KCl_salt_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)
car::Anova(model_trial_1, test.statistic="LR", type=3)

#Based upon what we're using for NaCl data:

model_trial_2 <- glmer(PER ~ (Background+Allele+Sex)^2 + salt_conc + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl_salt_reshaped, family=binomial)

model_trial_2   

Anova(model_trial_2, test = "chisq")

model_trial_3 <- glmer(PER ~ (Background+Allele+Sex)^2 + salt_conc + Background:salt_conc + Sugar_before_salt + time + (1 + time|subject) + 
                         (1|date_test), data=KCl_salt_reshaped, family=binomial)

model_trial_3   

AIC(model_trial_2, model_trial_3)
BIC(model_trial_2, model_trial_3)

Anova(model_trial_3, test = "chisq")