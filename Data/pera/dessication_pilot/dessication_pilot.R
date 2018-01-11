#This is a script to look at the differences in response between flies run 
#with a wet flug vs. just agarose during their 24 h starvation.

setwd("~/Dropbox/Megan/Drosophila_work/PERA_data/dessication_pilot")

x <- c("lattice", "grid", "plyr", "data.table", "car", "lme4")

lapply(x, FUN = function(X) {
  do.call("require", list(X)) 
})

dessic_pilot <- read.table("PERA_dessic_pilot.txt", header = T)

dessic_pilot[21:22, "avg_h2o"] <- NA
dessic_pilot[22:23, "avg_sugar"] <- NA
dessic_pilot[23:24, "avg_salt"] <- NA

dessic_pilot <- transform(dessic_pilot, avg_h2o = rowMeans(dessic_pilot[, c(5, 11, 17)], na.rm = TRUE))
dessic_pilot <- transform(dessic_pilot, avg_sugar = rowMeans(dessic_pilot[, c(7,13)], na.rm = TRUE))
dessic_pilot <- transform(dessic_pilot, avg_salt = rowMeans(dessic_pilot[, c(9,15)], na.rm = TRUE))

dessic_pilot <- data.frame(dessic_pilot)

wetflug <- data.frame(subset(dessic_pilot, Wetted_Flug == "y"))
dryflug <- data.frame(subset(dessic_pilot, Wetted_Flug == "n"))

##looking at average responses by individual to each tastant
plot.new()
tiff("avg_h2o_dens.tif", unit = "px", width = 600, height = 800 )
plot1 <- densityplot(wetflug$avg_h2o, ylim = c(0,5))
plot2 <- densityplot(dryflug$avg_h2o, ylim = c(0,5))
print(plot1, position = c(0.0, 0, 1, 0.5), more = TRUE)
print(plot2, position = c(0, 0.5, 1, 1), more = TRUE)
dev.off()

plot.new()
tiff("avg_sugar_dens.tif", unit = "px", width = 600, height = 800)
plot1 <- densityplot(wetflug$avg_sugar, ylim = c(0,5))
plot2 <- densityplot(dryflug$avg_sugar, ylim = c(0,5))
print(plot1, position = c(0.0, 0, 1, 0.5), more = TRUE)
print(plot2, position = c(0, 0.5, 1, 1), more = TRUE)
dev.off()

plot.new()
tiff("avg_salt_dens.tif", unit = "px", width = 600, height = 800)
plot1 <- densityplot(wetflug$avg_salt, ylim = c(0,5))
plot2 <- densityplot(dryflug$avg_salt, ylim = c(0,5))
print(plot1, position = c(0.0, 0, 1, 0.5), more = TRUE)
print(plot2, position = c(0, 0.5, 1, 1), more = TRUE)
dev.off()

# May be easier to see like this
tiff("dessic_avg_h2o.tif")
plot.new()
par(mfrow = c(2,1))
stripchart(jitter(wetflug$avg_h2o) ~ wetflug$Geno) 
stripchart(jitter(dryflug$avg_h2o) ~ dryflug$Geno)
dev.off()

tiff("dessic_avg_sugar.tif")
plot.new()
par(mfrow = c(2, 1))
stripchart(jitter(wetflug$avg_sugar) ~ wetflug$Geno) 
stripchart(jitter(dryflug$avg_sugar) ~ dryflug$Geno)
dev.off()

tiff("dessic_avg_salt.tif")
plot.new()
par(mfrow = c(2, 1))
stripchart(jitter(wetflug$avg_salt) ~ wetflug$Geno) 
stripchart(jitter(dryflug$avg_salt) ~ dryflug$Geno)
dev.off()

#reformatting to fit a model for water response

dessic_pilot_h2o <- dessic_pilot[, c(1,2,3,4,5,11,17,19,20)]


# reshaping into long format
dp_reshaped <- reshape(dessic_pilot_h2o, 
                                  varying = list(5:7), 
                                  direction="long", 
                                  idvar="subject", 
                                  v.names="PER")

# sorting data by subject    
dp_reshaped<- dp_reshaped[order(dp_reshaped$subject),]

# check that reshaped produced appropriate data set
head(dp_reshaped)    
head(dessic_pilot_h2o)
tail(dp_reshaped)
tail(dessic_pilot_h2o)


dp_reshaped$time <- as.factor(dp_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ (Geno+Sex+Wetted_Flug)^2 + time, 
                     data=dp_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)


#reformatting to fit a model for sugar

dessic_pilot_sug <- dessic_pilot[, c(1,2,3,4,7,13,19,20)]


# reshaping into long format
dpsug_reshaped <- reshape(dessic_pilot_sug, 
                       varying = list(5:6), 
                       direction="long", 
                       idvar="subject", 
                       v.names="PER")

# sorting data by subject    
dpsug_reshaped<- dpsug_reshaped[order(dpsug_reshaped$subject),]

# check that reshaped produced appropriate data set
head(dpsug_reshaped)    
head(dessic_pilot_sug)
tail(dpsug_reshaped)
tail(dessic_pilot_sug)

dpsug_reshaped$time <- as.factor(dpsug_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Geno*Sex+Wetted_Flug + time, 
                     data=dpsug_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)


#now reformatting for a quick look at salt responses

dessic_pilot_kcl <- dessic_pilot[, c(1,2,3,4,9,15,19,20)]


# reshaping into long format
dpkcl_reshaped <- reshape(dessic_pilot_kcl, 
                          varying = list(5:6), 
                          direction="long", 
                          idvar="subject", 
                          v.names="PER")

# sorting data by subject    
dpkcl_reshaped<- dpkcl_reshaped[order(dpkcl_reshaped$subject),]

# check that reshaped produced appropriate data set
head(dpkcl_reshaped)    
head(dessic_pilot_kcl)
tail(dpkcl_reshaped)
tail(dessic_pilot_kcl)

dpkcl_reshaped$time <- as.factor(dpkcl_reshaped$time) # only 2 "levels": first and second, so treat as factor.

model_trial_1 <- glm(PER ~ Geno*Sex+Wetted_Flug + time, 
                     data=dpkcl_reshaped, family=binomial, x=T)

summary(model_trial_1)
anova(model_trial_1, test="LRT")
car::Anova(model_trial_1, test.statistic="LR", type=2)

exp(cbind(OR = coef(model_trial_1), confint(model_trial_1)))

