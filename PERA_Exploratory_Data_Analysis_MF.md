PERA\_Exploratory\_Data\_Analysis
================
Megan Fritz
December 19, 2017

Loading in data sets
====================

``` r
#use this to specify the path to your data file

setwd("~/Dropbox/Megan/Drosophila_work") 

#Loading data sets
NaCl_data <- read.table("./Data/pera/PERA_NaCl_edited.txt", header = T)
NaCl_data <- data.frame(NaCl_data)

KCl_data <- read.table("./Data/pera/PERA_KCl.txt", header = T)
KCl_data <- data.frame(KCl_data)

#loading libraries
library(sciplot)
```

Because changes in barometric pressure are known to modulate insect behavioral responses (see Ankney 1984, Pellegrino et al. 2013), we explored the potential it had to influence gustation responses in our assays.

``` r
NaCl_barom <- read.table("./Data/pera/NaCl_barom_readings.txt", header = T)
NaCl_barom$delta_pres <- (NaCl_barom$pressure_6am - NaCl_barom$pressure_9am)

KCl_barom <- read.table("./Data/pera/KCl_barom_readings.txt", header = T)
KCl_barom$delta_pres <- (KCl_barom$pressure_6am - KCl_barom$pressure_9am)

#merging barometric pressure into data sets.
NaCl_data_pres <- merge(NaCl_data, NaCl_barom, by="date_test")
KCl_data_pres <- merge(KCl_data, KCl_barom, by="date_test")
```

Our goal was to examine whether gustation responses were influenced by genetic background, scalloped allele, sex, etc. We first visualized these relationships by plotting average fly responses to each tastant (water, sugar, salt) by each of these factors.

``` r
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
probos_N <- subset(NaCl_data_pres, receptor == "prob")  
probos_K <- subset(KCl_data_pres, receptor == "prob")
```

Average Per Fly Responses to water applied to the labellar sensilla.

``` r
tiff("./Output/mean_h2o_bybackground&geno2.tif", unit = "px", height = 800, width = 600)
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
```

    ## png 
    ##   2

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.