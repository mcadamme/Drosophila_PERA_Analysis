PERA\_Full\_Model\_Selection
================
Megan Fritz
February 20, 2018

Preparing data sets and loading libraries.
------------------------------------------

``` r
#use this to specify the path to your data file

setwd("~/Dropbox/Megan/Drosophila_work") 

#Setting up my data for statistical modelling
NaCl_data <- read.table("./Data/pera/PERA_NaCl_edited.txt", header = T)
NaCl_data <- data.frame(NaCl_data)

KCl_data <- read.table("./Data/pera/PERA_KCl.txt", header = T)
KCl_data <- data.frame(KCl_data)

#QC - checking for duplicate rows in each data file
sum(duplicated(NaCl_data))
```

    ## [1] 0

``` r
sum(duplicated(KCl_data))
```

    ## [1] 0

``` r
#adding in barometric pressure data - see PERA_Exploratory_Data_Analysis script for details
NaCl_barom <- read.table("./Data/pera/NaCl_barom_readings.txt", header = T)
NaCl_barom$delta_pres <- (NaCl_barom$pressure_6am - NaCl_barom$pressure_9am)

KCl_barom <- read.table("./Data/pera/KCl_barom_readings.txt", header = T)
KCl_barom$delta_pres <- (KCl_barom$pressure_6am - KCl_barom$pressure_9am)

#merging barometric pressure into data sets.
NaCl_data_pres <- merge(NaCl_data, NaCl_barom, by="date_test")
KCl_data_pres <- merge(KCl_data, KCl_barom, by="date_test")

#Looking at each sensory structure separately for each tastant
probos_N <- subset(NaCl_data_pres, receptor == "prob")
tars_N <- subset(NaCl_data_pres, receptor == "tars")
probos_K <- subset(KCl_data_pres, receptor == "prob")

#generating the subject levels for use as a random effect
probos_N$subject <- with(probos_N, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))
tars_N$subject <- with(tars_N, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))
probos_K$subject <- with(probos_K, interaction(Allele, Background, Sex, date_test, Fly, salt_conc, drop=T, sep="_"))
```

``` r
#loading relevant libraries
x <- c("sciplot", "lme4", "car", "effects", "bbmle", "ggplot2", "gridExtra", "RCurl", "optimx", "MCMCglmm") 

lapply(x, FUN = function(X) {
  do.call("library", list(X)) 
})
```
