PERA\_Exploratory\_Data\_Analysis
================
Megan Fritz
December 19, 2017

Loading in data sets and libraries.
===================================

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

Adding barometric pressure data
===============================

Because changes in barometric pressure are known to modulate insect behavioral responses (see Ankney 1984, Pellegrino et al. 2013), we explored the potential it had to influence gustation responses in our assays. These data were obtained from the weather underground historical data archive (<http://www.wunderground.com/history>) for Lansing, MI.

``` r
NaCl_barom <- read.table("./Data/pera/NaCl_barom_readings.txt", header = T)
NaCl_barom$delta_pres <- (NaCl_barom$pressure_6am - NaCl_barom$pressure_9am)

KCl_barom <- read.table("./Data/pera/KCl_barom_readings.txt", header = T)
KCl_barom$delta_pres <- (KCl_barom$pressure_6am - KCl_barom$pressure_9am)

#merging barometric pressure into data sets.
NaCl_data_pres <- merge(NaCl_data, NaCl_barom, by="date_test")
KCl_data_pres <- merge(KCl_data, KCl_barom, by="date_test")
```

Calculating average per fly responses
=====================================

Our goal was to examine whether gustation responses were influenced by genetic background, scalloped allele, sex, etc. We first visualized these relationships by plotting average fly responses to each tastant (water, sugar, salt) by each of these factors. We first had to generate these average responses per fly for each of two experiments targeting the labellar sensilla - one where the aversive salt tastant was NaCl and another where the salt was KCl.

``` r
#Getting average response per fly to each stimulus
#NaCl
NaCl_data_pres[24:25, "avg_h2o"] <- NA
NaCl_data_pres[25:26, "avg_sugar"] <- NA
NaCl_data_pres[26:27, "avg_salt"] <- NA

NaCl_data_pres <- transform(NaCl_data_pres, avg_h2o = rowMeans(NaCl_data_pres[, c(7, 10, 13)], na.rm = TRUE))
NaCl_data_pres <- transform(NaCl_data_pres, avg_sugar = rowMeans(NaCl_data_pres[, c(8,11)], na.rm = TRUE))
NaCl_data_pres <- transform(NaCl_data_pres, avg_salt = rowMeans(NaCl_data_pres[, c(9,12)], na.rm = TRUE))

#KCl
KCl_data_pres[24:25, "avg_h2o"] <- NA
KCl_data_pres[25:26, "avg_sugar"] <- NA
KCl_data_pres[26:27, "avg_salt"] <- NA

KCl_data_pres <- transform(KCl_data_pres, avg_h2o = rowMeans(KCl_data_pres[, c(7, 10, 13)], na.rm = TRUE))
KCl_data_pres <- transform(KCl_data_pres, avg_sugar = rowMeans(KCl_data_pres[, c(8,11)], na.rm = TRUE))
KCl_data_pres <- transform(KCl_data_pres, avg_salt = rowMeans(KCl_data_pres[, c(9,12)], na.rm = TRUE))


###First looking at avg responses by individuals tested using labellar receptors
probos_N <- subset(NaCl_data_pres, receptor == "prob")  
probos_K <- subset(KCl_data_pres, receptor == "prob")
```

The bootstrap function used to generate CIs around average responses to each tastant (water, sugar, and salt) in the following plots.

``` r
#writing bootstrapped 95% CI function
boot.fn <- function(x, N=5000) {
  Int.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  Int.CI <- quantile(Int.1, probs=c(0.025,0.975))
  Int.CI
}
```

Plotting average per fly responses
==================================

Then I applied the mean and bootstrap functions to my dataset to calculate and plot average per fly response to each tastant (water, sugar, then salt) when applied to the labellar sensilla. Note that the plots allow us to visualize whether there is an interaction between allele and genetic background. Responses to each tastant are plotted separately for NaCl and KCl experiments.

Mean water response rates (2.5, 97.5% CIs)
==========================================

``` r
#printing overall means and CIs
mean.avg_h2o.NaCl <- tapply(probos_N$avg_h2o, probos_N$Allele, mean)
CI.avg_h2o.NaCl <- tapply(probos_N$avg_h2o, probos_N$Allele, boot.fn)
mean.avg_h2o.NaCl
```

    ##        58d       etx4        sd1       sde3         wt 
    ## 0.11666667 0.09251101 0.16666667 0.11484099 0.06666667

``` r
CI.avg_h2o.NaCl
```

    ## $`58d`
    ##       2.5%      97.5% 
    ## 0.06666667 0.17222222 
    ## 
    ## $etx4
    ##       2.5%      97.5% 
    ## 0.07048458 0.11747430 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.1329966 0.2037037 
    ## 
    ## $sde3
    ##       2.5%      97.5% 
    ## 0.08892815 0.14369847 
    ## 
    ## $wt
    ##       2.5%      97.5% 
    ## 0.05305556 0.08056250

``` r
#Plot to view interactions between background and allele for water responses in NaCl assay
lineplot.CI(Allele, avg_h2o, group = Background, data = probos_N, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_h2o_resp", ylim = c(0, 0.3), cex.lab = 1.5, 
                     col = c("blue",  "red"), 
                     pch = c(16,16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_water_responses-1.png)

``` r
#printing overall means and CIs
mean.avg_h2o.KCl <- tapply(probos_K$avg_h2o, probos_K$Allele, mean)
CI.avg_h2o.KCl <- tapply(probos_K$avg_h2o, probos_K$Allele, boot.fn)
mean.avg_h2o.KCl
```

    ##        58d       etx4        sd1       sde3         wt 
    ## 0.03875969 0.13787879 0.13356974 0.15248227 0.05366922

``` r
CI.avg_h2o.KCl
```

    ## $`58d`
    ##       2.5%      97.5% 
    ## 0.01937984 0.06201550 
    ## 
    ## $etx4
    ##      2.5%     97.5% 
    ## 0.1060606 0.1696970 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.1039894 0.1654846 
    ## 
    ## $sde3
    ##      2.5%     97.5% 
    ## 0.1205674 0.1855792 
    ## 
    ## $wt
    ##       2.5%      97.5% 
    ## 0.04381161 0.06389193

``` r
#Plot to view interactions between background and allele for water responses in KCl assay
lineplot.CI(Allele, avg_h2o, group = Background, data = probos_K, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_h2o_resp", ylim = c(0, 0.3), cex.lab = 1.5, 
                     col = c("blue", "red"), pch = c(16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_water_responses-2.png)

Mean sugar response rates (2.5, 97.5% CIs)
==========================================

``` r
#printing overall means and CIs
mean.avg_sug.NaCl <- tapply(probos_N$avg_sugar, probos_N$Allele, mean)
CI.avg_sug.NaCl <- tapply(probos_N$avg_sugar, probos_N$Allele, boot.fn)
mean.avg_sug.NaCl
```

    ##       58d      etx4       sd1      sde3        wt 
    ## 0.3500000 0.3920705 0.5984848 0.5636042 0.3241667

``` r
CI.avg_sug.NaCl
```

    ## $`58d`
    ##      2.5%     97.5% 
    ## 0.2583333 0.4500000 
    ## 
    ## $etx4
    ##      2.5%     97.5% 
    ## 0.3370044 0.4471366 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.5429293 0.6540404 
    ## 
    ## $sde3
    ##      2.5%     97.5% 
    ## 0.5123675 0.6148410 
    ## 
    ## $wt
    ##      2.5%     97.5% 
    ## 0.2925000 0.3566667

``` r
#Plot to view interactions between background and allele for sugar responses in NaCl assay
lineplot.CI(Allele, avg_sugar, group = Background, data = probos_N, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_sugar_resp", ylim = c(0, 0.8), cex.lab = 1.5, 
                     col = c("blue",  "red"), 
                     pch = c(16,16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_sugar_responses-1.png)

``` r
#printing overall means and CIs
mean.avg_sug.KCl <- tapply(probos_K$avg_sugar, probos_K$Allele, mean)
CI.avg_sug.KCl <- tapply(probos_K$avg_sugar, probos_K$Allele, boot.fn)
mean.avg_sug.KCl
```

    ##       58d      etx4       sd1      sde3        wt 
    ## 0.3546512 0.6272727 0.5195035 0.5886525 0.3203724

``` r
CI.avg_sug.KCl
```

    ## $`58d`
    ##      2.5%     97.5% 
    ## 0.2732558 0.4418605 
    ## 
    ## $etx4
    ##      2.5%     97.5% 
    ## 0.5727273 0.6795455 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.4698582 0.5691489 
    ## 
    ## $sde3
    ##      2.5%     97.5% 
    ## 0.5407801 0.6365248 
    ## 
    ## $wt
    ##      2.5%     97.5% 
    ## 0.2940854 0.3472070

``` r
#Plot to view interactions between background and allele for sugar responses in KCl assay
lineplot.CI(Allele, avg_sugar, group = Background, data = probos_K, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_sugar_resp", ylim = c(0, 0.8), cex.lab = 1.5, 
                     col = c("blue", "red"), pch = c(16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_sugar_responses-2.png)

Mean salt response rates (2.5, 97.5% CIs)
=========================================

``` r
#printing overall means and CIs
mean.avg_salt.NaCl <- tapply(probos_N$avg_salt, probos_N$Allele, mean)
CI.avg_salt.NaCl <- tapply(probos_N$avg_salt, probos_N$Allele, boot.fn)
mean.avg_salt.NaCl
```

    ##       58d      etx4       sd1      sde3        wt 
    ## 0.1583333 0.2665198 0.4015152 0.3780919 0.1308333

``` r
CI.avg_salt.NaCl
```

    ## $`58d`
    ##  2.5% 97.5% 
    ## 0.075 0.250 
    ## 
    ## $etx4
    ##      2.5%     97.5% 
    ## 0.2202643 0.3149780 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.3434343 0.4621212 
    ## 
    ## $sde3
    ##      2.5%     97.5% 
    ## 0.3286219 0.4293286 
    ## 
    ## $wt
    ##      2.5%     97.5% 
    ## 0.1083333 0.1541667

``` r
#Plot to view interactions between background and allele for salt responses in NaCl assay
lineplot.CI(Allele, avg_salt, group = Background, data = probos_N, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_salt_resp", ylim = c(0, 0.8), cex.lab = 1.5, 
                     col = c("blue",  "red"), 
                     pch = c(16,16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_salt_responses-1.png)

``` r
#printing overall means and CIs
mean.avg_salt.KCl <- tapply(probos_K$avg_salt, probos_K$Allele, mean)
CI.avg_salt.KCl <- tapply(probos_K$avg_salt, probos_K$Allele, boot.fn)
mean.avg_salt.KCl
```

    ##        58d       etx4        sd1       sde3         wt 
    ## 0.18604651 0.31818182 0.15425532 0.33687943 0.05805038

``` r
CI.avg_salt.KCl
```

    ## $`58d`
    ##     2.5%    97.5% 
    ## 0.127907 0.250000 
    ## 
    ## $etx4
    ##      2.5%     97.5% 
    ## 0.2636364 0.3750000 
    ## 
    ## $sd1
    ##      2.5%     97.5% 
    ## 0.1187943 0.1914894 
    ## 
    ## $sde3
    ##      2.5%     97.5% 
    ## 0.2854610 0.3865248 
    ## 
    ## $wt
    ##       2.5%      97.5% 
    ## 0.04600219 0.07064622

``` r
#Plot to view interactions between background and allele for salt responses in KCl assay
lineplot.CI(Allele, avg_salt, group = Background, data = probos_K, cex = 1.5, xlab = "Allele", 
                     ylab = "avg_salt_resp", ylim = c(0, 0.8), cex.lab = 1.5, 
                     col = c("blue", "red"), pch = c(16,16,16,16),
                     ci.fun= boot.fn)
```

![](PERA_Exploratory_Data_Analysis_MF_files/figure-markdown_github/avg_salt_responses-2.png)

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
