###############################################################################
###############################################################################
# Replication code for figures and tables presented in:

# Authors: Buntaine,Mark T.; Hamilton,Stuart E.; Millones,Marco
# Title: Titling community land to prevent deforestation: An evaluation of 
#a best-case program in Morona-Santiago,Ecuador
# Journal: Global Environmental Change,vol. 33,p. 32-43.

# Replication code for statistical analysis and 
# matching prepared by Mark Buntaine
# Mark Buntaine contact (as of May 2015): buntaine@bren.ucsb.edu / 805-893-4075

# Current for R Version 3.1.3 (version "Smooth Sidewalk")
###############################################################################
###############################################################################

# TDR: Note that there are a number of comments in the text from the authors 
# themselves. To differentiate my comments, I will put my initials at the start
# of the comment block, as I have done here.

###############################################################################

# Setting up
# TDR: Loading the packages. arm is a modeling package, MatchIt is a package 
# (by Kosuke!) for matching treated-untreated cases for causal inference (we
# will see below that they use genetic matching), coefplot is for plotting 
# coefficients, ggplot2 we know well, and rgenoud is for genetic optimization;
# not sure how it's used here yet. 

require(arm) #Version 1.7-07
require(MatchIt) #Version 2.4-21
require(coefplot) #Version 1.2.0
require(ggplot2) #Version 1.0.1
require(rgenoud)

# TDR: I've changed this to work from the repo's top level directory

data <- read.csv("TitlingDeforestation/BHM_sample20.csv")


# TDR: these notes are from the author
# Note: the current sample20 is modified to move first two columns 
# from raw GIS export
# Note: A352 moved after A460 to match earlier versions and keep order of cols
# Note: all NULL changed to NA with find and replace in Excel

## Treatment variables 
# TDR: this creates two new variables to indicate treatment. Var A352 indicates
# intervention year. 0 is no intervention, 2 is 2002, 3 is 2003, etc. Therefore
# PML indicates whether or not treatment happened at all. PMI indicates whether 
# or not the intervention occurred in 2003 or later. At least, according to the
# metadata. In the supporting information, PML indicates parcels that were only
# given legal tenure, while PMI indicates ones given legal tenure and a 
# complementary community management plan. Need to clarify. 

data$PML <- ifelse(data$A352 == 0,0,1)
data$PMI <- ifelse(data$A352 >= 3,1,0) 

# Note: all PSUR intervention had PMI in 2003 and after
# TDR: Just renaming cols, and then converting the "distance to Peruvian 
# border" variable from meters to km. The renamed cols are ones corresponding
# to yearly forest loss data

names(data)[287:298] <- c("fl.2001","fl.2002","fl.2003","fl.2004",
                          "fl.2005","fl.2006", "fl.2007","fl.2008",
                          "fl.2009","fl.2010","fl.2011","fl.2012")
data$toPeru.km <- data$A486 / 1000

# This function remakes the data frame with only specified columns & NAs removed:
# TDR: this defines a function, make(), which reworks the raw dataframe into a 
# separate, model-specific one. 

make <- function (data, names) {
  datasub <- as.data.frame(data)
  sub = substitute(names)
  names <- as.character(sub)
  data.make <- subset(datasub,select = names(datasub) %in% names)
  data.make <<- na.omit(data.make)
  data.make <- na.omit(data.make)
  print(c("ob =",nrow(data.make)),quote=F)
}

# Example> make(c(monit2other,sa1all))

###############################################################################

# Creating "flat" dataset that allows matching across pre-treatment years
# TDR: for roughly the next 100 lines they go over each year in the data (each
# year is 17 lines) creating a subset and then creating new, year-specific
# variables. I comment line-by-line for the first year:

# TDR: first, subset/select for observations that either were not treated or
# were treated in the relevant year (eg 2002)
wave.2002 <- subset(data, (A352 == 2 | A352 == 0))
# TDR: new var total deforestation for five years after treatment 
wave.2002$fl.5yr.postPML <- wave.2002$fl.2003 + wave.2002$fl.2004 + 
                            wave.2002$fl.2005 + wave.2002$fl.2006 + 
                            wave.2002$fl.2007
# TDR: these rename annual forest loss to a standard form
wave.2002$fl.yr1 <- wave.2002$fl.2003
wave.2002$fl.yr2 <- wave.2002$fl.2004
wave.2002$fl.yr3 <- wave.2002$fl.2005
wave.2002$fl.yr4 <- wave.2002$fl.2006
wave.2002$fl.yr5 <- wave.2002$fl.2007
# TDR: this var captures pre-treatment forest loss
wave.2002$fl.pre.5cell <- wave.2002$fl.2001 + wave.2002$A393 + 
                          wave.2002$A461
# TDR: this var is the diff of A408 (area of pixel with forest cover in 2000)
# and fl.2001, forest loss 
wave.2002$forest.pre <- wave.2002[, 314] - wave.2002[, 287]
# TDR: this is "donut population density" from LandScan, a population density
# dataset made by Oak Ridge Laboratory
wave.2002$pop.den.1km.pre <- wave.2002$A173
# TDR: same thing, different res
wave.2002$pop.den.5km.pre <- wave.2002$A176
# TDR: ditto
wave.2002$pop.den.10km.pre <- wave.2002$A174
# TDR: this var is based on MODIS data (MODIS is a specific sensor carried
# by a pair of satellites, more information on request). One of the standard
# products distributed by MODIS is MCD12Q, a land use classification scheme. 
# They categorize some uses as "disturbed," such as agriculture or built 
# environment. This var is the distance of a given pixel from a pixel that
# was disturbed in 2002. 
wave.2002$dist.disturb.pre <- wave.2002$A116
# TDR: not sure yet why they need a column of the word two 
wave.2002$wave = "two"
# TDR: subsetting the df for pixels with more than 450000m2 of forest and 
# which had formal land rights by 2008
wave.2002 <- subset(wave.2002, forest.pre > 450000 & A346 == 0)

wave.2003 <- subset(data,(A352 == 3 | A352 == 0))
wave.2003$fl.5yr.postPML <- wave.2003$fl.2004 + wave.2003$fl.2005 + 
                            wave.2003$fl.2006 + wave.2003$fl.2007 + 
                            wave.2003$fl.2008
wave.2003$fl.yr1 <- wave.2003$fl.2004
wave.2003$fl.yr2 <- wave.2003$fl.2005
wave.2003$fl.yr3 <- wave.2003$fl.2006
wave.2003$fl.yr4 <- wave.2003$fl.2007
wave.2003$fl.yr5 <- wave.2003$fl.2008
wave.2003$fl.pre.5cell <- wave.2003$fl.2002 + wave.2003$A394 + wave.2003$A462
wave.2003$forest.pre <- wave.2003[,314] - rowSums(wave.2003[,287:288])
wave.2003$pop.den.1km.pre <- wave.2003$A178
wave.2003$pop.den.5km.pre <- wave.2003$A181
wave.2003$pop.den.10km.pre <- wave.2003$A179
wave.2003$dist.disturb.pre <- wave.2003$A117
wave.2003$wave = "three"
wave.2003 <- subset(wave.2003, forest.pre > 450000 & A346 == 0)

wave.2004 <- subset(data, A352 == 0)
wave.2004$fl.5yr.postPML <- wave.2004$fl.2005 + wave.2004$fl.2006 +
                            wave.2004$fl.2007 + wave.2004$fl.2008 + 
                            wave.2004$fl.2009
wave.2004$fl.yr1 <- wave.2004$fl.2005
wave.2004$fl.yr2 <- wave.2004$fl.2006
wave.2004$fl.yr3 <- wave.2004$fl.2007
wave.2004$fl.yr4 <- wave.2004$fl.2008
wave.2004$fl.yr5 <- wave.2004$fl.2009
wave.2004$fl.pre.5cell <- wave.2004$fl.2003 + wave.2004$A395 + wave.2004$A463
wave.2004$forest.pre <- wave.2004[, 314] - rowSums(wave.2004[, 287:289])
wave.2004$pop.den.1km.pre <- wave.2004$A183
wave.2004$pop.den.5km.pre <- wave.2004$A186
wave.2004$pop.den.10km.pre <- wave.2004$A184
wave.2004$dist.disturb.pre <- wave.2004$A118
wave.2004$wave = "four"
wave.2004 <- subset(wave.2004, forest.pre > 450000 & A346 == 0)

wave.2005 <- subset(data, (A352 == 5 | A352 == 0))
wave.2005$fl.5yr.postPML <- wave.2005$fl.2006 + wave.2005$fl.2007 + 
                            wave.2005$fl.2008 + wave.2005$fl.2009 +
                            wave.2005$fl.2010
wave.2005$fl.yr1 <- wave.2005$fl.2006
wave.2005$fl.yr2 <- wave.2005$fl.2007
wave.2005$fl.yr3 <- wave.2005$fl.2008
wave.2005$fl.yr4 <- wave.2005$fl.2009
wave.2005$fl.yr5 <- wave.2005$fl.2010
wave.2005$fl.pre.5cell <- wave.2005$fl.2004 + wave.2005$A396 + wave.2005$A464
wave.2005$forest.pre <- wave.2005[,314] - rowSums(wave.2005[,287:290])
wave.2005$pop.den.1km.pre <- wave.2005$A188
wave.2005$pop.den.5km.pre <- wave.2005$A191
wave.2005$pop.den.10km.pre <- wave.2005$A189
wave.2005$dist.disturb.pre <- wave.2005$A119
wave.2005$wave = "five"
wave.2005 <- subset(wave.2005, forest.pre > 450000 & A346 == 0)

wave.2006 <- subset(data,(A352 == 6 | A352 == 0))
wave.2006$fl.5yr.postPML <- wave.2006$fl.2007 + wave.2006$fl.2008 + 
                            wave.2006$fl.2009 + wave.2006$fl.2010 + 
                            wave.2006$fl.2011
wave.2006$fl.yr1 <- wave.2006$fl.2007
wave.2006$fl.yr2 <- wave.2006$fl.2008
wave.2006$fl.yr3 <- wave.2006$fl.2009
wave.2006$fl.yr4 <- wave.2006$fl.2010
wave.2006$fl.yr5 <- wave.2006$fl.2011
wave.2006$fl.pre.5cell <- wave.2006$fl.2005 + wave.2006$A397 + wave.2006$A465
wave.2006$forest.pre <- wave.2006[, 314] - rowSums(wave.2006[, 287:291])
wave.2006$pop.den.1km.pre <- wave.2006$A193
wave.2006$pop.den.5km.pre <- wave.2006$A196
wave.2006$pop.den.10km.pre <- wave.2006$A194
wave.2006$dist.disturb.pre <- wave.2006$A120
wave.2006$wave = "six"
wave.2006 <- subset(wave.2006, forest.pre > 450000 & A346 == 0)

wave.2007 <- subset(data,(A352 == 7 | A352 == 0))
wave.2007$fl.5yr.postPML <- wave.2007$fl.2008 + wave.2007$fl.2009 + 
                            wave.2007$fl.2010 + wave.2007$fl.2011 + 
                            wave.2007$fl.2012
wave.2007$fl.yr1 <- wave.2007$fl.2008
wave.2007$fl.yr2 <- wave.2007$fl.2009
wave.2007$fl.yr3 <- wave.2007$fl.2010
wave.2007$fl.yr4 <- wave.2007$fl.2011
wave.2007$fl.yr5 <- wave.2007$fl.2012
wave.2007$fl.pre.5cell <- wave.2007$fl.2006 + wave.2007$A398 + wave.2007$A466
wave.2007$forest.pre <- wave.2007[, 314] - rowSums(wave.2007[, 287:292])
wave.2007$pop.den.1km.pre <- wave.2007$A198
wave.2007$pop.den.5km.pre <- wave.2007$A201
wave.2007$pop.den.10km.pre <- wave.2007$A199
wave.2007$dist.disturb.pre <- wave.2007$A121
wave.2007$wave = "seven"
wave.2007 <- subset(wave.2007, forest.pre > 450000 & A346 == 0)

data.flat <- rbind(wave.2002, wave.2003, wave.2004, wave.2005, wave.2006,
                   wave.2007)
data.flat$id <- paste(data.flat$wave, data.flat$OBJECTID, sep="")

###############################################################################
# Figure 4

# Excluding PMI plots from control group

data.PMLonly.flat <- subset(data.flat, PML == 0 | (PML == 1 & PMI == 0)) 

# Models for PML

fig4.a <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + forest.pre +
  pop.den.5km.pre + dist.disturb.pre + A128 + A129 + 
  A130 + A154 + A347 + A361 + A362,data=data.PMLonly.flat)
summary(fig4.a)

make(data.PMLonly.flat, c(fl.5yr.postPML,fl.yr1,fl.yr2,fl.yr3,fl.yr4, 
  fl.yr5,PML,fl.pre.5cell,forest.pre,
  pop.den.5km.pre,dist.disturb.pre,A128,A129,
  A130,A154,A347,A361,A362,id))
out.PMLonly.flat <- matchit(PML ~ fl.pre.5cell + forest.pre + 
  pop.den.5km.pre + dist.disturb.pre +
  A128 + A129 + A130 + A154 + A347 +
  A361 + A362, method = "genetic",
  int.seed = 104, unif.seed = 104,
  pop.size = 1500, wait.generations = 20,
  MemoryMatrix = FALSE, ties = FALSE,
  discard = "control", data = data.make)

# Note: matching takes approximately 1 week on standard-build desktop 
# computer at time of writing

# Information used for SI Table B1

summary(out.PMLonly.flat, interactions = T, standardize = T) 

data.PMLonly.flat.matched <- match.data(out.PMLonly.flat)

# Note: data.PMLonly.flat.matched saved as "PMLonly_flat_matched_140414.csv"
# and available from authors on request

fig4.b <- lm(fl.5yr.postPML ~ PML, weights = weights, 
             data = data.PMLonly.flat.matched)
summary(fig4.b)

fig4.c <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + f
             orest.pre + pop.den.5km.pre + dist.disturb.pre + 
             A128 + A129 + A130 + A154 + A347 + A361 + A362, 
             weights = weights, data = data.PMLonly.flat.matched)
summary(fig4.c)

fig4 <- multiplot(fig4.a, fig4.b, fig4.c, coefficients = "PML", 
                  names = c("(a)", "(b)", "(c)"), title=NULL,
                  ylab=NULL, xlab = expression(paste(
                    "Effect on deforestation over 5 years (",
                    m^{2} , " per " , km^{2}, ")")), pointSize = 5, 
                    lwdOuter = 2, horizontal = TRUE, secret.weapon = TRUE)
fig4 + theme_bw(base_size=18) 

###############################################################################

# Figure 5

# Model 1 is year-on-year decomposition of Model (a) in Figure 4

m1.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362, data = data.PMLonly.flat)
summary(m1.PML1)

m1.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362,data = data.PMLonly.flat)
summary(m1.PML2)

m1.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362,data = data.PMLonly.flat)
summary(m1.PML3)

m1.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362,data = data.PMLonly.flat)
summary(m1.PML4)

m1.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362,data = data.PMLonly.flat)
summary(m1.PML5)

# Model 2 is year-on-year decomposition of Model (b) in Figure 4

m2.PML1 <- lm(fl.yr1 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML1)

m2.PML2 <- lm(fl.yr2 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML2)

m2.PML3 <- lm(fl.yr3 ~ PML, weights = weights,
              data = data.PMLonly.flat.matched)
summary(m2.PML3)

m2.PML4 <- lm(fl.yr4 ~ PML, weights = weights,
              data = data.PMLonly.flat.matched)
summary(m2.PML4)

m2.PML5 <- lm(fl.yr5 ~ PML,weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML5)

# Model 3 is year-on-year decomposition of Model (c) in Figure 6

m3.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362, weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML1)

m3.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362, weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML2)

m3.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362, weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML3)

m3.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362, weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML4)

m3.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362, weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML5)

##Year by year effects of PML

year <- seq(1:5)
m1.coef <- c(m1.PML1$coef[2], m1.PML2$coef[2], m1.PML3$coef[2], 
             m1.PML4$coef[2], m1.PML5$coef[2])
m1.se <- c(coef(summary(m1.PML1))[2, 2], coef(summary(m1.PML2))[2, 2], 
           coef(summary(m1.PML3))[2, 2], coef(summary(m1.PML4))[2, 2], 
           coef(summary(m1.PML5))[2, 2])

m2.coef <- c(m2.PML1$coef[2], m2.PML2$coef[2], m2.PML3$coef[2], 
             m2.PML4$coef[2], m2.PML5$coef[2])
m2.se <- c(coef(summary(m2.PML1))[2, 2], coef(summary(m2.PML2))[2, 2], 
           coef(summary(m2.PML3))[2, 2], coef(summary(m2.PML4))[2, 2], 
           coef(summary(m2.PML5))[2, 2])

m3.coef <- c(m3.PML1$coef[2], m3.PML2$coef[2], m3.PML3$coef[2],
             m3.PML4$coef[2], m3.PML5$coef[2])
m3.se <- c(coef(summary(m3.PML1))[2, 2], coef(summary(m3.PML2))[2, 2], 
           coef(summary(m3.PML3))[2, 2], coef(summary(m3.PML4))[2, 2], 
           coef(summary(m3.PML5))[2, 2])

par(mar = c(4.1,5.1,1.1,1.1))
plot(1, type = "n", xlim = c(0.8, 5.2), ylim = c(-600, 400), 
     ylab = expression(paste("Effect on Deforestation (",
                       m^{2} , " per " , km^{2}, ")")), 
     xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
abline(h = 0,lty = 2)
polygon(x=c(year[1:5], year[5:1]), y = c(m1.coef + 2 * m1.se, m1.coef[5:1] - 2
                                         * m1.se[5:1]), 
                                         col = rgb(.6, .6, .6, 0.5), 
                                         border = NA)
polygon(x=c(year[1:5], year[5:1]), y = c(m2.coef + 2 * m2.se, m2.coef[5:1] - 2
                                         * m2.se[5:1]), 
                                         col = rgb(0, 0, 0.8, 0.5), 
                                         border = NA)
polygon(x=c(year[1:5], year[5:1]), y = c(m3.coef + 2 * m3.se, m3.coef[5:1] - 2
                                         * m3.se[5:1]), 
                                         col = rgb(0, 0, 1, 0.5),
                                         border = NA)
lines(year, m1.coef, type = "l", col = "black", lwd = 5)
lines(year, m2.coef, type = "l", col = "blue", lwd = 5)
lines(year, m3.coef, type = "l", col = "blue", lwd = 5)


###############################################################################

# Figure 6

# Excluding PML-only plots from control group

data.PMI.flat <- subset(data.flat, PML == 0 | (PML == 1 & PMI == 1)) 

# Models

fig6.a <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + 
             pop.den.5km.pre + dist.disturb.pre + A128 + A129 + A130 +
             A154 + A347 + A361 + A362, data = data.PMI.flat)
summary(fig6.a)

make(data.PMI.flat, c(fl.5yr.postPML, fl.yr1, fl.yr2, fl.yr3, fl.yr4, 
                      fl.yr5, PMI, fl.pre.5cell, forest.pre, pop.den.5km.pre,
                      dist.disturb.pre, A128, A129, A130, A154, A347, A361, 
                      A362, id))
out.PMI.flat <- matchit(PMI ~ fl.pre.5cell + forest.pre + pop.den.5km.pre +
                        dist.disturb.pre + A128 + A129 + A130 + A154 + A347 +
                        A361 + A362, method = "genetic", int.seed = 101, 
                        unif.seed = 101, pop.size = 500, wait.generations = 20,
                        MemoryMatrix = FALSE, ties = FALSE, 
                        discard = "control", data=data.make)

# Note: matching takes approximately 48 hours on standard-build desktop
# computer at time of writing

# Information used for SI Table B1

summary(out.PMI.flat, interactions = T, standardize = T) 
data.PMI.flat.matched <- match.data(out.PMI.flat)

# Note: data.PMI.flat.matched saved as "PMI_flat_matched_140414.csv" and 
# available from authors on request

fig6.b <- lm(fl.5yr.postPML ~ PMI, weights = weights, 
             data = data.PMI.flat.matched)
summary(fig6.b)

fig6.c <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + 
             pop.den.5km.pre + dist.disturb.pre + A128 + A129 + 
             A130 + A154 + A347 + A361 + A362, weights = weights, 
             data = data.PMI.flat.matched)
summary(fig6.c)

fig6 <- multiplot(fig6.a, fig6.b, fig6.c, coefficients = "PMI", 
                  names = c("(a)", "(b)", "(c)"), title = NULL, 
                  ylab = NULL, 
                  xlab = expression(paste("Effect on deforestation over 5 years (",
                         m^{2} , " per " , km^{2}, ")")),
                  pointSize = 5, lwdOuter = 2, horizontal = TRUE, 
                  secret.weapon = TRUE)
fig6 + theme_bw(base_size=18)


###############################################################################

# Figure 7

# Model 1 is year-on-year decomposition of Model (a) in Figure 6

m1.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, data = data.PMI.flat)
summary(m1.1)

m1.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, data = data.PMI.flat)
summary(m1.2)

m1.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, data = data.PMI.flat)
summary(m1.3)

m1.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, data = data.PMI.flat)
summary(m1.4)

m1.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, data = data.PMI.flat)
summary(m1.5)

# Model 2 is year-on-year decomposition of Model (b) in Figure 6

m2.1 <- lm(fl.yr1 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.1)

m2.2 <- lm(fl.yr2 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.2)

m2.3 <- lm(fl.yr3 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.3)

m2.4 <- lm(fl.yr4 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.4)

m2.5 <- lm(fl.yr5 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.5)

# Model 3 is year-on-year decomposition of Model (c) in Figure 6

m3.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, weights=weights,data=data.PMI.flat.matched)
summary(m3.1)

m3.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, weights = weights, data = data.PMI.flat.matched)
summary(m3.2)

m3.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, weights = weights, data = data.PMI.flat.matched)
summary(m3.3)

m3.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
           A362, weights = weights, data = data.PMI.flat.matched)
summary(m3.4)

m3.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
           A362, weights = weights, data = data.PMI.flat.matched)
summary(m3.5)

##Year by year effects of PMI

year <- seq(1:5)
m1.coef <- c(m1.1$coef[2], m1.2$coef[2], m1.3$coef[2], m1.4$coef[2], 
             m1.5$coef[2])
m1.se <- c(coef(summary(m1.1))[2, 2], coef(summary(m1.2))[2, 2], 
           coef(summary(m1.3))[2, 2], coef(summary(m1.4))[2, 2], 
           coef(summary(m1.5))[2, 2])
m2.coef <- c(m2.1$coef[2], m2.2$coef[2], m2.3$coef[2], m2.4$coef[2], 
             m2.5$coef[2])
m2.se <- c(coef(summary(m2.1))[2, 2], coef(summary(m2.2))[2, 2], 
           coef(summary(m2.3))[2, 2], coef(summary(m2.4))[2, 2], 
           coef(summary(m2.5))[2, 2])
m3.coef <- c(m3.1$coef[2], m3.2$coef[2], m3.3$coef[2], m3.4$coef[2],
             m3.5$coef[2])
m3.se <- c(coef(summary(m3.1))[2, 2], coef(summary(m3.2))[2, 2], 
           coef(summary(m3.3))[2, 2], coef(summary(m3.4))[2, 2], 
           coef(summary(m3.5))[2, 2])

par(mar = c(4.1, 5.1, 1.1,1.1))
plot(1, type = "n", xlim = c(0.8, 5.2), ylim = c(-500, 400), 
     ylab = expression(paste("Effect on Deforestation (", m^{2} , 
                             " per ", km^{2}, ")")), 
     xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
abline(h = 0, lty = 2)
polygon(x = c(year[1:5], year[5:1]),
        y = c(m1.coef + 2 * m1.se, m1.coef[5:1] - 2 * m1.se[5:1]),
        col = rgb(.6, .6, .6, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m2.coef + 2 * m2.se, m2.coef[5:1] - 2 * m2.se[5:1]), 
        col = rgb(0, 0, 0.8, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m3.coef + 2 * m3.se, m3.coef[5:1] - 2 * m3.se[5:1]), 
        col = rgb(0, 0, 1, 0.5), border = NA)
lines(year, m1.coef, type = "l", col = "black", lwd = 5)
lines(year, m2.coef, type = "l", col = "blue", lwd =  5)
lines(year, m3.coef, type = "l", col = "blue", lwd = 5)


###############################################################################

# Figure B1 (Supporting Information)

data.PML.flat.Shaur <- subset(data.flat, (PML == 0 | (PML == 1 & PMI == 0)) 
                              & A154 == 1)

m1s.PML <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + forest.pre + 
             pop.den.5km.pre + dist.disturb.pre + A128 + A129 + 
             A130 + A347 + A361 + A362, data = data.PML.flat.Shaur)
summary(m1s.PML)

m1s.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + 
               pop.den.5km.pre + dist.disturb.pre + A128 + 
               A129 + A130 + A347 + A361 + A362, data = data.PML.flat.Shaur)
summary(m1s.PML1)

m1s.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
               dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362, 
               data = data.PML.flat.Shaur)
summary(m1s.PML2)

m1s.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
               dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362, 
               data = data.PML.flat.Shaur)
summary(m1s.PML3)

m1s.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
               dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
               data = data.PML.flat.Shaur)
summary(m1s.PML4)

m1s.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
               dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
               data = data.PML.flat.Shaur)
summary(m1s.PML5)

make(data.PML.flat.Shaur, c(fl.5yr.postPML, fl.yr1, fl.yr2, fl.yr3, fl.yr4, 
                            fl.yr5, PML, fl.pre.5cell, forest.pre, 
                            pop.den.5km.pre, dist.disturb.pre, A128, A129, 
                            A130, A347, A361, A362, id, x, y))
out.PML.s.flat <- matchit(PML ~ fl.pre.5cell + forest.pre + pop.den.5km.pre + 
                          dist.disturb.pre + A128 + A129 + A130 + A347 + A361 +
                          A362, method = "genetic", int.seed = 102, 
                          unif.seed = 102, pop.size = 2000, 
                          wait.generations = 20, MemoryMatrix = FALSE, 
                          ties = FALSE, discard = "control", data = data.make)
sum.s.PML <- summary(out.PML.s.flat, interactions = T, standardize = T)
data.PMLonly.s.flat.matched <- match.data(out.PML.s.flat)

# write.csv(data.PMLonly.s.flat.matched,quote=F, 
#   file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PMLonly_Shaur_flat_matched_140422.csv", row.names=F)

m2s.PML <- lm(fl.5yr.postPML ~ PML, weights = weights, 
              data = data.PMLonly.s.flat.matched)
summary(m2s.PML)

m2s.PML1 <- lm(fl.yr1 ~ PML, weights = weights, 
               data = data.PMLonly.s.flat.matched)
summary(m2s.PML1)

m2s.PML2 <- lm(fl.yr2 ~ PML, weights = weights, 
               data = data.PMLonly.s.flat.matched)
summary(m2s.PML2)

m2s.PML3 <- lm(fl.yr3 ~ PML, weights = weights, 
               data = data.PMLonly.s.flat.matched)
summary(m2s.PML3)

m2s.PML4 <- lm(fl.yr4 ~ PML, weights = weights,
               data = data.PMLonly.s.flat.matched)
summary(m2s.PML4)

m2s.PML5 <- lm(fl.yr5 ~ PML, weights = weights,
               data = data.PMLonly.s.flat.matched)
summary(m2s.PML5)

m3s.PML <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + forest.pre +
              pop.den.5km.pre + dist.disturb.pre + A128 + A129 +
              A130 + A347 + A361 + A362,
              weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML)

m3s.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + 
               pop.den.5km.pre + dist.disturb.pre + A128 +
               A129 + A130 + A347 + A361 + A362, 
               weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML1)

m3s.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre +
               pop.den.5km.pre + dist.disturb.pre + A128 + 
               A129 + A130 + A347 + A361 + A362, 
               weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML2)

m3s.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + 
               pop.den.5km.pre + dist.disturb.pre + A128 +
               A129 + A130 + A347 + A361 + A362, 
               weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML3)

m3s.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + 
               pop.den.5km.pre + dist.disturb.pre + A128 + 
               A129 + A130 + A347 + A361 + A362,
               weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML4)

m3s.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + 
               pop.den.5km.pre + dist.disturb.pre + A128 + 
               A129 + A130 + A347 + A361 + A362,
               weights = weights, data = data.PMLonly.s.flat.matched)
summary(m3s.PML5)

B1 <- multiplot(m1s.PML, m2s.PML, m3s.PML, coefficients = "PML",
                names = c("(a)", "(b)", "(c)"), title = NULL, 
                ylab = NULL, 
                xlab = expression(paste("Effect on deforestation over 5 years (",
                                        m^{2} , " per " , km^{2}, ")")),
                pointSize = 5, lwdOuter = 2, horizontal = TRUE, 
                secret.weapon = TRUE)
B1 + theme_bw(base_size = 18) 

# Note: save as 600X400


###############################################################################

# Figure B2 (Supporting Information)

year <- seq(1:5)
m1s.coef <- c(m1s.PML1$coef[2], m1s.PML2$coef[2], m1s.PML3$coef[2], 
              m1s.PML4$coef[2], m1s.PML5$coef[2])
m1s.se <- c(coef(summary(m1s.PML1))[2, 2], coef(summary(m1s.PML2))[2, 2],
            coef(summary(m1s.PML3))[2, 2], coef(summary(m1s.PML4))[2, 2],
            coef(summary(m1s.PML5))[2, 2])
m2s.coef <- c(m2s.PML1$coef[2], m2s.PML2$coef[2], m2s.PML3$coef[2], 
              m2s.PML4$coef[2], m2s.PML5$coef[2])
m2s.se <- c(coef(summary(m2s.PML1))[2, 2], coef(summary(m2s.PML2))[2, 2], 
            coef(summary(m2s.PML3))[2, 2], coef(summary(m2s.PML4))[2, 2], 
            coef(summary(m2s.PML5))[2, 2])
m3s.coef <- c(m3s.PML1$coef[2], m3s.PML2$coef[2], m3s.PML3$coef[2], 
              m3s.PML4$coef[2], m3s.PML5$coef[2])
m3s.se <- c(coef(summary(m3s.PML1))[2, 2], coef(summary(m3s.PML2))[2, 2], 
            coef(summary(m3s.PML3))[2, 2], coef(summary(m3s.PML4))[2, 2], 
            coef(summary(m3s.PML5))[2, 2])

par(mar = c(4.1, 5.1, 1.1, 1.1))
plot(1, type = "n", xlim = c(0.8, 5.2), ylim = c(-350, 650), 
     ylab = expression(paste("Effect on Deforestation (", 
                             m^{2} , " per " , km^{2}, ")")),
     xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
abline(h = 0, lty = 2)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m1s.coef + 2 * m1s.se, m1s.coef[5:1] - 2 * m1s.se[5:1]), 
        col = rgb(.6, .6, .6, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m2s.coef + 2 * m2s.se, m2s.coef[5:1] - 2 * m2s.se[5:1]), 
        col = rgb(0, 0, 0.8, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m3s.coef + 2 * m3s.se, m3s.coef[5:1] - 2 * m3s.se[5:1]), 
        col = rgb(0, 0, 1, 0.5), border = NA)
lines(year, m1s.coef, type = "l", col = "black", lwd = 5)
lines(year, m2s.coef, type = "l", col = "blue", lwd = 5)
lines(year, m3s.coef, type = "l", col = "blue", lwd = 5)

# Save as 800x600 PNG file


###############################################################################

# Figure B3 (Supporting Information)

data.PMI.flat.Shaur <- subset(data.flat, 
                              (PML == 0 | (PML == 1 & PMI == 1)) & A154 == 1)

# Models

m1s <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
          dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
          data = data.PMI.flat.Shaur)
summary(m1s)

m1s.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            data = data.PMI.flat.Shaur)
summary(m1s.1)

m1s.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            data = data.PMI.flat.Shaur)
summary(m1s.2)

m1s.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            data = data.PMI.flat.Shaur)
summary(m1s.3)

m1s.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            data = data.PMI.flat.Shaur)
summary(m1s.4)

m1s.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            data = data.PMI.flat.Shaur)
summary(m1s.5)

make(data.PMI.flat.Shaur, 
     c(fl.5yr.postPML, fl.yr1, fl.yr2, fl.yr3, fl.yr4, fl.yr5, PMI, 
       fl.pre.5cell, forest.pre, pop.den.5km.pre, dist.disturb.pre, 
       A128, A129, A130, A347, A361, A362, id, x, y))
out.PMI.s.flat <- matchit(PMI ~ fl.pre.5cell + forest.pre + pop.den.5km.pre +
                          dist.disturb.pre + A128 + A129 + A130 + A347 + A361 +
                          A362, method = "genetic", int.seed = 101, 
                          unif.seed = 101, pop.size = 500, 
                          wait.generations = 20,
                          MemoryMatrix = FALSE, ties = FALSE, 
                          discard = "control", data = data.make)

# Note: approximately 48 hours on iMac to reach matching solution with this
# setup, 2-Mar-2014

sum.s.PMI <- summary(out.PMI.s.flat, interactions = T, standardize = T)
data.PMI.s.flat.matched <- match.data(out.PMI.s.flat)

# write.csv(data.PMI.s.flat.matched,quote = F, 
# file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PMI_Shaur_flat_matched_140422.csv", row.names = F)

m2s <- lm(fl.5yr.postPML ~ PMI, weights = weights, 
          data = data.PMI.s.flat.matched)
summary(m2s)

m2s.1 <- lm(fl.yr1 ~ PMI,weights = weights, data = data.PMI.s.flat.matched)
summary(m2s.1)

m2s.2 <- lm(fl.yr2 ~ PMI,weights = weights, data = data.PMI.s.flat.matched)
summary(m2s.2)

m2s.3 <- lm(fl.yr3 ~ PMI,weights = weights, data = data.PMI.s.flat.matched)
summary(m2s.3)

m2s.4 <- lm(fl.yr4 ~ PMI,weights = weights, data = data.PMI.s.flat.matched)
summary(m2s.4)

m2s.5 <- lm(fl.yr5 ~ PMI,weights = weights, data = data.PMI.s.flat.matched)
summary(m2s.5)

m3s <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
          dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
          weights = weights, data = data.PMI.s.flat.matched)
summary(m3s)

m3s.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            weights = weights, data = data.PMI.s.flat.matched)
summary(m3s.1)

m3s.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            weights = weights, data = data.PMI.s.flat.matched)
summary(m3s.2)

m3s.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            weights = weights, data = data.PMI.s.flat.matched)
summary(m3s.3)

m3s.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            weights = weights, data = data.PMI.s.flat.matched)
summary(m3s.4)

m3s.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
            dist.disturb.pre + A128 + A129 + A130 + A347 + A361 + A362,
            weights = weights, data = data.PMI.s.flat.matched)
summary(m3s.5)

B3 <- multiplot(m1s, m2s, m3s, coefficients = "PMI", 
                names = c("(a)", "(b)", "(c)"), 
                title = NULL, ylab = NULL, 
                xlab = expression(paste("Effect on deforestation over 5 years (",
                                        m^{2} , " per " , km^{2}, ")")),
                pointSize = 5, lwdOuter = 2, horizontal = TRUE, 
                secret.weapon = TRUE)
B3 + theme_bw(base_size = 18) #Note: save as 600X400


###############################################################################

# Figure B4 (Supporting Information)

year <- seq(1:5)
m1s.coef <- c(m1s.1$coef[2], m1s.2$coef[2], m1s.3$coef[2], m1s.4$coef[2],
              m1s.5$coef[2])
m1s.se <- c(coef(summary(m1s.1))[2, 2], coef(summary(m1s.2))[2, 2], 
            coef(summary(m1s.3))[2, 2], coef(summary(m1s.4))[2, 2], 
            coef(summary(m1s.5))[2, 2])
m2s.coef <- c(m2s.1$coef[2], m2s.2$coef[2], m2s.3$coef[2], m2s.4$coef[2],
              m2s.5$coef[2])
m2s.se <- c(coef(summary(m2s.1))[2, 2], coef(summary(m2s.2))[2, 2], 
            coef(summary(m2s.3))[2, 2], coef(summary(m2s.4))[2, 2], 
            coef(summary(m2s.5))[2, 2])
m3s.coef <- c(m3s.1$coef[2], m3s.2$coef[2], m3s.3$coef[2], m3s.4$coef[2],
              m3s.5$coef[2])
m3s.se <- c(coef(summary(m3s.1))[2, 2], coef(summary(m3s.2))[2, 2], 
            coef(summary(m3s.3))[2, 2], coef(summary(m3s.4))[2, 2], 
            coef(summary(m3s.5))[2, 2])

par(mar = c(4.1, 5.1, 1.1, 1.1))
plot(1, type = "n", xlim = c(0.8, 5.2), ylim = c(-500, 550),
     ylab = expression(paste("Effect on Deforestation (",
                             m^{2} , " per " , km^{2}, ")")),
     xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
abline(h = 0,lty = 2)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m1s.coef + 2* m1s.se, m1s.coef[5:1] - 2 * m1s.se[5:1]),
        col = rgb(.6, .6, .6, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m2s.coef + 2* m2s.se, m2s.coef[5:1] - 2 * m2s.se[5:1]), 
        col = rgb(0, 0, 0.8, 0.5), border = NA)
polygon(x = c(year[1:5], year[5:1]), 
        y = c(m3s.coef  +2* m3s.se, m3s.coef[5:1] - 2 * m3s.se[5:1]), 
        col = rgb(0, 0, 1, 0.5), border = NA)
lines(year, m1s.coef, type = "l", col = "black", lwd = 5)
lines(year, m2s.coef, type = "l", col = "blue", lwd = 5)
lines(year, m3s.coef, type = "l", col = "blue", lwd = 5)


###############################################################################

# Setting up for EVI analysis

names(data)[39:51] <- c("EVI_2000", "EVI_2001", "EVI_2002", "EVI_2003", 
                        "EVI_2004", "EVI_2005", "EVI_2006", "EVI_2007", 
                        "EVI_2008", "EVI_2009", "EVI_2010", "EVI_2011", 
                        "EVI_2012")

data.evi <- subset(data, A360 == 6 | A360 == 12 | A360 == 13)

# Flattening waves for analysis of 5-year effects

wave02 <- subset(data.evi, (A352 == 2 | A352 == 0))
wave02$EVI_pre <- wave02$EVI_2001
wave02$EVI_t0 <- wave02$EVI_2002
wave02$EVI_t5 <- wave02$EVI_2007
wave02$EVI_5yr_diff <- wave02$EVI_t5 - wave02$EVI_t0
wave02$EVI.D.yr1 <- wave02$EVI_2003 - wave02$EVI_2002
wave02$EVI.D.yr2 <- wave02$EVI_2004 - wave02$EVI_2003
wave02$EVI.D.yr3 <- wave02$EVI_2005 - wave02$EVI_2004
wave02$EVI.D.yr4 <- wave02$EVI_2006 - wave02$EVI_2005
wave02$EVI.D.yr5 <- wave02$EVI_2007 - wave02$EVI_2006
wave02$EVI.D.pre.5cell <- wave02$A474 - wave02$A473
wave02$forest.pre <- wave02[, 314] - wave02[, 287]
wave02$pop.den.5km.pre <- wave02$A176
wave02$dist.disturb.pre <- wave02$A116
wave02$wave = "two"
wave02$rain_5yr_diff <- wave02$A373 - wave02$A368
wave02$rain_yr1_diff <- wave02$A369 - wave02$A368
wave02$rain_yr2_diff <- wave02$A370 - wave02$A369
wave02$rain_yr3_diff <- wave02$A371 - wave02$A370
wave02$rain_yr4_diff <- wave02$A372 - wave02$A371
wave02$rain_yr5_diff <- wave02$A373 - wave02$A372
wave02 <- subset(wave02, forest.pre > 450000 & A346 == 0) 

wave03 <- subset(data.evi,(A352 == 3 | A352 == 0))
wave03$EVI_pre <- wave03$EVI_2002
wave03$EVI_t0 <- wave03$EVI_2003
wave03$EVI_t5 <- wave03$EVI_2008
wave03$EVI_5yr_diff <- wave03$EVI_t5 - wave03$EVI_t0
wave03$EVI.D.yr1 <- wave03$EVI_2004 - wave03$EVI_2003
wave03$EVI.D.yr2 <- wave03$EVI_2005 - wave03$EVI_2004
wave03$EVI.D.yr3 <- wave03$EVI_2006 - wave03$EVI_2005
wave03$EVI.D.yr4 <- wave03$EVI_2007 - wave03$EVI_2006
wave03$EVI.D.yr5 <- wave03$EVI_2008 - wave03$EVI_2007
wave03$EVI.D.pre.5cell <- wave03$A475 - wave03$A474
wave03$forest.pre <- wave03[, 314] - rowSums(wave03[, 287:288])
wave03$pop.den.5km.pre <- wave03$A181
wave03$dist.disturb.pre <- wave03$A117
wave03$wave = "three"
wave03$rain_5yr_diff <- wave03$A374 - wave03$A369
wave03$rain_yr1_diff <- wave03$A370 - wave03$A369
wave03$rain_yr2_diff <- wave03$A371 - wave03$A370
wave03$rain_yr3_diff <- wave03$A372 - wave03$A371
wave03$rain_yr4_diff <- wave03$A373 - wave03$A372
wave03$rain_yr5_diff <- wave03$A374 - wave03$A373
wave03 <- subset(wave03, forest.pre > 450000 & A346 == 0)

wave05 <- subset(data.evi,(A352 == 5 | A352 == 0))
wave05$EVI_pre <- wave05$EVI_2004
wave05$EVI_t0 <- wave05$EVI_2005
wave05$EVI_t5 <- wave05$EVI_2010
wave05$EVI_5yr_diff <- wave05$EVI_t5 - wave05$EVI_t0
wave05$EVI.D.yr1 <- wave05$EVI_2006 - wave05$EVI_2005
wave05$EVI.D.yr2 <- wave05$EVI_2007 - wave05$EVI_2006
wave05$EVI.D.yr3 <- wave05$EVI_2008 - wave05$EVI_2007
wave05$EVI.D.yr4 <- wave05$EVI_2009 - wave05$EVI_2008
wave05$EVI.D.yr5 <- wave05$EVI_2010 - wave05$EVI_2009
wave05$EVI.D.pre.5cell <- wave05$A477 - wave05$A476
wave05$forest.pre <- wave05[, 314] - rowSums(wave05[, 287:290])
wave05$pop.den.5km.pre <- wave05$A191
wave05$dist.disturb.pre <- wave05$A119
wave05$wave = "five"
wave05$rain_5yr_diff <- wave05$A376 - wave05$A371
wave05$rain_yr1_diff <- wave05$A372 - wave05$A371
wave05$rain_yr2_diff <- wave05$A373 - wave05$A372
wave05$rain_yr3_diff <- wave05$A374 - wave05$A373
wave05$rain_yr4_diff <- wave05$A375 - wave05$A374
wave05$rain_yr5_diff <- wave05$A376 - wave05$A375
wave05 <- subset(wave05, forest.pre > 450000 & A346 == 0)

wave06 <- subset(data.evi, (A352 == 6 | A352 == 0))
wave06$EVI_pre <- wave06$EVI_2005
wave06$EVI_t0 <- wave06$EVI_2006
wave06$EVI_t5 <- wave06$EVI_2011
wave06$EVI_5yr_diff <- wave06$EVI_t5 - wave06$EVI_t0
wave06$EVI.D.yr1 <- wave06$EVI_2007 - wave06$EVI_2006
wave06$EVI.D.yr2 <- wave06$EVI_2008 - wave06$EVI_2007
wave06$EVI.D.yr3 <- wave06$EVI_2009 - wave06$EVI_2008
wave06$EVI.D.yr4 <- wave06$EVI_2010 - wave06$EVI_2009
wave06$EVI.D.yr5 <- wave06$EVI_2011 - wave06$EVI_2010
wave06$EVI.D.pre.5cell <- wave06$A478 - wave06$A477
wave06$forest.pre <- wave06[, 314] - rowSums(wave06[, 287:291])
wave06$pop.den.5km.pre <- wave06$A196
wave06$dist.disturb.pre <- wave06$A120
wave06$wave ="six"
wave06$rain_5yr_diff <- wave06$A377 - wave06$A372
wave06$rain_yr1_diff <- wave06$A373 - wave06$A372
wave06$rain_yr2_diff <- wave06$A374 - wave06$A373
wave06$rain_yr3_diff <- wave06$A375 - wave06$A374
wave06$rain_yr4_diff <- wave06$A376 - wave06$A375
wave06$rain_yr5_diff <- wave06$A377 - wave06$A376
wave06 <- subset(wave06, forest.pre > 450000 & A346 == 0)

wave07 <- subset(data.evi,(A352 == 7 | A352 == 0))
wave07$EVI_pre <- wave07$EVI_2006
wave07$EVI_t0 <- wave07$EVI_2007
wave07$EVI_t5 <- wave07$EVI_2012
wave07$EVI_5yr_diff <- wave07$EVI_t5 - wave07$EVI_t0
wave07$EVI.D.yr1 <- wave07$EVI_2008 - wave07$EVI_2007
wave07$EVI.D.yr2 <- wave07$EVI_2009 - wave07$EVI_2008
wave07$EVI.D.yr3 <- wave07$EVI_2010 - wave07$EVI_2009
wave07$EVI.D.yr4 <- wave07$EVI_2011 - wave07$EVI_2010
wave07$EVI.D.yr5 <- wave07$EVI_2012 - wave07$EVI_2011
wave07$EVI.D.pre.5cell <- wave07$A479 - wave07$A478
wave07$forest.pre <- wave07[, 314] - rowSums(wave07[, 287:292])
wave07$pop.den.5km.pre <- wave07$A201
wave07$dist.disturb.pre <- wave07$A121
wave07$wave = "seven"
wave07$rain_5yr_diff <- wave07$A378 - wave07$A373
wave07$rain_yr1_diff <- wave07$A374 - wave07$A373
wave07$rain_yr2_diff <- wave07$A375 - wave07$A374
wave07$rain_yr3_diff <- wave07$A376 - wave07$A375
wave07$rain_yr4_diff <- wave07$A377 - wave07$A376
wave07$rain_yr5_diff <- wave07$A378 - wave07$A377
wave07 <- subset(wave07, forest.pre>450000 & A346 == 0) 

data.flat.evi <- rbind(wave02, wave03, wave05, wave06, wave07)


###############################################################################

# Figure B5 (Supporting Information)

# Exact matching on these ecoregions and protected status by subsetting

data.PML.flat <- subset(data.flat.evi,PMI == 0 & (A360 == 12 | A360 == 13) & 
                        wave == "two" & A347 == 0) 
data.PML.flat$EVI_5yr_down1000 <- 
  ifelse(data.PML.flat$EVI_5yr_diff <= -1000, 1, 0)
data.PML.flat$EVI_5yr_down1500 <- 
  ifelse(data.PML.flat$EVI_5yr_diff <= -1500, 1, 0)
data.PML.flat$EVI_5yr_down2000 <- 
  ifelse(data.PML.flat$EVI_5yr_diff <= -2000, 1, 0)
data.PML.flat$A360 <- as.factor(data.PML.flat$A360)

# Models

m1.pml <- glm(EVI_5yr_down1000 ~ PML + EVI.D.pre.5cell +
              EVI_t0 + pop.den.5km.pre + dist.disturb.pre +
              A128 + A129 + A130 + A154 + A360 + A361 + A362,
              family = binomial, data = data.PML.flat)
summary(m1.pml)

make(data.PML.flat, c(EVI_5yr_diff, EVI.D.yr1, EVI.D.yr2, 
                      EVI.D.yr3, EVI.D.yr4, EVI.D.yr5, PML, 
                      EVI.D.pre.5cell, EVI_t0, rain_5yr_diff,
                      rain_yr1_diff, rain_yr2_diff, rain_yr3_diff, 
                      rain_yr4_diff, rain_yr5_diff, forest.pre, 
                      pop.den.5km.pre, dist.disturb.pre, A128, 
                      A129, A130, A154, A347, A360, A361, A362, id))
out.PML.flat.evi <- matchit(PML ~ EVI.D.pre.5cell + EVI_t0 + 
                            pop.den.5km.pre + dist.disturb.pre +
                            A128 + A129 + A130 + A154 + A360 + A361 + 
                            A362, method = "genetic", int.seed = 101,
                            unif.seed = 101, pop.size = 500,
                            wait.generations = 20, MemoryMatrix = FALSE,
                            ties = FALSE, discard = "control", 
                            data = data.make)
summary(out.PML.flat.evi, interactions = T, standardize = T)
data.PML.flat.evi.matched <- match.data(out.PML.flat.evi)

# write.csv(data.PML.flat.evi.matched,
#           quote = F, 
#           file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PML_flatevi_matched_141102.csv", 
#           row.names = F)

data.PML.flat.evi.matched$EVI_5yr_down500 <- 
  ifelse(data.PML.flat.evi.matched$EVI_5yr_diff <= -500, 1, 0)
data.PML.flat.evi.matched$EVI_5yr_down1000 <- 
  ifelse(data.PML.flat.evi.matched$EVI_5yr_diff <= -1000, 1, 0)
data.PML.flat.evi.matched$EVI_5yr_down1500 <- 
  ifelse(data.PML.flat.evi.matched$EVI_5yr_diff <= -1500, 1, 0)
data.PML.flat.evi.matched$EVI_5yr_down2000 <- 
  ifelse(data.PML.flat.evi.matched$EVI_5yr_diff <= -2000, 1, 0)

m2.pml <- glm(EVI_5yr_down1000 ~ PML , family = binomial, 
              data = data.PML.flat.evi.matched, weights = weights)
summary(m2.pml)

m3.pml <- glm(EVI_5yr_down1000 ~ PML + EVI.D.pre.5cell + EVI_t0 + 
              pop.den.5km.pre + dist.disturb.pre + A128 + A129 + A130 +
              A154 + A360 + A361 + A362, family = binomial, 
              data = data.PML.flat.evi.matched, weights = weights)
summary(m3.pml)

B5 <- multiplot(m1.pml, m2.pml, m3.pml, 
                coefficients = "PML", 
                names = c("(a)", "(b)", "(c)"),
                title = NULL, ylab = NULL, 
                xlab = expression(paste("Log-Odds of 0.1 Decrease in EVI over 5 Years")),
                pointSize = 5, lwdOuter = 2, horizontal = TRUE, 
                secret.weapon = TRUE)
B5 + theme_bw(base_size = 18) # Note: save as 600X400


###############################################################################

# Figure B6 (Supporting Information)

# Exact matching on these ecoregions by subsetting

data.PMI.flat <- subset(data.flat.evi,(PML == 0 |PMI == 1) & 
                         (A360 == 6 | A360 == 12)) 
data.PMI.flat$EVI_5yr_down1000 <- 
  ifelse(data.PMI.flat$EVI_5yr_diff <= -1000, 1, 0)
data.PMI.flat$EVI_5yr_down1500 <- 
  ifelse(data.PMI.flat$EVI_5yr_diff <= -1500, 1, 0)
data.PMI.flat$EVI_5yr_down2000 <- 
  ifelse(data.PMI.flat$EVI_5yr_diff <= -2000, 1, 0)
data.PMI.flat$A360 <- as.factor(data.PMI.flat$A360)

m1 <- glm(EVI_5yr_down1000 ~ PMI + EVI.D.pre.5cell + EVI_t0 + pop.den.5km.pre +
          dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A360 + A361 + 
          A362 + wave, family = binomial, data = data.PMI.flat)
summary(m1)

wave03.pmi <- subset(wave03,(PML == 0 |PMI == 1) & (A360 == 6 | A360 == 12))
wave05.pmi <- subset(wave05,(PML == 0 |PMI == 1) & (A360 == 6 | A360 == 12))
wave06.pmi <- subset(wave06,(PML == 0 |PMI == 1) & (A360 == 6 | A360 == 12))
wave07.pmi <- subset(wave07,(PML == 0 |PMI == 1) & (A360 == 6 | A360 == 12))

make(wave03.pmi, c(EVI_5yr_diff, EVI.D.yr1, EVI.D.yr2, EVI.D.yr3, EVI.D.yr4, 
                   EVI.D.yr5, PMI, EVI.D.pre.5cell, EVI_t0, rain_5yr_diff, 
                   rain_yr1_diff, rain_yr2_diff, rain_yr3_diff, rain_yr4_diff, 
                   rain_yr5_diff, forest.pre, pop.den.5km.pre, 
                   dist.disturb.pre, A128, A129, A130, A154, A347, A360, A361,
                   A362, wave, id))
out.wave03.evi <- matchit(PMI ~ EVI.D.pre.5cell + EVI_t0 + pop.den.5km.pre + 
                          dist.disturb.pre + A128 + A129 + A130 + A154 + A347 +
                          A360 + A361 + A362, method = "genetic",
                          int.seed = 101, unif.seed = 101,pop.size = 500,
                          wait.generations = 20, MemoryMatrix = FALSE,
                          ties = FALSE, discard = "control", data = data.make)
summary(out.wave03.evi, interactions = T, standardize = T)
data.wave03.evi.matched <- match.data(out.wave03.evi)

make(wave05.pmi, c(EVI_5yr_diff, EVI.D.yr1, EVI.D.yr2, EVI.D.yr3, EVI.D.yr4, 
                   EVI.D.yr5, PMI, EVI.D.pre.5cell, EVI_t0, rain_5yr_diff, 
                   rain_yr1_diff, rain_yr2_diff, rain_yr3_diff, rain_yr4_diff,
                   rain_yr5_diff, forest.pre, pop.den.5km.pre, 
                   dist.disturb.pre, A128, A129, A130, A154, A347, A360, A361,
                   A362, wave, id))
out.wave05.evi <- matchit(PMI ~ EVI.D.pre.5cell + EVI_t0 + pop.den.5km.pre +
                          dist.disturb.pre + A128 + A129 + A130 + A154 + A347 +
                          A360 + A361 + A362, method = "genetic", 
                          int.seed = 101, unif.seed = 101, pop.size = 500, 
                          wait.generations = 20, MemoryMatrix = FALSE, 
                          ties = FALSE, discard = "control", data = data.make)
summary(out.wave05.evi, interactions = T, standardize = T)
data.wave05.evi.matched <- match.data(out.wave05.evi)

make(wave06.pmi, c(EVI_5yr_diff, EVI.D.yr1, EVI.D.yr2, EVI.D.yr3, EVI.D.yr4, 
                   EVI.D.yr5, PMI, EVI.D.pre.5cell, EVI_t0, rain_5yr_diff, 
                   rain_yr1_diff, rain_yr2_diff, rain_yr3_diff, rain_yr4_diff,
                   rain_yr5_diff, forest.pre, pop.den.5km.pre, 
                   dist.disturb.pre, A128, A129, A130, A154, A347, A360, A361, 
                   A362, wave, id))
out.wave06.evi <- matchit(PMI ~ EVI.D.pre.5cell + EVI_t0 + pop.den.5km.pre + 
                          dist.disturb.pre + A128 + A129 + A130 + A154 + A347 +
                            A360 + A361 + A362, method = "genetic", 
                          int.seed = 101, unif.seed = 101, pop.size = 500,
                          wait.generations = 20, MemoryMatrix = FALSE, 
                          ties = FALSE, discard = "control", data = data.make)
summary(out.wave06.evi, interactions = T, standardize = T)
data.wave06.evi.matched <- match.data(out.wave06.evi)

make(wave07.pmi, c(EVI_5yr_diff, EVI.D.yr1, EVI.D.yr2, EVI.D.yr3, EVI.D.yr4, 
                   EVI.D.yr5, PMI, EVI.D.pre.5cell, EVI_t0, rain_5yr_diff, 
                   rain_yr1_diff, rain_yr2_diff, rain_yr3_diff, rain_yr4_diff,
                   rain_yr5_diff, forest.pre, pop.den.5km.pre, 
                   dist.disturb.pre, A128, A129, A130, A154, A347, A360, A361,
                   A362, wave, id))
out.wave07.evi <- matchit(PMI ~ EVI.D.pre.5cell + EVI_t0 + pop.den.5km.pre + 
                          dist.disturb.pre + A128 + A129 + A130 + A154 + 
                          A347 + A360 + A361 + A362, method = "genetic",
                          int.seed = 101, unif.seed = 101, pop.size = 500,
                          wait.generations = 20, MemoryMatrix = FALSE,
                          ties = FALSE, discard = "control", data = data.make)
summary(out.wave07.evi, interactions = T, standardize = T)
data.wave07.evi.matched <- match.data(out.wave07.evi)

data.pmi.evi.matched <- rbind(data.wave03.evi.matched, data.wave05.evi.matched,
                              data.wave06.evi.matched, data.wave07.evi.matched)

# write.csv(data.pmi.evi.matched,
#           quote = F, 
#           file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PMI_wavesevi_matched_141110.csv", 
#           row.names = F)

data.pmi.evi.matched$A360 <- as.factor(data.pmi.evi.matched$A360)
data.pmi.evi.matched$EVI_5yr_down500 <- 
  ifelse(data.pmi.evi.matched$EVI_5yr_diff <= -500, 1, 0)
data.pmi.evi.matched$EVI_5yr_down1000 <- 
  ifelse(data.pmi.evi.matched$EVI_5yr_diff <= -1000, 1, 0)
data.pmi.evi.matched$EVI_5yr_down1500 <- 
  ifelse(data.pmi.evi.matched$EVI_5yr_diff <= -1500, 1, 0)
data.pmi.evi.matched$EVI_5yr_down2000 <- 
  ifelse(data.pmi.evi.matched$EVI_5yr_diff <= -2000, 1, 0)

m2.waves <- glm(EVI_5yr_down1000 ~ PMI,family = binomial, 
                data = data.pmi.evi.matched, weights = weights)
summary(m2.waves)

m3.waves <- glm(EVI_5yr_down1000 ~ PMI + EVI.D.pre.5cell + EVI_t0 +
                pop.den.5km.pre + dist.disturb.pre + A128 + A129 + 
                A130 + A154 + A347 + A360 + A361 + A362 + wave, 
                family = binomial, data = data.pmi.evi.matched,
                weights = weights)
summary(m3.waves)

B6 <- multiplot(m1, m2.waves, m3.waves, coefficients = "PMI", 
                names = c("(a)", "(b)", "(c)"), title = NULL,
                ylab = NULL,
                xlab = expression(paste("Log-Odds of 0.1 Decrease in EVI over 5 Years")),
                pointSize = 5, lwdOuter = 2, horizontal = TRUE, 
                secret.weapon = TRUE)
B6 + theme_bw(base_size = 18) #Note: save as 600X400


###############################################################################

# Figure B8 (Supporting Information)

data.PMLonly.flat <- subset(data.flat, PML == 0 | (PML == 1 & PMI == 0))
data.PMLonly.flat$toPeru.within10km <- 
  ifelse(data.PMLonly.flat$toPeru.km <= 10, 1, 0)

# Models for PML

m1.PML <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + forest.pre + 
             pop.den.5km.pre + dist.disturb.pre + A128 + A129 + A130 + A154 +
             A347 + A361 + A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML)
m1.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML1)
m1.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML2)
m1.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML3)
m1.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML4)
m1.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.km, data = data.PMLonly.flat)
summary(m1.PML5)

make(data.PMLonly.flat, c(fl.5yr.postPML, fl.yr1, fl.yr2, fl.yr3, fl.yr4,
                          fl.yr5, PML, fl.pre.5cell, forest.pre, 
                          pop.den.5km.pre, dist.disturb.pre, A128, 
                          A129, A130, A154, A347, A361, A362, 
                          toPeru.within10km, id))
out.PMLonly.flat <- matchit(PML ~ fl.pre.5cell + forest.pre + pop.den.5km.pre +
                            dist.disturb.pre + A128 + A129 + A130 + A154 + 
                            A347 + A361 + A362 + toPeru.within10km,
                            method = "genetic", int.seed = 102, 
                            unif.seed = 102, pop.size = 1000, 
                            wait.generations = 20, MemoryMatrix = FALSE,
                            ties = FALSE, discard = "control", 
                            data = data.make)

# Note: pop.size = 500, seed = 101 does not result in 0.1 Std. Mean Diff.
# on main covariates after 26 hours on iMac, March 2015

summary(out.PMLonly.flat, interactions = T, standardize = T)
data.PMLonly.flat.matched <- match.data(out.PMLonly.flat)

#write.csv(data.PMLonly.flat.matched,
#   quote = F, 
#   file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PMLonly_flat_matched_toPeruBin_150311.csv", 
#   row.names = F)

m2.PML <- lm(fl.5yr.postPML ~ PML, weights = weights,
             data = data.PMLonly.flat.matched)
summary(m2.PML)
m2.PML1 <- lm(fl.yr1 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML1)
m2.PML2 <- lm(fl.yr2 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML2)
m2.PML3 <- lm(fl.yr3 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML3)
m2.PML4 <- lm(fl.yr4 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML4)
m2.PML5 <- lm(fl.yr5 ~ PML, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m2.PML5)

m3.PML <- lm(fl.5yr.postPML ~ PML + fl.pre.5cell + forest.pre +
             pop.den.5km.pre + dist.disturb.pre + A128 + A129 + 
             A130 + A154 + A347 + A361 + A362 + toPeru.within10km,
             weights = weights, data = data.PMLonly.flat.matched)
summary(m3.PML)
m3.PML1 <- lm(fl.yr1 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362 + toPeru.within10km, weights = weights,
              data = data.PMLonly.flat.matched)
summary(m3.PML1)
m3.PML2 <- lm(fl.yr2 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
              A362 + toPeru.within10km, weights = weights,
              data = data.PMLonly.flat.matched)
summary(m3.PML2)
m3.PML3 <- lm(fl.yr3 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.within10km, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m3.PML3)
m3.PML4 <- lm(fl.yr4 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.within10km, weights = weights, 
              data = data.PMLonly.flat.matched)
summary(m3.PML4)
m3.PML5 <- lm(fl.yr5 ~ PML + fl.pre.5cell + forest.pre + pop.den.5km.pre +
              dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + 
              A362 + toPeru.within10km, weights = weights,
              data = data.PMLonly.flat.matched)
summary(m3.PML5)

B8 <- multiplot(m1.PML, m2.PML, m3.PML, coefficients = "PML", 
                names = c("(a)", "(b)", "(c)"), title = NULL,
                ylab = NULL,
                xlab = expression(paste("Effect on deforestation over 5 years (",
                                        m^{2} , " per " , km^{2}, ")")),
                pointSize = 5, lwdOuter = 2,
                horizontal = TRUE, secret.weapon = TRUE)
B8 + theme_bw(base_size = 18) # Note: save as 600X400


###############################################################################

# Figure B9 (Supporting Information)

data.PMI.flat <- subset(data.flat,PML == 0 | (PML == 1 & PMI == 1)) 

m1 <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + A362 +
           toPeru.km, data = data.PMI.flat)
summary(m1)
m1.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km, data = data.PMI.flat)
summary(m1.1)
m1.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km, data = data.PMI.flat)
summary(m1.2)
m1.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km, data = data.PMI.flat)
summary(m1.3)
m1.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km, data = data.PMI.flat)
summary(m1.4)
m1.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km, data = data.PMI.flat)
summary(m1.5)

make(data.PMI.flat, c(fl.5yr.postPML, fl.yr1, fl.yr2, fl.yr3, fl.yr4, 
                      fl.yr5, PMI, fl.pre.5cell, forest.pre, 
                      pop.den.5km.pre, dist.disturb.pre, A128, A129, 
                      A130, A154, A347, A361, A362, toPeru.km, id))
out.PMI.flat <- matchit(PMI ~ fl.pre.5cell + forest.pre + pop.den.5km.pre +
                        dist.disturb.pre + A128 + A129 + A130 + A154 + A347 +
                        A361 + A362 + toPeru.km,
                        method = "genetic", int.seed = 101, unif.seed = 101,
                        pop.size = 500, wait.generations = 20,
                        MemoryMatrix = FALSE, ties = FALSE, 
                        discard = "control", data = data.make)

# Note: approximately 100 hours on iMac to reach matching solution with this
# setup, 2-Mar-2015

summary(out.PMI.flat, interactions = T, standardize = T)
data.PMI.flat.matched <- match.data(out.PMI.flat)

#write.csv(data.PMI.flat.matched,
#  quote = F, 
#  file = "~/Dropbox/Remote Sensing & Aid/PSUR/PSUR/PMI_flat_matched_toPeru_150302.csv", 
#  row.names = F)

m2 <- lm(fl.5yr.postPML ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2)
m2.1 <- lm(fl.yr1 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.1)
m2.2 <- lm(fl.yr2 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.2)
m2.3 <- lm(fl.yr3 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.3)
m2.4 <- lm(fl.yr4 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.4) 
m2.5 <- lm(fl.yr5 ~ PMI, weights = weights, data = data.PMI.flat.matched)
summary(m2.5)

m3 <- lm(fl.5yr.postPML ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre +
           dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 + A362 +
           toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3)
m3.1 <- lm(fl.yr1 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3.1)
m3.2 <- lm(fl.yr2 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3.2)
m3.3 <- lm(fl.yr3 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3.3)
m3.4 <- lm(fl.yr4 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3.4)
m3.5 <- lm(fl.yr5 ~ PMI + fl.pre.5cell + forest.pre + pop.den.5km.pre + 
             dist.disturb.pre + A128 + A129 + A130 + A154 + A347 + A361 +
             A362 + toPeru.km,weights = weights,data = data.PMI.flat.matched)
summary(m3.5)

B9 <- multiplot(m1, m2, m3, coefficients = "PMI",
                names = c("(a)", "(b)", "(c)"),
                title = NULL, ylab = NULL,
                xlab = expression(paste("Effect on deforestation over 5 years (",
                                        m^{2} , " per " , km^{2}, ")")),
                pointSize = 5, lwdOuter = 2,
                horizontal = TRUE, secret.weapon = TRUE)
B9 + theme_bw(base_size = 18) # Note: save as 600X400
