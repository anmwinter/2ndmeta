# Frequenist to Bayesian Modeling of Growth Rates and Maximum Body Sizes
# Or an introduction to modeling 

# All original text is released under the CC-BY 4.0, Ara Kooser

# References
# If you got this from me or github then it comes with some data. If you use the data in a talk, workshop, or paper
# Please cite:
# Grady, John M., et al. "Evidence for mesothermy in neoaurs." Science 344.6189 (2014): 1268-1272.

# If you use the Venditti prior please cite:
# Sibly, Richard M., et al. "How body mass and lifestyle affect juvenile biomass production in placental mammals."
# Proceedings of the Royal Society of London B: Biological Sciences 281.1777 (2014): 20132818.

library(ggplot2)
library(lme4)
library(BayesFactor)
library(BayesianFirstAid)
library(MCMCglmm)
library(plyr)
library(MCMCpack)
library(scapeMCMC)

setwd("~/Desktop/phylo_glmm/")
neo <- read.csv("~/Desktop/phylo_glmm/neornithes_growth_data.csv", sep=",",na.strings = "")

head(neo)
tail(neo)
summary(neo)

# Scatterplot
# Make some notes on what each part of the does? Why did I choose Clade_3? 
# Is a scatter plot appropriate for your data? 
# Darker outer circle, lighter inside
maxG3.est_G3M.est <- ggplot(neo,aes(x=log(G3M.est), log(y=maxGG3))) + geom_point(aes(color=Clade_3)) +
  ggtitle("Log plot of G3M.est and maxGG3 for Neornithes") +
  xlab("log G3M")+
  ylab("log maxGG3")+
  scale_fill_hue(l=50)+
  guides(col = guide_legend(nrow = 16))+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )
maxG3.est_G3M.est

# What about a box plot?
# 

# Frequentist
# "Most ecologists use the frequentist approach. This approach focuses on P(D|H), the 
# probability of the data, given the hypothesis. That is, this approach treats 
# data as random (if you repeated the study, the data might come 
# out differently), and hypotheses as fixed (the hypothesis is either true or 
# false, and so has a probability of either 1 or 0, you 
# just don’t know for sure which it is)."

# “The world is a certain way, but I don’t know how it is. Further, I can’t necessarily tell how the world is
# just by collecting data, because data are always finite and noisy. So I’ll use statistics to line up the
# alternative possibilities, and see which ones the data more or less rule out.”
# -Brian McGill 

# CORRELATION

# Is there a correlation between GG3max and G3M?

# Correlation measures the degree to which two variables are related. Correlation does NOT fit a line through 
# your data points. You are computing a correlation coefficient (r) that tells you how much one variable 
# tends to change when the other one does. When r is 0.0, there is no relationship. When r is positive, there is 
# a trend that one variable goes up as the other one goes up. When r is negative, there is a trend that one 
# variable goes up as the other one goes down.

# Correlation computes the value of the Pearson correlation coefficient, r. Its value ranges from -1 to +1
# Assumes a bivariate normal distribution.

cor.test( ~ log(maxGG3) + log(G3M.est), data = neo)

# To dig further into correlations you need to meet the following assumptions:
# X and Y are measured
# Both are sampled from Gaussian distributions

# REGULAR LINEAR REGRESSION
# Linear regression finds the best line that predicts Y from X.
# Linear regression quantifies goodness of fit with r2, sometimes shown in uppercase as R2. If you put the 
# same data into correlation (which is rarely appropriate), the square of r from correlation will 
# equal r2 from regression. 

# The vertical distances of the points from the best-fit line (the residuals) are assumed to follow a Gaussian 
# distribution, with the SD of the scatter not related to the X or Y values.

# The model below is in R linear model language. It reads:
# Does G3M.est predict maxGG3 which are found the the neo dataframe

lm_neo = lm(log(maxGG3) ~ log(G3M.est), neo)

# NOTE!!! The above actually says this log(maxGG3) ~ log(G3M.est) + Epsilon
# ~ log(G3M.est) is our fixed term
# Epsilon is our error term that accounts for everything else
# So this reads log(maxGG3) is predicted by log(G3M.est) + some random effects

summary(lm_neo)

# Graphical check
maxG3.est_G3M.est + geom_abline(intercept = -1.1825, slope = 0.59605)

# Explanation of the output
# Multiple R-squared = R^2, the measure of variance accounted for
# Adjusted R-squared includes fixed effects and can be lower if you pile on fixed effects!

# Call:
#   lm(formula = log(maxGG3) ~ log(G3M.est), data = neo)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.93555 -0.37635  0.08668  0.46264  1.01321 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.18250    0.24872  -4.754 1.39e-05 ***
# log(G3M.est)  0.59605    0.03756  15.867  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6415 on 57 degrees of freedom
# Multiple R-squared:  0.8154,  Adjusted R-squared:  0.8122 
# F-statistic: 251.8 on 1 and 57 DF,  p-value: < 2.2e-16

# CHECKS!!!
# Check linearity see ggplot log vs log
# Check heteroskedasticity
plot(fitted(lm_neo),residuals(lm_neo))
# Slight heteroskedasticity
# Check heteroskedasticity against this: 
# plot(rnorm(100),rnorm(100))

# Normality of residuals
hist(residuals(lm_neo))
# Right skewed
qqnorm(residuals(lm_neo))
# Curved

# Absence of influential data points?
dfbeta(lm_neo)

# Independence? Check!

# LINEAR MIXED MODEL
# Here we are curious about modeling phylogenetics (as a categorical variable) added as a random effect
# What's a random effect? Well it's complicated. 

# Why a mixed model? Well, phylogenetics! 
# We can model phylogenetic differences by assuming different random intercepts for each clade.
# Each clade is assigned a different intercept value and the mixed model estimates these intercepts. 

# So in lm's we only have fixed effects and the error term is usually lacking structure (i.e. boring)
# In a lmm random effects essentially give structure to the error term epsilon. 

lmm <- lmer(log(maxGG3) ~ log(G3M.est) + (1|Clade_3), data=neo)
# Notice something weird here?
# (1|Clade_3), yeah that. This says assume an intercept that’s different for each category in Clade_3
# You might want to look at different clade levels

lmm
# Notice, no summary here! Just call the model

# Linear mixed model fit by REML ['lmerMod']
# Formula: log(maxGG3) ~ log(G3M.est) + (1 | Clade_3)
# Data: neo
# REML criterion at convergence: 70.0405
# Random effects:
# Groups   Name        Std.Dev.
# Clade_3  (Intercept) 0.5954  
# Residual             0.3348  
# Number of obs: 59, groups:  Clade_3, 10
# Fixed Effects:
# (Intercept)  log(G3M.est)  
# -2.0479        0.7543 

# Let's break this down. Check out the random effects first.
# Random effects:
# Groups   Name        Std.Dev.
# Clade_3  (Intercept) 0.5954  
# Residual             0.3348  

# The standard deviation is a measure of variability. Don't freak! Remember log scale!
# Residuals are the "random" deviations we haven't thought of

# Next let's look at the fixed effects
# Fixed Effects:
# (Intercept)  log(G3M.est)  
# -2.0479        0.7543 

# What's important here is the number under log(G3M.est). This is the slope for the categorical effects
# of phylogenetics. Compare this to the lm model which only looked at G3M.est.

# CHECKS!!!
# Check linearity see ggplot log vs log
# Check heteroskedasticity
plot(fitted(lmm),residuals(lmm))

# Normality of residuals
hist(residuals(lmm))

qqnorm(residuals(lmm))

# The profile function systematically varies the parameters in a model, assessing
# the best possible fit that can be obtained with one parameter fixed
# at a specific value and comparing this fit to the globally optimal fit, which is
# the original model fit that allowed all the parameters to vary.

pr01 <- profile(lmm)
# Can can make lots of fun plots here if you want to from the pr01
pr01


#Check the lmm line
maxG3.est_G3M.est + geom_abline(intercept = -2.0479, slope = 0.7543)

# COMPARING lm and lmer models
# Let's test the two models against each other.

model <- lmer(log(maxGG3) ~ log(G3M.est) + (1|Clade_3), data=neo,REML=FALSE) # with taxonomy
null_model <- lm(log(maxGG3) ~ log(G3M.est), neo) #no taxonomy, no phylogenetics

AIC(model,null_model)
# AIC lower is generally better but with multiple models of similar AIC you will want to do
# model averaging 

# Delta AIC
bbmle::AICtab(model,null_model)


# BAYESIAN
# "Bayesian statistics focuses on P(H|D), the probability of the hypothesis, given the data. That is, this approach
# treats the data as fixed (these are the only data you have) and hypotheses as random (the hypothesis might be
# true or false, with some probability between 0 and 1). This approach is called Bayesian because you need to use
# Bayes’ Theorem to calculate P(H|D)."

# A Bayesian basically says, “I don’t know how the world is. All I have to go on is finite data. So I’ll use
# statistics to infer something from those data about how probable different possible states of the world are.”
# -Brian McGill

# A note on priors. Prior distrubitions do NOT accurately reflect your state of knowledge of belief.
# The prior distrubtion is actually part of the model and can be checked (see our multiple models).
# Prior are in fact a mixture of "substantive knowledge, scientific conjectures, statistical properties,
# analytical convenience, disciplinary tradition, and computational tractability."
# Gelman and Shalizi, British Journal of Mathematical and Statustical Psychology, 2013, 66

# Your models should be checked against your data not necessarily against each other. 

# CORRELATION
# Bayesian First Aid Pearson's Correlation Coefficient Test
# Bivariate t distribution
# You get wider tails so it copes better with outliers 
# Built on a flat prior for details see: 
# http://www.sumsar.net/blog/2014/03/bayesian-first-aid-pearson-correlation-test/

fit_neo <- bayes.cor.test( ~ log(maxGG3) + log(G3M.est), data = neo)
fit_neo
plot(fit_neo)
# Shows you the distribution we would expect a new data point to have
# Light blue oval is the 95% density distrubution region
# Dark blue oval is the 50% density distribution region


## BAYESIAN LINEAR REGRESSION - no phylogenetics or taxonomy
model1 <- MCMCregress(log(maxGG3) ~ log(G3M.est), data = neo, b0 = 0, B0 = 0.1, c0 = 2, d0 = 0.11,
                      marginal.likelihood = "Chib95")
summary(model1)

# What this is saying:

# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
 
#   Mean      SD  Naive SE Time-series SE
# (Intercept)  -1.1725 0.24849 0.0024849      0.0024849
# log(G3M.est)  0.5946 0.03772 0.0003772      0.0003772
# sigma2        0.4143 0.07968 0.0007968      0.0008238

# I looked at 10,000 different possible regression lines and scored how likely each 
# one is. Then I averaged over all those to come up with the mean parameters, and I took
# std devs over them all to get the standard deviation of parameters (and so on for the 
# other statistics)

# So the single most probable line is given by the mean and can be read like this:
# y = 0.5946x + 0.4143 with a square residual of 0.4143
# The quartiles tell us that 97.5% of the lines tested have a postive slope. 
# In other words given a random line going through the data set it will be positive 
# 97.5% of the time. Our mean line occurs between 50-75% of the time.


# Check the mean Bayesian Regression Line
maxG3.est_G3M.est + geom_abline(intercept = -1.1725, slope = 0.5946)


## MCMCglmm - Bayesian mixed model with taxonomy

# Weakly informed Inverse Gamma prior with taxonomy as a randon effect
prior0<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
# G-structure is the covariance matrix of the random effects.
# R-structure is the covariance matrix of the residuals.
# nu is a degreee of belief, nu = 0 is an improper prior
# So what does a V of 1 mean in the above prior?

# NOTE: We did not specify B! What assumptions are being made here?
# B is for the fixed effects and is a list containing V and mu
# V is the expected covariance
# mu is expected value
# combined these give you a degree of belief 

model_simple0<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Clade_3,
                       family="gaussian",prior=prior0,
                       data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple0)

allChains <- as.mcmc(cbind(model_simple0$Sol,model_simple0$VCV))
plotTrace(allChains,axes=TRUE,las=1)
# You want this to look like white noise. 

# Check the Bayesian line with a weak prior with taxonomy as the random effect
maxG3.est_G3M.est + geom_abline(intercept = -2.0514, slope = 0.7545)

# MCMCglmm with phylogenetic tree
library(ape)

# Reads in the newick tree
phylo<-read.tree("bird.tre")
phylo

head(neo)

# Some magic happens here, what is it?
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)

# Ara's starting simple model with weakly informed Inverse Gamma prior with phylogenetics as a random effect
prior1<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
# Many thanks to Jarrod Hadfield for breaking this down

# We have a single random effect for each of species (n). They will have a correlation structure A. 
# The correlation structure A is an nun matrix which represents the proportion of evolution time shared 
# between each species. So the covariance structure is Va*A where Va is some scalar variance estimated 
# from the data. 

# units (see summary) is a unique factor in each row of the data frame. This will have a standard correlation structure I
# for the residuals. I is an identity matrix. This is because we believe that the residual are independent 
# and have the same finite variance. If there is only one observation per species then I is nxn. 
# The covariance structure is Ve*I where Ve is some scalar variance estimated from the data.  

model_simple1<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Species_binomial,
                       family="gaussian",ginverse=list(Species_binomial=inv.phylo$Ainv),prior=prior1,
                       data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple1)

allChains <- as.mcmc(cbind(model_simple1$Sol,model_simple1$VCV))
plotTrace(allChains,axes=TRUE,las=1)

# Check the Bayesian model_simple1
maxG3.est_G3M.est + geom_abline(intercept = -2.2931, slope = 0.7323)

# Model with phylogenetics and parameter expansion for G (phylogenetics)
# Weakly informed Inverse Gamma prior for R

prior2<-list(G=list(G1=list(V=1,nu=0.02,alpha.mu=0.01)),R=list(V=1,nu=0.02))
# alpha.mu is code for parameter-expanded priors.
# Main advantage of parameter expansion is algorithmic.
# It increases the efficiency of the MCMC algorithm when the posterior for G contains values close to zero.
# alpha.mu is the mean vector

# Parameter expansion
# Jarrod Hadfield - based on de Villemereuil et al. (2013) Methods in Ecology and Evolution 4:260-275

model_simple2<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Species_binomial,
                        family="gaussian",ginverse=list(Species_binomial=inv.phylo$Ainv),prior=prior2,
                        data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple2)

allChains <- as.mcmc(cbind(model_simple2$Sol,model_simple2$VCV))
plotTrace(allChains,axes=TRUE,las=1)

# Check the Bayesian model_simple2
maxG3.est_G3M.est + geom_abline(intercept = -2.2793, slope = 0.7289)

## With aplha.v specified 
## Model with phylogenetics and parameter expansion for G (phylogenetics)

prior3<-list(G=list(G1=list(V=1,nu=0.02,alpha.mu=0.02,alpha.V=1000)),R=list(V=1,nu=0.02))
# alpha.V is a covariance matrix 
# What values of alpha.V are reasonable? 

model_simple3<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Species_binomial,
                        family="gaussian",ginverse=list(Species_binomial=inv.phylo$Ainv),prior=prior3,
                        data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple3)

allChains <- as.mcmc(cbind(model_simple3$Sol,model_simple3$VCV))
plotTrace(allChains,axes=TRUE,las=1)

# Check the Bayesian model_simple3
maxG3.est_G3M.est + geom_abline(intercept = -2.2890, slope = 0.7301)

## Chris Venditti's Prior
# This should be fine for almost all comparative work that assumes a Gaussian error structure
prior4<-list(R=list(V=1,nu=0.02),B=list(mu=(rep(0,2)),V=diag(2)*1e10),
             G=list(G1=list(V=1,nu=1,alpha.mu=1000,alpha.V=1000)))
# B is fixed effects with mu and V specified

model_simple4<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Species_binomial,
                        family="gaussian",ginverse=list(Species_binomial=inv.phylo$Ainv),prior=prior4,
                        data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple4)

allChains <- as.mcmc(cbind(model_simple3$Sol,model_simple3$VCV))
plotTrace(allChains,axes=TRUE,las=1)

# Check the Bayesian model_simple4
maxG3.est_G3M.est + geom_abline(intercept = -2.3095, slope = 0.7305)

# Deviance information criterion
# Look for lower values to justify inclusions of terms and priors
model_simple0$DIC # weak prior with taxonomy
model_simple1$DIC # weak prior with phylogenetics
model_simple2$DIC # weak prior with parameter expansion and phylogenetics
model_simple3$DIC # weak prior with parameter expansion, alpha.V, and phylogenetics
model_simple4$DIC # Venditti's Prior
model_simple5$DIC # Combined Gelman and Venditti Prior

# EXPERIMENTAL MODEL
# Gelman's V prior and Chris Venditti's Priors for R,B, and G
neo$log.maxGG3 <- log(neo$maxGG3)
neo$log.maxGG3

prior5<-list(R=list(V=1,nu=0.02),B=list(mu=(rep(0,2)),V=gelman.prior(~log.maxGG3, data=neo, scale=1+pi^2/3)),
             G=list(G1=list(V=1,nu=1,alpha.mu=1000,alpha.V=1000)))
# nu varying degree of belief parameter
# G is the covariance matrix of the random effects
# R is the covariance matrix of the residuals
# B is fixed effects
# V is the covariance matrix of the fixed effects
# alpha.mu is code for parameter-expanded priors

model_simple5<-MCMCglmm(log(maxGG3)~log(G3M.est),random=~Species_binomial,
                        family="gaussian",ginverse=list(Species_binomial=inv.phylo$Ainv),prior=prior5,
                        data=neo,nitt=500000,burnin=1000,thin=500)

summary(model_simple5)

allChains <- as.mcmc(cbind(model_simple5$Sol,model_simple5$VCV))
plotTrace(allChains,axes=TRUE,las=1)

# Check the Bayesian model_simple5
maxG3.est_G3M.est + geom_abline(intercept = -5.4808, slope = 0.7934)
