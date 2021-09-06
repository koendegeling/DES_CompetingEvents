# Supplementary to the manuscript:
#
#   Comparing modeling approaches for discrete event simulations with competing risks based on censored individual
#   patient data: a simulation study and illustration in colorectal cancer
#
#   by K Degeling, MJ IJzerman, CGM Groothuis-Oudshoorn, MD Franken, M Koopman, MS Clements, H Koffijberg
# 
# This code illustrates how the event-specific distributions (ESD) and unimodal distribution and regression model
# (UDR) modeling approaches can be implemented to model competing events in discrete event simulations using 
# censored individual patient data (IPD) in R. It does so using a dummy dataset that includes hypothetical
# observations for 312 individuals who where subject to two competing risks: progression and death. The dataset,
# represented by the df_IPD data.frame that has been made available through the df_IPD.RData file together with
# this script, includes the following variables/columns:
# - id      unique identifier for each hypothetical individual
# - event   what event was observed: progression (Progression), death (Death), or no event (Censored)
# - time    time corresponding to the event, in case of progression or death, or time of censoring
#
# The aim of this script is to illustrate how the modeling approaches can be implemented, not to serve as a 
# survival analysis tutorial. Therefore, different parametric distribution types are not explored to model time-to-
# events, not are transformations of the time variable considered in the regression model used in the UDR approach.
# The scripts does demonstrate how the time-to-event distributions are fitted using the commonly used 'flexsurv'
# package, as well as the 'fitdistrplus' package, which is considered more userfriendly by some.
#
# The script is structured as follows. First, some administrative tasks are performed, which includes loading the
# data and required packages. Second, the data is analyzed according to the different approach. Third, simulation
# are performed using the outputs of the data analysis. Fourth and finally, the results of the simulation are
# compared against the original data for internal validation.
#
# Updates of and extensions to this script may be provided online: [BLINDED FOR PEER REVIEW]
#
# Koen Degeling, Friday 5 February 2021, R version 3.6.1




### 1. INITIALIZATION ----

# This first section of code performs some administrative tasks

# Uncomment to clear the environment
#rm(list = ls()); gc();


# Loading the dummy dataset
# - Update the path or your working directory if necessary, or load the file through the Files explorer
load(file = "df_IPD.RData");


# Uncomment to install packages if they are not already installed
#install.packages(pkgs = "flexsurv");
#install.packages(pkgs = "fitdistrplus");
#install.packages(pkgs = "dplyr");

# Loading packages
library(flexsurv);      # version 1.1.1, for the survival analyses
library(fitdistrplus);  # version 1.1-3, for the survival analyses
library(dplyr);         # version 1.0.4, for elegant data manipulation




### 2. DATA ANALYSIS ----

# This second sections performs the data analysis according to the two modeling approaches. The time-to-event
# analyses, or survival analyses if you prefer, are illustrated according to two commonly used packages: the
# flexsurv package and the fitdistrplus package, which have a quite substantially different syntax. Some background
# information on the syntax is provided, but please refer to the help files for further information.

# Before performing the data analyses for the two different modeling approaches, have a look at the dataset
head(df_IPD);


## 2.1 Event-specific distributions approach (ESD) ----

# To perform the data analysis according to the ESD approach, the first step is to define the censoring indicator 
# for both events, i.e. progression (c_progr) and death (c_death)
# - Note that 1 indicates that the event of interest occured and 0 indicates the other event or censoring
df_IPD <- df_IPD %>% mutate(
  c_progr = if_else(event %in% c("Death", "Censored"), 0, 1),
  c_death = if_else(event %in% c("Progression", "Censored"), 0, 1)
);


# Using this information, the time-to-event analysis can be performed using the flexsurvreg function of the 
# flexsurv package that requires a Surv() object as input, which can be specified using the time variable and 
# event-specific censoring indicator. Weibull distributions will be used throughout this illustration.
# - Note that both distributions can be conveniently fit independently
# - By wrapping the code in curved brackets, the result is printed to the console
(fit_flexsurv_ESD_progr <- flexsurvreg(formula = Surv(time = time, event = c_progr) ~ 1, 
                                       data    = df_IPD, 
                                       dist    = "weibull"));

(fit_flexsurv_ESD_death <- flexsurvreg(formula = Surv(time = time, event = c_death) ~ 1, 
                                       data    = df_IPD, 
                                       dist    = "weibull"));


# To fit the Weibull distribution using the fitdistcens function of the fitdistrplus package, some pre-processing
# of the data is required, because the function requires a data.frame with two columns: left and right, to define
# the type of censoring. In short, for right-censored data, column left contains the times for all individuals, 
# whereas the right column should be NA for censored individuals.
# - Note that individuals who experienced the competing event are also considered censored in the ESD approach
df_fitdistcens_progr <- data.frame(
  left  = df_IPD$time,
  right = ifelse(df_IPD$event %in% c("Death", "Censored"), NA, df_IPD$time)
);

df_fitdistcens_death <- data.frame(
  left  = df_IPD$time,
  right = ifelse(df_IPD$event %in% c("Progression", "Censored"), NA, df_IPD$time)
);

(fit_fitdist_ESD_progr <- fitdistcens(censdata = df_fitdistcens_progr, 
                                      distr    = "weibull"));

(fit_fitdist_ESD_death <- fitdistcens(censdata = df_fitdistcens_death, 
                                      distr    = "weibull"));


# As final step the estimated parameters will be extracted from the fit object so that they can easily be used to 
# perform simulations later on
# - Note that the flexsurv fit object includes coefficients that are not always on the real scale - in case of the
#   Weibull distribution both the shape and scale coefficients need to be exponentiated to obtain the parameters
(pars_flexsurv_ESD_progr <- exp(fit_flexsurv_ESD_progr$coefficients));
(pars_fitdist_ESD_progr  <- fit_fitdist_ESD_progr$estimate);

(pars_flexsurv_ESD_death <- exp(fit_flexsurv_ESD_death$coefficients));
(pars_fitdist_ESD_death  <- fit_fitdist_ESD_death$estimate);


## 2.2 Unimodal distribution and regression model (UDR) ----

# To perform the data analysis according to the UDR approach, the first step is to define the censoring indicator 
# for both events combined
# - Note that 1 indicates that any event occured and 0 indicates censoring
df_IPD <- df_IPD %>% mutate(
  c = if_else(event == "Censored", 0, 1)
);


# Using the flexsurvreg function of the flexsurv package
(fit_flexsurv_UDR <- flexsurvreg(formula = Surv(time = time, event = c) ~ 1, 
                                 data    = df_IPD, 
                                 dist    = "weibull"));


# Using the fitdistcens function of the fitdistrplus package
# - Now the right column only is NA if the individual was censored
df_fitdistcens <- data.frame(
  left  = df_IPD$time,
  right = ifelse(df_IPD$event == "Censored", NA, df_IPD$time)
);

(fit_fitdist_UDR <- fitdistcens(censdata = df_fitdistcens, 
                                distr    = "weibull"));


# Extracting the fitted parameters
(pars_flexsurv_UDR <- exp(fit_flexsurv_UDR$coefficients));
(pars_fitdist_UDR  <- fit_fitdist_UDR$estimate);


# The final step in the data analysis for the UDR approach is fitting the logistic regression model that predicts
# the type of event based on the time-to-event
# - Note that only individuals for whom an event was observed are included in the df_logreg data.frame
# - Note that the logistic regression model was defined to predict progression rather than death
df_logreg <- filter(df_IPD, event != "Censoring");

(fit_logreg_UDR <- glm(formula = (event == "Progression") ~ time, 
                       family  = binomial(link = "logit"), 
                       data    = df_logreg));




### 3. SIMULATION ----

# This third section illustrates how simulation are performed according the two modeling approaches. Although two 
# sets of parameters have been estimated for the time-to-event distributions, only one set will be used in this
# illustration, because the only difference would be to change the name of the objects that contain the parameters.

# Set the random number seed for reproducibility
set.seed(1);

# Define the number of individuals to be simulated
n_sim <- 10*10^3;


## 3.1 Event-specific distributions approach (ESD) ----

# The first step in simulating according to the ESD approach is to draw a time to each competing event from the
# corresponding time-to-event distributions
df_sim_ESD <- data.frame(
  time_progr = rweibull(n     = n_sim, 
                        shape = pars_fitdist_ESD_progr["shape"], 
                        scale = pars_fitdist_ESD_progr["scale"]),
  time_death = rweibull(n     = n_sim, 
                        shape = pars_fitdist_ESD_death["shape"], 
                        scale = pars_fitdist_ESD_death["scale"])
);


# The second and final step is to select the event that is the first to occur
df_sim_ESD <- df_sim_ESD %>% mutate(
  time  = pmin(time_progr, time_death),
  event = if_else(time_progr < time_death, "Progression", "Death")
);


## 3.2 Unimodal distribution and regression model (UDR) ----

# The first step in simulating for to the UDR approach is to sample from the combined time-to-event distributions
df_sim_UDR <- data.frame(
  time = rweibull(n     = n_sim, 
                  shape = pars_fitdist_UDR["shape"], 
                  scale = pars_fitdist_UDR["scale"])
);


# The second and final step is to sample the event corresponding to the sampled time
# - Note that the logistic regression was defined to predict progression
# - Note that probabilities are predicted using the predict function, by setting the type to "response", and 
#   subsequently the event is sampled using a random number from a Uniform distribution
df_sim_UDR <- df_sim_UDR %>% mutate(
  p_progr = predict(object = fit_logreg_UDR, newdata = data.frame(time = time), type = "response"),
  event   = if_else(runif(n = n()) < p_progr, "Progression", "Death")
);




### 4. VALIDATION ----

# This fourth and final section illustrates the results of this demonstration as a form of validation. Note that
# this is just an illustration using Weibull distributions for all time-to-event distributions and without 
# exploring transformations for the time variable in the logistic regression model used in the UDR approach.

# Inspecting the observed vs. simulated event incidences
# - Note the substantial differences between the two approaches: 72% progression for ESD vs 61% according to UDR
prop.table(table(df_IPD$event));      # Data
prop.table(table(df_sim_ESD$event));  # ESD approach
prop.table(table(df_sim_UDR$event));  # UDR approach


# To inspect the time-to-event curves, Kaplan-Meier curves are fitted and subsequently plotted
# - Again, note the differences between the methods: ESD better represents progression and UDR (probably) death
km_data_progr <- survfit(formula = Surv(time = time, event = (event == "Progression")) ~ 1, 
                         data    = df_IPD);
km_data_death <- survfit(formula = Surv(time = time, event = (event == "Death")) ~ 1, 
                         data    = df_IPD);

km_ESD_progr <- survfit(formula = Surv(time = time, event = (event == "Progression")) ~ 1, 
                        data    = df_sim_ESD);
km_ESD_death <- survfit(formula = Surv(time = time, event = (event == "Death")) ~ 1, 
                        data    = df_sim_ESD);

km_UDR_progr <- survfit(formula = Surv(time = time, event = (event == "Progression")) ~ 1, 
                        data    = df_sim_UDR);
km_UDR_death <- survfit(formula = Surv(time = time, event = (event == "Death")) ~ 1, 
                        data    = df_sim_UDR);

plot(km_data_progr, main = "Time to PROGRESSION", las = 1, xlab = "Time", ylab = "Proportion event-free");
lines(x = km_ESD_progr$time, y = km_ESD_progr$surv, col = "blue");
lines(x = km_UDR_progr$time, y = km_UDR_progr$surv, col = "purple");
legend("topright", legend = c("Data", "ESD", "UDR"), col = c("black", "blue", "purple"), lty = 1, bty = "n");

plot(km_data_death, main = "Time to DEATH", las = 1, xlab = "Time", ylab = "Proportion event-free");
lines(x = km_ESD_death$time, y = km_ESD_death$surv, col = "blue");
lines(x = km_UDR_death$time, y = km_UDR_death$surv, col = "purple");
legend("topright", legend = c("Data", "ESD", "UDR"), col = c("black", "blue", "purple"), lty = 1, bty = "n");



