################################################################################
#                                                                              #
#  MANUSCRIPT: Comparative Effectiveness of Alternative Intervals between      #
#              First and Second Doses of the mRNA COVID-19 Vaccines:           #
#              a Target Trial Emulation Approach                               #
#                                                                              #
#    CODED BY: Kayoko Shioda, PhD, DVM, MPH                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Description of this script
#------------------------------------------------------------------------------#

# Main analysis of the paper.

#------------------------------------------------------------------------------#
# Set up HPC tasks
#------------------------------------------------------------------------------#

# Set up a task ID
bt <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if(is.na(bt)) bt<-1 
print(paste("My task is", bt))

#------------------------------------------------------------------------------#
# Set up input parameters
#------------------------------------------------------------------------------#

manufacturer     <- "MOD"   # MOD if Moderna. PFR if Pfizer.
analysis.title   <- "main"  

#------------------------------------------------------------------------------#
# Load functions
#------------------------------------------------------------------------------#

# Load functions
library(tidyr)
library(feather)
library(survival)
library(dplyr)
library(data.table)
library("scales")

#------------------------------------------------------------------------------#
# Set up and load the dataset
#------------------------------------------------------------------------------#

# Set up the folder name 
folder.name <- paste0("/rprojectnb/shiolab/TTE/Results/", analysis.title, "/", manufacturer) 

# Load functions for cloning
source("func.create.clones.R")

# Load the original dataset
orig.dt <- read_feather(path="./Data/d.wide.final.feather") 

#------------------------------------------------------------------------------#
# Data processing
#------------------------------------------------------------------------------#

# Remove those who had the 2nd dose >180 days after the 1st dose
d0 <- orig.dt[which(is.na(orig.dt$days.btw.doses.12) |  
                      (orig.dt$days.btw.doses.12>=4 & orig.dt$days.btw.doses.12<=180)),] 

# Keep the relevant columns and drop the rest
d0 <- subset(d0, select=c(RECIP_ID,
                          age, sex, race, ethnicity, recip_district, 
                          type.vaxx.1, dose.sch, dose.sch.ver2, 
                          date.pos.aft.vaxx,
                          days.btw.doses.12, 
                          days.btw.dose1.mar16,
                          days.btw.dose1.inf, 
                          prior.inf, post.inf, 
                          inf.moyr, dose1.moyr,
                          time.inf, time.censor))

# Select a subset for each manufacturer to be analyzed
d0 <- d0[which(d0$type.vaxx.1==manufacturer),]

#------------------------------------------------------------------------------#
# Set up for the clone censor weight analysis
#------------------------------------------------------------------------------#

# Create a formula for the Cox PH model for censor weights
fm <- as.formula(paste("Surv(time, status.censor) ~ age + sex + race + ethnicity + prior.inf + recip_district + dose1.moyr"))

# Create a vector with the name of the protocols to be evaluated
name_protocols <- c("recommended", "allowable", "late")

# Set the end of the evaluation period
end.eval <- 180 # days after the 1st dose

# Create data frames to store TTE results
res.rec  <- data.frame(time=seq(1,end.eval,1), cumrisk=rep(NA,end.eval))
res.allw <- data.frame(time=seq(1,end.eval,1), cumrisk=rep(NA,end.eval))
res.late <- data.frame(time=seq(1,end.eval,1), cumrisk=rep(NA,end.eval))

#------------------------------------------------------------------------------#
# Run the clone censor weight analysis
#------------------------------------------------------------------------------#

set.seed(123456 + bt)

print(paste("Bootstrap No.",bt))

# Bootstrap
if (bt==1) {
  d <- d0
} else if (bt>1) {
  d <- d0 %>% # Select a clone you want to work with
    slice_sample(prop = 1, replace = TRUE) %>%
    group_by(RECIP_ID) %>%
    mutate(replicate = row_number()) %>%
    mutate(RECIP_ID = paste0(RECIP_ID, "_", replicate))
}

# No one was uncensored in March 2022, so exclude that month from analysis
d <- d[which(d$dose1.moyr!="2022-03"),]

# Set the number of total individuals in the dataset
n_subjects <- nrow(d)

for (prt in 1:length(name_protocols)) { # For each protocol...
  
  print(name_protocols[prt]) # Which protocol are we running?
  
  # Create a protocol-specific clone using the function that was built in a 
  # separate R script (func.create.clones.R)
  if (name_protocols[prt] == "recommended") {
    clone <- func.create.clone.rec(d)
  } else if (name_protocols[prt] == "allowable") {
    clone <- func.create.clone.allw(d)
  } else if (name_protocols[prt] == "late") {
    clone <- func.create.clone.late(d)
  }
  
  # Separately for each set of clones, fit a Cox Model for censoring
  # (i.e., status=1 if censored) conditional on covariates.
  coxph.censor <- coxph(fm, data = clone)
  
  # Subset the population to just those who had events 
  # (i.e., COVID infection in this case)
  cases <- clone[which(clone$status==1), ]
  
  # Predict the probability of remaining uncensored at each person’s event time
  # (each day of COVID infection), using the survival package with broom.
  cases <- broom::augment(coxph.censor, 
                          newdata=cases, 
                          type.predict = "expected") %>%
    dplyr::mutate(prob = exp(-.fitted))
  
  # Order the subset data by the event time 
  # (i.e., number of days between the first dose and COVID infection)
  cases <- cases[order(cases$time.inf),]
  
  # Compute the cumulative sum of 1/predictions at each of the times, which is
  # our cumulative incidence curve
  cases$wt   <- 1/cases$prob
  cases$risk <- cumsum(cases$wt) / n_subjects
  
  # Extract the maximum risk (i.e., cumulative risk) on each event day
  res <- cases[which(cases$time<=180),] %>%
    group_by(time) %>%
    summarise(cumrisk = max(risk, na.rm=TRUE))
  
  # Check if any days are missing in the event time (Day 1-180 after the 1st
  # dose administration), and if so, add it.
  if (length(unique(res$time)) != 180) {
    comp.days <- data.frame(time=seq(1, 180, 1))
    res <- merge(comp.days, res, by = "time", all.x = T)
  }
  
  # If any event days were missing, the cumulative risk should be the same as 
  # that on the previous day.
  if (is.na(res$cumrisk[1])) {res$cumrisk[1] <- 0}
  for (k in 2:nrow(res)) {
    if (is.na(res$cumrisk[k])) {
      res$cumrisk[k] <- res$cumrisk[k-1]
    }
  }
  
  # Save the outputs
  if (name_protocols[prt] == "recommended") {
    res.rec  <- res
  } else if (name_protocols[prt] == "allowable") {
    res.allw <- res
  } else if (name_protocols[prt] == "late") {
    res.late <- res
  }
  
  rm(clone, cases, res) 
}

# Save the results
save(res.rec,  file=paste0(folder.name, "/res.rec.",  bt, ".rda"))
save(res.allw, file=paste0(folder.name, "/res.allw.", bt, ".rda"))
save(res.late, file=paste0(folder.name, "/res.late.", bt, ".rda"))
