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

# Create functions for cloning

#------------------------------------------------------------------------------#
# Function for creating a clone for the recommended protocol
#------------------------------------------------------------------------------#

func.create.clone.rec <- function (d) {
  
  p <- d[which(d$type.vaxx.1=="PFR"),]
  m <- d[which(d$type.vaxx.1=="MOD"),]
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Pfizer)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # PFIZER recommended
  p1 <- p[which(p$dose.sch=="recommended"),] 
  p1$status <- ifelse(p1$post.inf==1, 1, 0)
  p1$time <- ifelse(p1$status==1 , 
                    p1$time.inf, 
                    p1$time.censor)
  p.rec <- p1
  
  # PFIZER early
  p1 <- p[which(p$dose.sch=="early"),] 
  p1$status <- ifelse(p1$post.inf==1 & 
                        (p1$time.inf < p1$days.btw.doses.12), 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$days.btw.doses.12)
  p.ear <- p1
  
  # PFIZER allowable or late
  p1 <- p[which(p$dose.sch=="allowable" | p$dose.sch=="late"),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < 26), 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    25)
  p.rest <- p1
  
  # PFIZER missing
  p1 <- p[which(is.na(p$dose.sch)),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < p1$days.btw.dose1.mar16) &
                        (p1$time.inf < 26),  1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    ifelse(p1$days.btw.dose1.mar16<26, p1$days.btw.dose1.mar16, 25))
  p.na <- p1
  p.c.rec <- rbind(p.rec, p.ear, p.rest, p.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Moderna)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # MODERNA recommended
  m1 <- m[which(m$dose.sch=="recommended"),] 
  m1$status <- ifelse(m1$post.inf==1, 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$time.censor)
  m.rec <- m1
  
  # MODERNA early
  m1 <- m[which(m$dose.sch=="early"),] 
  m1$status <- ifelse(m1$post.inf==1 & 
                        (m1$time.inf < m1$days.btw.doses.12), 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$days.btw.doses.12)
  m.ear <- m1
  
  # MODERNA allowable or late
  m1 <- m[which(m$dose.sch=="allowable" | m$dose.sch=="late"),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < 33), 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    32)
  m.rest <- m1
  
  # MODERNA missing
  m1 <- m[which(is.na(m$dose.sch)),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < m1$days.btw.dose1.mar16) &
                        (m1$time.inf < 33),  1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    ifelse(m1$days.btw.dose1.mar16<33, m1$days.btw.dose1.mar16, 32))
  
  m.na <- m1
  m.c.rec <- rbind(m.rec, m.ear, m.rest, m.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Final cloned datasets
  #-----*-----*-----*-----*-----*-----*-----#
  
  # Merge
  c.rec  <- rbind(p.c.rec,  m.c.rec)
  
  # Create a status variable for censoring
  c.rec$status.censor  <- ifelse(c.rec$status ==1, 0, 1)
  
  return(c.rec)
}

#------------------------------------------------------------------------------#
# Function for creating a clone for the allowable protocol
#------------------------------------------------------------------------------#

func.create.clone.allw <- function (d) {
  
  p <- d[which(d$type.vaxx.1=="PFR"),]
  m <- d[which(d$type.vaxx.1=="MOD"),]
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Pfizer)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # PFIZER recommended or early
  p1 <- p[which(p$dose.sch=="recommended" | p$dose.sch=="early"),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < p1$days.btw.doses.12), 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$days.btw.doses.12)
  p.rest <- p1
  
  # PFIZER allowable
  p1 <- p[which(p$dose.sch=="allowable"),] 
  p1$status <- ifelse(p1$post.inf==1, 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$time.censor)
  p.allw <- p1
  
  # PFIZER late
  p1 <- p[which(p$dose.sch=="late"),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < 43), 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    42)
  p.late <- p1
  
  # PFIZER missing
  p1 <- p[which(is.na(p$dose.sch)),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < p1$days.btw.dose1.mar16) &
                        (p1$time.inf < 43),  1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    ifelse(p1$days.btw.dose1.mar16<43, p1$days.btw.dose1.mar16, 42))
  p.na <- p1
  p.c.allw <- rbind(p.allw, p.late, p.rest, p.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Moderna)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # MODERNA recommended or early
  m1 <- m[which(m$dose.sch=="recommended" | m$dose.sch=="early"),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < m1$days.btw.doses.12), 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$days.btw.doses.12)
  m.rest <- m1

  # MODERNA allowable
  m1 <- m[which(m$dose.sch=="allowable"),] 
  m1$status <- ifelse(m1$post.inf==1, 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$time.censor)
  m.allw <- m1
  
  # MODERNA late
  m1 <- m[which(m$dose.sch=="late"),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < 50), 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    49)
  m.late <- m1
  
  # MODERNA missing
  m1 <- m[which(is.na(m$dose.sch)),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < m1$days.btw.dose1.mar16) &
                        (m1$time.inf < 50),  1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    ifelse(m1$days.btw.dose1.mar16<50, m1$days.btw.dose1.mar16, 49))
  m.na <- m1
  m.c.allw <- rbind(m.allw, m.late, m.rest, m.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Final cloned datasets
  #-----*-----*-----*-----*-----*-----*-----#
  
  # Merge
  c.allw <- rbind(p.c.allw, m.c.allw)
  
  # Create a status variable for censoring
  c.allw$status.censor <- ifelse(c.allw$status==1, 0, 1)
  
  return(c.allw)
}


#------------------------------------------------------------------------------#
# Function for creating a clone for the late protocol
#------------------------------------------------------------------------------#

func.create.clone.late <- function (d) {
  
  p <- d[which(d$type.vaxx.1=="PFR"),]
  m <- d[which(d$type.vaxx.1=="MOD"),]
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Pfizer)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # PFIZER recommended or early or allowable
  p1 <- p[which(p$dose.sch=="recommended" | p$dose.sch=="early" | p$dose.sch=="allowable"),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < p1$days.btw.doses.12), 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$days.btw.doses.12)
  p.rest <- p1

  # PFIZER late
  p1 <- p[which(p$dose.sch=="late"),] 
  p1$status <- ifelse(p1$post.inf==1, 1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$time.censor)
  p.late <- p1
  
  # PFIZER missing
  p1 <- p[which(is.na(p$dose.sch)),] 
  p1$status <- ifelse(p1$post.inf==1 & (p1$time.inf < p1$days.btw.dose1.mar16),  1, 0)
  p1$time <- ifelse(p1$status==1, 
                    p1$time.inf, 
                    p1$days.btw.dose1.mar16)
  p.na <- p1
  p.c.late <- rbind( p.rest, p.late, p.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Create clones (Moderna)
  #-----*-----*-----*-----*-----*-----*-----#
  
  # MODERNA recommended or early or allowable
  m1 <- m[which(m$dose.sch=="recommended" | m$dose.sch=="early" | m$dose.sch=="allowable"),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < m1$days.btw.doses.12), 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$days.btw.doses.12)
  m.rest <- m1
  
  # MODERNA late
  m1 <- m[which(m$dose.sch=="late"),] 
  m1$status <- ifelse(m1$post.inf==1, 1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$time.censor)
  m.late <- m1
  
  # MODERNA missing
  m1 <- m[which(is.na(m$dose.sch)),] 
  m1$status <- ifelse(m1$post.inf==1 & (m1$time.inf < m1$days.btw.dose1.mar16),  1, 0)
  m1$time <- ifelse(m1$status==1, 
                    m1$time.inf, 
                    m1$days.btw.dose1.mar16)
  m.na <- m1
  m.c.late <- rbind( m.rest, m.late, m.na)
  
  #-----*-----*-----*-----*-----*-----*-----#
  # Final cloned datasets
  #-----*-----*-----*-----*-----*-----*-----#
  
  # Merge
  c.late <- rbind(p.c.late, m.c.late)
  
  # Create a status variable for censoring
  c.late$status.censor <- ifelse(c.late$status==1, 0, 1)
  
  return(c.late)
}