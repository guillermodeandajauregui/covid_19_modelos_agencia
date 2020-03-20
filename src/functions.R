################################################################################
#functions for SIR model. 
################################################################################

library(tidyverse)
library(vroom)

################################################################################
#distribution functions
################################################################################
gamma <- function(z = 12, b=0.143, a0=1/1.5, t){
  # piecewise function
  # default parameters z = 12, b=1/7, a0=1/1.5
  #    z: time at start of intervention (notionally March 12)
  #    b: intercept (positive)
  #    a0: post intervention isolation ratae
  #    t: time in the model
  
  gamma <- ifelse(t<=z, gamma <- b, gamma <- a0)
  return(gamma)
}

eta <- function(t, w=12) ifelse(t<=w,1/3,1/3)

q <- function(t, w=12, q0=1, q1=1) ifelse(t<=w,q0,q1)

beta <- function(t, w=12, beta0=0.6584, beta.factor=2) ifelse(t<=w,beta0,beta0/beta.factor)
################################################################################
#model step 
################################################################################

onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  
  S <- x[2]                            #local variable for susceptibles
  
  E1 <- x[3]                           #exposed classes
  E2 <- x[4]
  E3 <- x[5]
  E4 <- x[6]
  E5 <- x[7]
  E6 <- x[8]
  
  I1 <- x[9]                           #detected infectious classes
  I2 <- x[10]
  I3 <- x[11]
  I4 <- x[12]
  
  Iu1 <- x[13]                           #undetected infectious classes
  Iu2 <- x[14]
  Iu3 <- x[15]
  Iu4 <- x[16]
  
  I.detected <-I1+I2+I3+I4
  I.undetected <- Iu1+Iu2+Iu3+Iu4
  I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infectious
  
  H <- x[17]                           #local variable for hospitalized
  Ru <- x[18]                          #local variable for undetected recovereds
  
  C <- x[19]    # local variable for notifications
  
  N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
  t <- x[20]                          #get current time
  
  with(                                #use with to simplify code
    as.list(params), 
    {
      gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
      sigmai <- 6*sigma  # multiplier 6 for pseudo stages
      etat <- eta(t,w)     # case notification rate
      betat <- beta(t,w)   # time dependent transmissibility, presymptomatic=1 causes this transmissibility to apply to late stage latent cases as well
      
      rates <- as.numeric(c(betat*I.detected/N+betat*c*I.undetected/N+presymptomatic*betat*c*E6/N,                           # movements out of S
                            sigmai, sigmai, sigmai, sigmai, sigmai, sigmai,   # movements out of E
                            gammai, gammai, gammai, gammai,                   # movements out of I (detected)
                            b, b, b, b,                                       # movements out of I (undetected)
                            etat))
      
      states0 <- x[2:(length(x)-1)]
      
      # transition probabilities
      
      p <- matrix(0, nrow=length(rates),ncol=length(states0))                                     # matrix to hold transitions probs
      
      p[1,]  <- c(exp(-rates[1]*dt),1-exp(-rates[1]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # S-> E
      
      p[2,]  <- c(0, exp(-rates[2]*dt), 1-exp(-rates[2]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[3,]  <- c(0, 0, exp(-rates[3]*dt), 1-exp(-rates[3]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[4,]  <- c(0, 0, 0, exp(-rates[4]*dt), 1-exp(-rates[4]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[5,]  <- c(0, 0, 0, 0, exp(-rates[5]*dt), 1-exp(-rates[5]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[6,]  <- c(0, 0, 0, 0, 0, exp(-rates[6]*dt), 1-exp(-rates[6]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[7,]  <- c(0, 0, 0, 0, 0, 0, exp(-rates[7]*dt), (1-exp(-rates[7]*dt))*q(w), 0, 0, 0, (1-exp(-rates[7]*dt))*(1-q(w)), 0, 0, 0, 0, 0, 0)    # Transitions out of E
      
      p[8,]  <- c(0, 0, 0, 0, 0, 0, 0, exp(-rates[8]*dt), 1-exp(-rates[8]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I (detected)
      p[9,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[9]*dt), 1-exp(-rates[9]*dt), 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I
      p[10,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[10]*dt), 1-exp(-rates[10]*dt), 0, 0, 0, 0, 0, 0, 0)  # Transitions out of I
      p[11,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[11]*dt), 0, 0, 0, 0, 1-exp(-rates[11]*dt), 0, 0)  # Transitions out of I -> H
      
      p[12,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[12]*dt), 1-exp(-rates[12]*dt), 0, 0, 0, 0, 0)    # Transitions out of I (undetected)
      p[13,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[13]*dt), 1-exp(-rates[13]*dt), 0, 0, 0, 0)    # Transitions out of I
      p[14,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[14]*dt), 1-exp(-rates[14]*dt), 0, 0, 0)  # Transitions out of I
      p[15,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[15]*dt), 0, 1-exp(-rates[15]*dt), 0)  # Transitions out of I -> R_u
      
      
      p[16,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[16]*dt), 0, 1-exp(-rates[16]*dt))  # Transitions R_d -> to C (notification)
      
      # update states
      
      states1 <- matrix(0, nrow=length(rates),ncol=length(states0))                                # matrix to hold future states
      
      for(i in 1:length(rates)){
        states1[i,] <- t(rmultinom(1, states0[i], p[i,])) 
      }
      
      states1 <- colSums(states1)
      states1[17] <- states1[17]+Ru  #add formerly Recovered undetected cases
      states1[18] <- states1[18]+C  #add formerly notified cases
      
      return(x <- c(dt, states1, tail(x,1)+dt))
    }
  )
}

################################################################################
#run model 
################################################################################

model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,length(x)))         #set up array to store results
  colnames(output) <- c("time",
                        "S",
                        "E1", "E2", "E3", "E4", "E5", "E6",
                        "I1", "I2", "I3", "I4", 
                        "Iu1", "Iu2", "Iu3", "Iu4",
                        "H", 
                        "Ru", 
                        "C", 
                        "cum.time") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- as.numeric(onestep(x,params))
  }
  output                                    #return output
}

################################################################################
#evaluate model 
################################################################################

evaluate.model <- function(distanciamiento_social=100,
                           aislamiento_enfermos=TRUE,
                           iniciales = 1,
                           nsims=2, 
                           nstep=NULL, 
                           start=as.Date("2020-03-01"),
                           today=Sys.Date()
){
  
  #adapting parameters for backward compatibility with original Imperial College data
  #DO NOT CHANGE THESE PARAMETERS
  

  ########################
  #model parametrization
  #DO NOT MODIFY 
  ########################
  
  params=list(beta0=0.6584, 
              sigma=1/6.4, 
              z = ifelse(aislamiento_enfermos==TRUE, 12, 1200),
              b=0.143, 
              a0=1/1.5, 
              w = 1200*(distanciamiento_social/100),
              c=1, 
              presymptomatic=1, 
              dt=0.05)
  
  ########################
  #initial compartments
  #DO NOT MODIFY 
  ########################
  init = list(S=9000000, 
              E1=0, 
              E2=0, 
              E3=0, 
              E4=0, 
              E5=6, 
              E6=0,
              I1 = iniciales, 
              I2= 0, 
              I3=0, 
              I4=0, 
              Iu1=0, 
              Iu2=0, 
              Iu3=0, 
              Iu4=0,
              H=0, 
              Ru=0, 
              C=0)
  ##################
  if(is.null(nstep)) nstep <- (as.numeric(today-start)+1+28)/params$dt #run simulation from start to current time plus four weeks
  
  xstart <- c(time=0, unlist(init), cum.time = 0) #initial conditions
  
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  
  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    
  }
  
  return(data)
}

################################################################################
#extract model data 
################################################################################
extract.model.data <-
  function(data){
    #takes a model's output, and returns the simulations, and the means, in tidy format
    #CHECK THE COMMENTS TO UNDERSTAND THE VOLUMES
    
    #E  = Exposed group
    #I  = Infected group
    #Iu = Infected - unreported group
    #H  = HOSPITALIZED group  
    
    # process data
    nsims <- length(data)
    
    for(i in 1:nsims) data[[i]]$I <- data[[i]]$I1 + data[[i]]$I2 + data[[i]]$I3 +
        data[[i]]$I4
    for(i in 1:nsims) data[[i]]$Iu <- data[[i]]$Iu1 + data[[i]]$Iu2 + data[[i]]$Iu3 +
        data[[i]]$Iu4
    for(i in 1:nsims) data[[i]]$E <- data[[i]]$E1 + data[[i]]$E2 + data[[i]]$E3 +
        data[[i]]$E4 + data[[i]]$E5 + data[[i]]$E6
    
    max.time<-data[[1]]$cum.time[max(which(data[[1]]$I>0))] #maximum time in first simulation
    max.y<-max(data[[1]]$C)       #find max total confirmed cases for plotting range
    
    #cumulative time
    cum.time.df <- data.frame(cum.time = data[[1]]$cum.time)
    
    # calculate means
    m1 <- m2 <- m3 <- m4 <- m5 <- matrix(nrow=length(data[[1]]$I), ncol=nsims)
    for(i in 1:nsims){
      m1[,i] <- data[[i]]$E
      m2[,i] <- data[[i]]$I+data[[i]]$Iu
      m3[,i] <- data[[i]]$Iu
      m4[,i] <- data[[i]]$H
      m5[,i] <- data[[i]]$C
    }
    E.mean <- rowMeans(m1)
    I.mean <- rowMeans(m2)
    Iu.mean <- rowMeans(m3)
    H.mean <- rowMeans(m4)
    C.mean <- rowMeans(m5)
    
    my_sim_results <- 
      list(E = m1,
           I = m2,
           Iu = m3,
           H = m4,
           C = m5)
    
    my_sim_results <- 
      lapply(my_sim_results, as.data.frame) %>% 
      lapply(function(i){
        colnames(i) <-paste0("sim", 1:length(colnames(i)))
        i = cbind(cum.time.df, i)
        i <- 
          i %>% 
          mutate(cum.time.day = ceiling(cum.time)) %>% 
          mutate(cum.time.day = ifelse(cum.time.day==0, 1, cum.time.day))
        return(i)
      })
    
    my_sim_results <-
      my_sim_results %>% 
      lapply(FUN = function(i){
        pivot_longer(data = i,
                     cols = -c(cum.time, cum.time.day), 
                     names_to = "sim")
      })
    
    my_sim_results <-
      my_sim_results %>% 
      bind_rows(.id = "variable")
    
    my_sim_results <-
      my_sim_results %>% 
      select(cum.time, cum.time.day, sim, variable, value)
    
    my_sim_means <- 
      data.frame(
        cum.time = cum.time.df$cum.time,
        E = E.mean,
        I = I.mean,
        Iu = Iu.mean,
        H  = H.mean,
        C  = C.mean
      )
    
    my_sim_means <-
      my_sim_means %>% 
      mutate(cum.time.day = ceiling(cum.time)) %>% 
      mutate(cum.time.day = ifelse(cum.time.day==0, 1, cum.time.day))
    
    my_sim_means <-
      my_sim_means %>% 
      pivot_longer(cols = -c(cum.time, cum.time.day), 
                   names_to = "variable")
    
    my_sim_means <-
      my_sim_means %>% 
      mutate(sim = "mean") %>% 
      select(cum.time, cum.time.day, sim, variable, value)
    
    my_return <-
      list(my_sim_results = my_sim_results,
           my_sim_means = my_sim_means)
    
    return(my_return)
  }