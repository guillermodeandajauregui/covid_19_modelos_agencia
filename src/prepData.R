prepData <-
function(data){
  # The function `plot.model` provides automated visualization of model simulations
  
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
      colnames(i) <-paste0("sim", 1:25)
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