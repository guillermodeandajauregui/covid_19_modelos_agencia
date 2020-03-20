source("src/functions.R")

#run model
probando <- 
evaluate.model(distanciamiento_social = 100, 
               aislamiento_enfermos = T, 
               iniciales = 1)

#extract results
mis_resultados <- extract.model.data(probando)

#plot HOSPITALIZED mean vs time

mis_resultados$my_sim_means %>% 
  filter(variable == "H") %>% 
  group_by(cum.time.day) %>% 
  summarise(hospitalized = sum(value)) %>% 
  ggplot(mapping = aes(x = cum.time.day, y = hospitalized)) +
  geom_line()
