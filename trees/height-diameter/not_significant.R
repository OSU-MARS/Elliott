

#get the models that are not significant from base model

temp<-heightDiameterResults %>% group_by(responseVariable,species) %>% filter(significant=="FALSE") %>% reframe(unique(name))
write.csv(temp,"not_significant.csv")
