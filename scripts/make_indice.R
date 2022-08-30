# Objetivo Crear un diccionario de variable


# abreviaturas


# nombres extensos
delitos_abrev <- c("hurto", "robo", "asalto", "amen", "abus", "crim", "viol")
delitos_extense <- c("Minor property offenses","Robberies","Burglaries",
                     "Drunkennes, damages and disorders",
                     "Family violence and aggression",
                     "Injuries, drugs and weapons","High violence and murder")

dic_del <- data.frame(delitos_abrev, delitos_extense)
# View(dic_del)

variables_abrev <- c("pobdens", "tamhog", "educ", "vuln", "mobil", "flux",
                     "ofcom")
variables_extense <- c("Pop.density", "Household.size","Education.av",
                       "Neigh.disadvantage","Neigh.mobility",
                       "flux", "Facilities")

dic_variables <- data.frame(variables_abrev, variables_extense)

names(dic_del) <- c("abrev", "Variables")
names(dic_variables) <- c("abrev", "Variables")
indice <- rbind(dic_del, dic_variables)

library(openxlsx)
write.csv(indice, "docs/indice_var.csv", row.names = F)
