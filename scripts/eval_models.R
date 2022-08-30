# Directorio de Trabajo ---------------------------------------------------
setwd("~/OneDrive - Universidad Adolfo Ibanez/FONDECYT/Crime_Cellphones")

# Librerías y Recursos ----------------------------------------------------
source("scripts/recursos.R")


# Funciones ---------------------------------------------------------------
source("scripts/funciones.R")



# List of crimes/delitos and cities/ciudades
crimes=c("Minor property offenses","Robberies","Burglaries","Drunkennes, damages and disorders",
         "Family violence and aggression","Injuries, drugs and weapons","High violence and murder")
delitos = c("hurto", "robo", "asalto", "amen", "abus", "crim", "viol")
# delitos = c("hurto", "robo", "asalto", "amen")
# crimes=c("Minor property offenses","Robberies","Burglaries","Drunkennes, damages and disorders")
cities <- c("Coquimbo-Serena")
ciudades <- c("urb_4_18")


### Test and compilation functions

# Coefficients standardization function
betam <-  function (model)
{
  if (class(model)[1] == "Sarlm")
  {
    b <- summary(model)$coefficients[-1]
    sx <- apply(data.frame(model$X)[2:ncol(model$X)],2,FUN=sd)
    sy <- sd(model$y)
  }
  if (class(model)[1] == "glmerMod")
  {
    b <- summary(model)$coefficients[-1,1]
    sx <- apply(model@frame[2:(ncol(model@frame)-1)],2,FUN=sd)
    sy <- sd(unlist(model@frame[1]))
  }
  beta <- b * sx/sy
  return(beta)
}

# Spatial prediction error (SPE) function
ks.spat <-  function (cityp, nb, model, delito, val_ks=F)
{ 
  cityp$rows <- rownames(cityp@data) <- rownames(cityp@coords)
  dato <- cityp[,c("rows",delito)]
  
  # model predictions
  if (class(model)[1] == "Sarlm")
  {
    mod <- "sar"
    ypred <- data.frame(ypred=as.numeric(predict(model, listw = nb, pred.type = "TS"))**3)
    rownames(ypred) <- cityp$rows
  }
  if (class(model)[1] == "glmerMod")
  {
    mod <- "mlb"
    ypred <- data.frame(ypred=exp(predict(model)))
    rownames(ypred) <- cityp$rows
  }
  
  # classification error, based on KS method
  set.seed(1)
  deciles <- 10
  dato$ydcl <- cut(dato@data[,2], c(unique(getJenksBreaks(dato@data[,2], deciles)), Inf), labels = F, include.lowest = T, right = F)
  ypred$ypcl <- cut(ypred[,1], c(unique(getJenksBreaks(ypred[,1], deciles)), Inf), labels = F, include.lowest = T, right = F)
  ypred$rows <- rownames(ypred)
  dato <- sp::merge(dato,ypred,by="rows")
  dato$ycfit <- dato$ypcl - dato$ydcl
  fits <- dato@data[,c(1,4,3,5,6)]
  value_ks = nrow(fits[fits$ycfit==0,])/nrow(fits)
  names(fits) <- gsub("y",paste0(mod,delito),names(fits),fixed=T)
  out <- if(val_ks == T){
    tab_ks <- data.frame(
      Vars = "KS.spat",
             valor = value_ks)
    colnames(tab_ks) = c("Vars", delito)
    rownames(tab_ks) = "KS.spat"
    tab_ks
  }else{
    fits
  }
  return(out)
}

# Results compilation function
tabula.coef <- function(cityp, nb, model, delito)
{
  # data extraction
  if (class(model)[1] == "Sarlm")  {
    tab <- data.frame(summary(model)$Coef)
  }
  
  if (class(model)[1] == "glmerMod")  { 
    tab <- data.frame(summary(model)$coefficients)
    row.names(tab)[2] <- "lag_crime"
  }
  tab[-1,1] <- betam(model)
  tab$significancia <- ifelse(tab[,4]<0.01,"***",ifelse(tab[,4]<0.05,"**",ifelse(tab[,4]<0.1,"*","")))
  tab$var <- row.names(tab)
  tab$resumen <- paste0(round(tab$Estimate,2),tab$significancia)
  tab <- tab[,c("var","resumen")]
  colnames(tab) <- c("Vars", delito)
  
  # model fit tabulation
  if (class(model)[1] == "Sarlm")
  {
    errores <- data.frame(AIC(model))
  }
  if (class(model)[1] == "glmerMod")
  { 
    errores <- data.frame(summary(model)$AICtab[1])
  }
  errores$var <- "AIC"
  colnames(errores) = c(delito, "Vars")
  errores <- errores[,c("Vars", delito)]
  
  # spatial prediction error
  SKS <- ks.spat(cityp, nb, model, delito,val_ks = T)
  
  # compile
  tab <- rbind(tab, errores, SKS)
  names(tab)[2] <- delito
  return(tab)
}

### Parallell compilation of statistics and regression models

# Parameters
iter <- length(ciudades)

# Models compilation
foreach(obs=1:iter, .packages=c("spatialreg","MASS","BAMMtools","lme4","spdep")) %do% {
  
  city <- ciudades[obs] # city=ciudades[1]
  cityp <- readRDS(paste0("data/",city,".rds"))
  nb <- nb <- spatial_weights(city_df = cityp, nvec= 12)
  
  # multilevel lagged binomial
  regressions <- data.frame(Vars=NA)
  for (delito in delitos) {
    # delito=delitos[1]
    mlb <- readRDS(paste0("output/mlb_",city,"_",delito,".rds"))
    regresion <- tabula.coef(cityp = cityp, nb = nb, 
                             model = mlb, delito = delito)
    regressions <- merge(regressions, regresion, by = "Vars",
                         all.x = T, all.y = T, sort = F)    
  }
  # names(regressions) <- c("Variables", crimes)
  saveRDS(regressions, paste0("output_reg/mlb_",city,".rds"))
  
  # spatial autorregressive (slow prediction)
  regressions <- data.frame(Vars=NA)
  for (delito in delitos)  {
    esar <- readRDS(paste0("output/sar_",city,"_",delito,".rds"))
    regresion <- tabula.coef(cityp, nb, esar, delito)
    regressions <- merge(regressions, regresion, by = "Vars", 
                         all.x = T, all.y = T, sort = F) 
  }
  # names(regressions) <- c("Variables", crimes)
  saveRDS(regressions, paste0("output_reg/esar_",city,".rds"))
  
}

### Organize and export to Excel
statistics <- data.frame(Variables = c(crimes,"Facilities","Pop.density",
                                       "Household.size","Education.av",
                                       "Neigh.disadvantage","Neigh.mobility"))
# delitos_out <- c("viol","asalto")
delitos_out <- c("abus","crim", "viol")
indices  <- read.csv2(file = "docs/indice_var.csv", header = T, sep = ",")


for (city in ciudades) {
  
  # load city stats
  stats <- readRDS(paste0("data/",gsub("urb","stat",city),".rds"))
  stats <- stats %>% filter(!row.names(.) %in% delitos_out)
  stats[nrow(stats)+1,] <- rep(city, ncol(stats))
  
  indices <- indices %>% rbind(c("City", city))
  statistics <- left_join(tibble::rownames_to_column(stats),indices, 
            by=c("rowname" = "abrev")) %>% 
    select(Variables, everything(), -rowname) %>% 
    mutate(Variables = ifelse(is.na(Variables), "City", Variables))
  
  # stats$Variables <- c(crimes,"Pop.density","Household.size","Facilities",
  #                      "Education.av","Neigh.disadvantage",
  #                      "Neigh.mobility","City")
  # statistics <- merge(statistics, stats, by = "Variables", 
  #                     all.x=T, all.y=T, sort = F)
  
  
  # load model summaries
  sars <- readRDS(paste0("output_reg/esar_",city,".rds"))
  sars <- sars[-nrow(sars),]
  sars$Variables <- c("(Intercept)","Facilities","Pop.density",
                      "Household.size","Education.av.","Lag.Fac.","Lag.p.d.",
                      "Lag.hh.s.","Lag.e.a","AIC","Spatial.KS")
  mlbs <- readRDS(paste0("output_reg/mlb_",city,".rds"))
  mlbs <- mlbs[-nrow(mlbs),]
  mlbs$Variables <- c("(Intercept)","Lag.Y","Facilities","Pop.density",
                      "Household.size","Education.av.","Neigh.disadvantage",
                      "Neigh.mobility","AIC","Spatial.KS")
  summaries <- rbind(sars, mlbs)
  
  # save results
  write.xlsx(summaries, file="output_reg/results.xlsx", sheetName=city,
             append=T, rowNames=F, colNames=T)
}
write.xlsx(statistics, file="output_reg/results_stats.xlsx", 
           sheetName="Statistics", append=T, rowNames=F, colNames=T)

