####################################
###### Crime models comparison #####
####################################

####### 1. Setup #######

options(scipen = 999, stringsAsFactors=F, encoding = "UTF-8", width = 1000)

## Check, download and load packages
list.packages <- c("sp","leaflet","DescTools","spatialreg","MASS","BAMMtools","lme4","xlsx","spdep","rstudioapi","doSNOW","doParallel")
new.packages <- list.packages[!(list.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(sp))
suppressMessages(library(leaflet))
suppressMessages(library(DescTools))
suppressMessages(library(spatialreg))
suppressMessages(library(MASS))
suppressMessages(library(BAMMtools))
suppressMessages(library(lme4))
suppressMessages(library(xlsx))
suppressMessages(library(spdep))
suppressMessages(library(rstudioapi))
suppressMessages(library(doSNOW))
suppressMessages(library(doParallel))

## Set path to script's location
setwd(dirname(getActiveDocumentContext()$path)) 

## List of crimes/delitos and cities/ciudades
crimes <- c("Minor property offenses", "Robberies", "Drunkennes, damages and disorders",
         "Family violence and aggression","Injuries, drugs and weapons")
delitos <- c("hurto", "robo", "amen", "abus", "crim")

cities <- c("Antofagasta", "Coquimbo-Serena", "Temuco")
ciudades <- c("urb_2_5", "urb_4_18", "urb_9_25")

# ciudades <- sub(".rds","",list.files("data",pattern="urb"))
# ciudades <- ciudades[!ciudades=="urb_13_3"]
# cities <- ciudades

### Test and compilation functions

# Coefficients standardization function
betam <-  function (model)
{
   if (class(model)[1] == "sarlm")
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
ks.spat <-  function (cityp, nb, model, delito)
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
  names(fits) <- gsub("y",paste0(mod,delito),names(fits),fixed=T)
  # nrow(fits[fits$ycfit==0,])/nrow(fits)
  return(fits)
}

# Results compilation function
tabula.coef <- function(cityp, nb, model, delito)
{
  # data extraction
  if (class(model)[1] == "sarlm")
  {
    tab <- data.frame(summary(model)$Coef)
  }
  if (class(model)[1] == "glmerMod")
  { 
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
  if (class(model)[1] == "sarlm")
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
  SKS <- ks.spat(cityp, nb, model, delito)

  # compile
  tab <- rbind(tab, errores, SKS)
  names(tab)[2] <- delito
  return(tab)
}

### Parallell compilation of statistics and regression models

# Parameters
nvec <- 12
mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
ncores <- 4
cl <- makeCluster(ncores,type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
iter <- length(ciudades)

# Models compilation
models <- foreach(obs=1:iter,.options.multicore=mcoptions,.packages=c("spatialreg","MASS","BAMMtools","lme4","spdep")) %dopar% {

  city <- ciudades[obs] # city=ciudades[1]
  cityp <- readRDS(paste0("../data/",city,".rds"))
  nb <- nb2listw(neighbours = knn2nb(
    knn = knearneigh( x = cityp, k = nvec, longlat = F)), 
    style = "W")

  # multilevel lagged binomial
  regressions <- data.frame(Vars=NA)
  for (delito in delitos) {
    # delito=delitos[1]
    mlb <- readRDS(paste0("../output/mlb_",city,"_",delito,".rds"))
    regresion <- tabula.coef(cityp, nb, mlb, delito)
    regressions <- merge(regressions, regresion, by = "Vars", all.x = T, all.y = T, sort = F)    
  }
  names(regressions) <- c("Variables", crimes)
  saveRDS(regressions, paste0("../output/mlb_",city,".rds"))
  
  # spatial autorregressive (slow prediction)
  regressions <- data.frame(Vars=NA)
  for (delito in delitos)  {
    sar <- readRDS(paste0("../output/sar_",city,"_",delito,".rds"))
    regresion <- tabula.coef(cityp, nb, sar, delito)
    regressions <- merge(regressions, regresion, by = "Vars", all.x = T, all.y = T, sort = F) 
  }
  names(regressions) <- c("Variables", crimes)
  saveRDS(regressions, paste0("../output/sar_",city,".rds"))

}
# End parallel
stopCluster(cl)

### Organize and export to Excel
statistics <- data.frame(matrix(c(crimes,"Facilities","Pop.density","Household.size","Education.av","Neigh.disadvantage","Neigh.mobility"),11,1))
names(statistics) <- "Variables"
for (city in ciudades) {
  
  # load city stats
  stats <- readRDS(paste0("../data/",gsub("urb","stat",city),".rds"))
  stats <- stats[-c(1,5),]
  stats[nrow(stats)+1,] <- rep(city, ncol(stats))
  stats$Variables <- c(crimes,"Pop.density","Household.size","Facilities","Education.av","Neigh.disadvantage","Neigh.mobility","City")
  statistics <- merge(statistics, stats, by = "Variables", all.x=T, all.y=T, sort = F)
  
  # load model summaries
  sars <- readRDS(paste0("../output/sar_",city,".rds"))
  sars <- sars[-nrow(sars),]
  sars$Variables <- c("(Intercept)","Facilities","Pop.density","Household.size","Education.av.","Lag.Fac.","Lag.p.d.","Lag.hh.s.","Lag.e.a","AIC","Spatial.KS")
  mlbs <- readRDS(paste0("../output/mlb_",city,".rds"))
  mlbs <- mlbs[-nrow(mlbs),]
  mlbs$Variables <- c("(Intercept)","Lag.Y","Facilities","Pop.density","Household.size","Education.av.","Neigh.disadvantage","Neigh.mobility","AIC","Spatial.KS")
  summaries <- rbind(sars, mlbs)
  
  # save results
  write.xlsx(summaries, file="../output/results.xlsx", sheetName=city, append=T, row.names=F, col.names=T)
}
write.xlsx(statistics, file="../output/results.xlsx", sheetName="Statistics", append=T, row.names=F, col.names=T)

