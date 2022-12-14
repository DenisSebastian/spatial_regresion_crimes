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
cities <- c("Antofagasta", "Temuco")
ciudades <- c("urb_2_5", "urb_9_25")
cities <- c("Antofagasta", "Coquimbo-Serena", "Temuco")
ciudades <- c("urb_2_5", "urb_4_18", "urb_9_25")

# ciudades <- sub(".rds","",list.files("data",pattern="urb"))
# ciudades <- ciudades[!ciudades=="urb_13_3"]
# cities <- ciudades


####### 2. Crime Analysis and Models by city #######

###  Parallelization parameters
mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
ncores <- 4
cl <- makeCluster(ncores,type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
iter <- length(ciudades)

### Parallel analyses
models <- foreach(obs=1:iter,.options.multicore=mcoptions,.packages=c("spatialreg","MASS","lme4","spdep")) %dopar% {
  
city <- ciudades[obs]
cityp <- readRDS(paste0("../data/",city,".rds"))
cityp$rows <- rownames(cityp@data) <- rownames(cityp@coords)
saveRDS(cityp, paste0("../data/",city,".rds")) 

##### statistical checks #####

# # Definir variables
# colnames <- c("viol","crim","abus","amen","robo","hurto","asalto")
# 
# # Ciclo para graficos agrupados
# par(mfrow=c(3, 3))
# 
# for (v in colnames) {
#   hist(ciudad@data[,v], breaks=100, main=v, probability=TRUE, col="gray")
# }
# 
# # Definir variables
# colnames <- c("ofcom","educ","tamhog","denspob","vulzon","mobil")
# 
# # Ciclo para graficos agrupados
# par(mfrow=c(3, 3))
# for (v in colnames) {
#   hist(ciudad@data[,v], breaks=100, main=v, probability=TRUE, col="gray")
# }
# 
# # Definir variables
# colnames <- c("lpflot","lofcom","educ","tamhog","ldenspob","vulzon","mobil")
# 
# # Ciclo para graficos agrupados
# par(mfrow=c(3, 3))
# for (v in colnames) {
#   hist(ciudad@data[,v], breaks=100, main=v, probability=TRUE, col="gray")
# }

##### spatial regression #####

## spatial weights matrix
nvec <- 12
nb <- nb2listw(neighbours = knn2nb(
  knn = knearneigh( x = cityp, k = nvec, longlat = F)), 
  style = "W")

# cubic transformation
cityp$cubcrim   <- (cityp$crim)**(1/3)
cityp$cubabus   <- (cityp$abus)**(1/3)
cityp$cubamen   <- (cityp$amen)**(1/3)
cityp$cubrobo   <- (cityp$robo)**(1/3)
cityp$cubhurto  <- (cityp$hurto)**(1/3)
# logarithmic transformation
cityp$logcrim   <- log1p(cityp$crim)
cityp$logabus   <- log1p(cityp$abus)
cityp$logamen   <- log1p(cityp$amen)
cityp$logrobo   <- log1p(cityp$robo)
cityp$loghurto  <- log1p(cityp$hurto)
# lagged variables
cityp$lagcrim   <- as.numeric(lag.listw(x=nb, var=cityp$crim))
cityp$lagabus   <- as.numeric(lag.listw(x=nb, var=cityp$abus))
cityp$lagamen   <- as.numeric(lag.listw(x=nb, var=cityp$amen))
cityp$lagrobo   <- as.numeric(lag.listw(x=nb, var=cityp$robo))
cityp$laghurto  <- as.numeric(lag.listw(x=nb, var=cityp$hurto))

## Iteration by crime types
  
for (delito in delitos)
  {
    # spatial error and lag model
    sar <- spatialreg::errorsarlm(
        formula = formula(paste0("cub",delito,"~lofcom+ldenspob+tamhog+educ")),
        listw = nb,
        etype = "emixed",
        data = cityp@data)
  
    # multilevel binomial lagged model
    mlb <- glmer.nb(
      formula = formula(paste0(delito,"~lag",delito,"+lofcom+ldenspob+tamhog+educ+vulzon+mobil+(1|zona)")),
      data = cityp@data, 
      verbose=F)
    
    # save results
    saveRDS(sar, paste0("../output/sar_",city,"_",delito,".rds"))
    saveRDS(mlb, paste0("../output/mlb_",city,"_",delito,".rds"))
  }
}
## End parallel
stopCluster(cl)
