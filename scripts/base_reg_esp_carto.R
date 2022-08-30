####################################
###### Crime models comparison #####
####################################

####### 1. Setup #######

options(scipen = 999, stringsAsFactors=F, encoding = "UTF-8", width = 1000)

## Check, download and load packages
list.packages <- c("data.table","sp","raster","rgdal","ggplot2","ggmap","RColorBrewer","scales","gridExtra","leaflet","DescTools","spatialreg","MASS","BAMMtools","lme4","xlsx","spdep","rstudioapi","doSNOW","doParallel")
new.packages <- list.packages[!(list.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(data.table))
suppressMessages(library(sp))
suppressMessages(library(raster))
suppressMessages(library(rgdal))
suppressMessages(library(ggplot2))
suppressMessages(library(ggmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(gridExtra))
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

## geographic projections
crs_latlon <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## List of crimes/delitos and cities/ciudades
crimes <- c("Minor property offenses", "Robberies", "Drunkennes, damages and disorders",
         "Family violence and aggression","Injuries, drugs and weapons")
delitos <- c("hurto", "robo", "amen", "abus", "crim")

cities <- c("Antofagasta", "Coquimbo-Serena", "Temuco")
ciudades <- c("urb_2_5", "urb_4_18", "urb_9_25")

# ciudades <- sub(".rds","",list.files("data",pattern="urb"))
# ciudades <- ciudades[!ciudades=="urb_13_3"]
# cities <- ciudades

# Spatial prediction error (SPE) function
ks.spat.preds <-  function (cityp, nb, model, delito)
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


####### 2. Compile spatial predictions #######

# Parameters
nvec <- 12
# mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
# ncores <- 4
# cl <- makeCluster(ncores,type="SOCK")
# registerDoParallel(cl)
# getDoParWorkers()
# iter <- length(ciudades)

# Models compilation
# models <- foreach(obs=1:iter,.options.multicore=mcoptions,.packages=c("spatialreg","MASS","BAMMtools","lme4","spdep")) %dopar% {

  city <- ciudades[obs] # city=ciudades[3]

for (city in ciudades){
  cityp <- readRDS(paste0("../data/",city,".rds"))
  cityp$rows <- rownames(cityp@data) <- rownames(cityp@coords)
  nb <- nb2listw(neighbours = knn2nb(
    knn = knearneigh( x = cityp, k = nvec, longlat = F)), 
    style = "W")

  # multilevel lagged binomial
  for (delito in delitos) {
    # delito=delitos[1]
    mlb <- readRDS(paste0("../output/mlb_",city,"_",delito,".rds"))
    cityp <- sp::merge(cityp, ks.spat.preds(cityp, nb, mlb, delito), by="rows")
  }
  
  # spatial autorregressive 
  for (delito in delitos)  {
    sar <- readRDS(paste0("../output/sar_",city,"_",delito,".rds"))
    cityp <- sp::merge(cityp, ks.spat.preds(cityp, nb, sar, delito), by="rows")
  }
  saveRDS(cityp, paste0("../carto/",city,"_preds.rds"))
}
# End parallel
# stopCluster(cl)


####### 3. Cartographic analysis #####

for (city in ciudades) { 
  # city <- ciudades[1]
  
  cityp <- readRDS(paste0("../carto/",city,"_preds.rds"))
  cityp <- spTransform(cityp, crs_latlon)
    
  # slight zoom out
  box <- cityp@bbox 
  box[1,1] <- cityp@bbox[1,1] - abs(cityp@bbox[1,1]-cityp@bbox[1,2])/20
  box[1,2] <- cityp@bbox[1,2] + abs(cityp@bbox[1,1]-cityp@bbox[1,2])/20
  box[2,1] <- cityp@bbox[2,1] - abs(cityp@bbox[2,1]-cityp@bbox[2,2])/20
  box[2,2] <- cityp@bbox[2,2] + abs(cityp@bbox[2,1]-cityp@bbox[2,2])/20
  
  fondo <- get_stamenmap( bbox = box, zoom = 13, maptype = "terrain", color="bw", force=T)
  mapa <- ggmap(fondo,  extent = "device", legend = "topleft")

  num <- 0
  
  if (city == ciudades[1]){ 
    ancho = 36
    alto = 44
    }
  if (city == ciudades[2]){ 
    ancho = 44
    alto = 34
    }
  if (city == ciudades[3]){ 
    ancho = 48
    alto = 32
    }

for (delito in delitos) {
  # delito <- delitos[1]
  num <- num+1
  
  cityp$delito <- cityp@data[,delito]
  cityp$errsar <- cityp@data[,paste0("sar",delito,"cfit")]
  cityp$errmlb <- cityp@data[,paste0("mlb",delito,"cfit")]
  
  errorsar <- round(nrow(cityp@data[cityp$errsar==0,])/nrow(cityp@data),2)
  errormlb <- round(nrow(cityp@data[cityp$errmlb==0,])/nrow(cityp@data),2)
    
  crimerep <- data.table(cityp@data)[,list(x=rep(x,delito),y=rep(y,delito))]
  sarrep <- cityp@data[cityp$errsar!=0,c("x","y","errsar")]
  mlbrep <- cityp@data[cityp$errmlb!=0,c("x","y","errmlb")]

  mapcrim <- mapa +
    stat_density2d(aes(x = x, y = y, fill = ..level../100, alpha = ..level..),
                 bins = 9, geom = "polygon", adjust = 1,
                 data = crimerep) +
    guides(fill = guide_legend("Crime density"), alpha = F) +
    scale_alpha_continuous(range=c(0.5,0.8)) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    theme_void() + 
    ggtitle(paste0(crimes[num])) + 
    theme(
      plot.title = element_text(size = 18), 
      legend.text=element_text(size = 11),
      panel.border = element_rect(colour = "grey", fill=NA, size=1)
    )

  mapsar <- mapa +
    geom_point(aes(x = x, y = y, colour = factor(round(errsar/3,0)), size = round(abs(errsar)/9,0)), data = sarrep, alpha = 0.7) +
    guides(colour=guide_legend("Sp Err SAR"), size = F) +
    scale_colour_brewer(palette = "Spectral") +
    theme_void() + 
    ggtitle(paste0("Hits SAR: ", errorsar)) + 
    theme(
      plot.title = element_text(size = 14),
      legend.text=element_text(size = 11),
      panel.border = element_rect(colour = "grey", fill=NA, size=1)
    )

  mapmlb <- mapa +
    geom_point(aes(x = x, y = y, colour = factor(round(errmlb/3,0)), size = round(abs(errmlb)/9,0)), data = mlbrep, alpha = 0.7) +
    guides(colour=guide_legend("Sp Err MLB"), size = F) +
    scale_colour_brewer(palette = "Spectral") +
    theme_void() + 
    ggtitle(paste0("Hits MLB: ", errormlb)) + 
    theme(
      plot.title = element_text(size = 14),
      legend.text=element_text(size = 11),
      panel.border = element_rect(colour = "grey", fill=NA, size=1)
    ) 
  
  p <- grid.arrange(mapcrim, mapsar, mapmlb, nrow = 1, heights = c(1), widths = c(1,1,1))
  ggsave( paste0("../carto/",city,delito,".png"),
          plot = p,
          width = ancho, height = alto, units = "cm", dpi = 150,)
  
  }
}

