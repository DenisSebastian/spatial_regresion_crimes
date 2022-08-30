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

## geographic projections
crs_latlon <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"



# Parameters
iter <- length(ciudades)


foreach(obs=1:iter, .packages=c("spatialreg","MASS","BAMMtools","lme4","spdep")) %do% {
  # obs <- 1
  
  
  city <- ciudades[obs] # city=ciudades[1]
  cityp <- readRDS(paste0("data/",city,".rds"))
  # cityp <- spTransform(cityp, CRS(crs_utm))
  # cityp$rows <- rownames(cityp@data) <- rownames(cityp@coords)
  # nb <- spatial_weights(city_df = cityp, nvec= 12)
  # coords <- cityp@coords %>% as.data.frame()
  
  #make windows
  ext <- extent(cityp)
  x_min <- ext[1] + 100
  x_max <- ext[2] - 100
  y_min <- ext[3] + 100
  y_max <- ext[4] - 100
  w <- as.owin(c(x_min,x_max, y_min, y_max))
  
  
  cityp <- cityp %>%
    st_as_sf(crs = 4326) %>% 
    mutate(across(.cols = all_of(delitos), .fns = t_cub, .names = "cub{.col}")) %>% 
    mutate(across(.cols = all_of(delitos), .fns = t_log1p, .names = "log{.col}")) %>% 
    mutate(across(.cols = all_of(delitos), .fns = lag_spatial, nb, .names = "lag{.col}")) %>%
    as("Spatial")
  
  foreach(delito=all_of(delitos), .export = c("cityp", "nb"), 
          packages=c("spatialreg","MASS","lme4","spdep")) %do% {
            
            # delito <-  delitos[1]
    print(paste0("Aplicando modelos a delito: ", delito)) 
  
    model_mlb <- readRDS(paste0("output/mlb_",city,"_",delito,".rds"))
    
    model <- model_mlb

    
   


  ppp_mbl_pred <- ppp_maker(model = model_mlb, city = cityp, 
                       delito = delito, ppp_predict = T, w = w)
  
  plot(density(ppp_mbl_pred, adjust=.25))
  
  ppp_mbl <- ppp_maker(model = model_mlb, city = cityp, 
                            delito = delito, ppp_predict = F, w= w)
  
  plot(density(ppp_mbl, adjust=.25))
  
  # st <- syrjala.test(ppp_mbl, ppp_mbl_pred, nsim = 999)
  cvm_mlb <- cdf.test(ppp_mbl, as.im(ppp_mbl), test="cvm")
  plot(cvm_mlb)


  model_sar <- readRDS(paste0("output/sar_",city,"_",delito,".rds"))
  ppp_sar_pred <- ppp_maker(model = model_sar, city = cityp, 
                            delito = delito, ppp_predict = T,
                            w = w, nb = nb)
  
  # plot(density(ppp_sar_pred, adjust=.1))
  
  ppp_sar <- ppp_maker(model = model_sar, city = cityp, 
                       delito = delito, ppp_predict = F, 
                       w = w, nb = nb)
  
  # plot(density(ppp_sar, adjust=.25))
  
  # st <- syrjala.test(ppp_mbl, ppp_mbl_pred, nsim = 999)
  cvm_sar <- cdf.test(ppp_sar, as.im(ppp_sar_pred), test="cvm")
  plot(cvm_sar)
  
  
  #guardar resultados
  saveRDS(cvm_mlb, "cvm_results/cdf.test_mlb.rds")
  saveRDS(cvm_sar, "cvm_results/cdf.test_sar.rds")
  
          }

}



