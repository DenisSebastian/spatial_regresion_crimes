
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

cities <- c("Coquimbo-Serena")
ciudades <- c("urb_4_18")

## geographic projections
crs_latlon <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"



# Parameters
iter <- length(ciudades)


foreach(obs=1:iter, .packages=c("spatialreg","MASS","BAMMtools","lme4","spdep")) %do% {
  obs <- 1
  
  
  city <- ciudades[obs] # city=ciudades[1]
  cityp <- readRDS(paste0("data/",city,".rds"))
  cityp <- spTransform(cityp, CRS(crs_utm))
  nb <- spatial_weights(city_df = cityp, nvec= 12)

  
  
  cityp <- cityp %>%
    st_as_sf(crs = 4326) %>% 
    mutate(across(.cols = all_of(delitos), .fns = t_cub, .names = "cub{.col}")) %>% 
    mutate(across(.cols = all_of(delitos), .fns = t_log1p, .names = "log{.col}")) %>% 
    mutate(across(.cols = all_of(delitos), .fns = lag_spatial, nb, .names = "lag{.col}")) %>%
    as("Spatial")
  
  foreach(delito=all_of(delitos), .export = c("cityp", "nb"), 
          packages=c("spatialreg","MASS","lme4","spdep")) %do% {
            
      delito <-  delitos[1]
      print(paste0("Aplicando modelos a delito: ", delito)) 
      
      
      # mlb
      model_mlb <- readRDS(paste0("output/mlb_",city,"_",delito,".rds"))
      mlb_df <- predict_df(model = model_mlb, city = cityp, longlat = F)
      
     
      res_mlb <- mlb_df %>% 
        slice(sample(1:nrow(.), 100)) %>% 
        mutate(
          cvm.sim = cvm_eval_spatial(coords_x = x, coords_y = y,
                                     var1 = delitos, var2 = ypred, 
                                     nperm = nrow(.)-1))
      
      head(res_mlb)


      
      
      #sar
      model_sar <- readRDS(paste0("output/sar_",city,"_",delito,".rds"))
      df_sar <- predict_df(model = model_mlb, city = cityp, longlat = F)
      
      res_sar <- df_sar %>% 
        slice(sample(1:nrow(.), 100)) %>% 
        mutate(
          cvm.sim = cvm_eval_spatial(coords_x = x, coords_y = y,
                                     var1 = delitos, var2 = ypred, 
                                     nperm = nrow(.)-1))

      
      #guardar resultados
      saveRDS(cvm_mlb, "cvm_results/cdf.test_mlb.rds")
      saveRDS(cvm_sar, "cvm_results/cdf.test_sar.rds")
            
          }
  
}



# res_mlb_sf <- st_as_sf(x = res_mlb,
#                        coords = c("x", "y"),
#                        crs = 4326)
#
# mapview::mapview(res_mlb_sf, zcol ="cvm.sim") +
# mapview::mapview(res_mlb_sf, zcol ="ypred_adj")+
#   mapview::mapview(res_mlb_sf, zcol ="y_adj")

# saveRDS(res_sar, "cvm_results/res_cvm_sar.rds")

