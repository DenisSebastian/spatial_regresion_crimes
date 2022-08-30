


# Directorio de Trabajo ---------------------------------------------------
  setwd("~/OneDrive - Universidad Adolfo Ibanez/FONDECYT/Crime_Cellphones")
  
# Librerías y Recursos ----------------------------------------------------
  source("scripts/recursos.R")
  
  
  # Funciones ---------------------------------------------------------------
  source("scripts/funciones.R")
  
  # Definición de Delitos ---------------------------------------------------
  delitos = c("hurto", "robo", "asalto", "amen", "abus", "crim", "viol")
  
 # Lectura de Insumos ------------------------------------------------------
  
  ciudades = gsub(".rds", "", list.files("data", pattern = "urb"))
  ciudades=ciudades[ciudades=="urb_4_18"]
  # cities = ciudades
  city <- ciudades
  cityp <- readRDS(paste0("data/",city,".rds"))
  cityp$rows <- rownames(cityp@data) <- rownames(cityp@coords)
  

# Matriz de Pesos Espaciales ----------------------------------------------

  nb <- spatial_weights(city_df = cityp, nvec= 12)
  

# Transformaciones  -------------------------------------------------------

  cityp <- cityp %>%
    st_as_sf(crs = 4326) %>% 
    mutate(across(.cols = delitos, .fns = t_cub, .names = "cub{.col}")) %>% 
    mutate(across(.cols = delitos, .fns = t_log1p, .names = "log{.col}")) %>% 
    mutate(across(.cols = delitos, .fns = lag_spatial, nb, .names = "lag{.col}")) %>%
    as("Spatial")
    
  

# Modelos -----------------------------------------------------------------

  ## Iteration by crime types
  # numCores <- detectCores()
  numCores <- length(delitos)
  numCores <- 3
  registerDoParallel(numCores) 
  
  foreach(delito=all_of(delitos), .export = c("cityp", "nb"),
           packages=c("spatialreg","MASS","lme4","spdep")) %dopar% {
             print(paste0("Aplicando modelos a delito: ", delito))
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
             saveRDS(sar, paste0("output/sar_",city,"_",delito,".rds"))
             saveRDS(mlb, paste0("output/mlb_",city,"_",delito,".rds"))
           }
  
  # When you're done, clean up the cluster
  stopImplicitCluster()
