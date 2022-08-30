
# 1.  Setup ---------------------------------------------------------------

# Librerías y Recursos ----------------------------------------------------
source("scripts/recursos.R")
library(ecespa)

# Funciones ---------------------------------------------------------------
source("scripts/funciones.R")
cvm_eval <- function(model, coords, nperm = 999){
  predic_model <- predict(model, model@frame[, -1]) %>% abs()%>% as.vector()
  model_cvm <- syrjala0(coords = coords,R = F,
                       var1 = model@frame[, 1],
                       var2 = predic_model, nsim = nperm)
  return(model_cvm)
}


ciudades=gsub(".rds","",list.files("data",pattern="urb"))
# ciudades=ciudades[!ciudades=="urb_13_3"]
cities=ciudades



# cvm_test <- goftest::cvm.test(x = pred, estimated = T)$statistic
# str(cvm_test)


#### List of crimes (delitos) and cities (ciudades)
crimes=c("Minor property offenses","Robberies","Burglaries","Drunkennes, damages and disorders",
         "Family violence and aggression","Injuries, drugs and weapons","High violence and murder")
delitos=c("hurto","robo","asalto","amen","abus","crim","viol")

### Creation of Data structures for results gathering
stats=data.frame(matrix(c("1","2","viol","crim","abus","amen","asalto","robo","hurto",
                          "pobdens","tamhog","flux","ofcom","educ","vuln","mobil"),16,1))
names(stats)="Variables"

syntstat=data.frame()

vars=data.frame(matrix(c("City","Significant Coef","Lagged Crime","Log(Flux Pop)","Log(Dens Pop)","Log(Dens Facilities)",
                         "Household Size","Education","Vulnerability (Zone)","New residents (Zone)"),20,1))
names(vars)="Variables"

regs=data.frame()

bests=data.frame(matrix(c("City","Crime Category","Model","Significant variables","BIC","logLik","KS statistic"),7,1))
names(bests)="Desc"

mls=data.frame(matrix(c("City","Crime Category","(Intercept)","Lagged Crime",
                        "Log(Flux Pop)","Log(Dens Pop)","Log(Dens Facilities)",
                        "Household Size","Education","Vulnerability (Zone)","New residents (Zone)",
                        "BIC","logLik","KS statistic") ,14,1))
names(mls)="Desc"


names_type <- c("Center", "Flux", "Center.Controls",
  "Flux.Controls", "Center Flux" )

ciudades_list <- list()

## for every city
for (ciud in 14) {#:length(ciudades
  # for (ciud in  c(1,2,12,17,24)) {
  ciudad=ciudades[ciud]
  city=cities[ciud]
  cityp=readRDS(paste0("data/",ciudad,".rds"))
  citystat=readRDS(paste0("data/",gsub("urb","stat",ciudad),".rds"))
  cityst=c(ciudad,citystat$n[8],as.integer(citystat$n[11]),as.integer(citystat$n[10]))
  citystat=rbind(rep(city,3),c("N","Mean","StD"),citystat)
  stats=cbind(stats,citystat)
  varpos=0
  varneg=0
  print(ciudad)
  coords <- spTransform(cityp, CRS("+init=epsg:32719")) 
  coords <- coords@coords %>% as.data.frame()
  ## for every crime
  delitos_list <- list()
  for (del in 1:7){
    delito=delitos[del]
    crime=crimes[del]
    print(delito)
    
    ### read all models for this crime type and city
    lagcent=readRDS(paste0("output/",ciudad,"_lgct_",delito,".rds"))
    lagflux=readRDS(paste0("output/",ciudad,"_lgfl_",delito,".rds"))
    cpcent=readRDS(paste0("output/",ciudad,"_ctcn_",delito,".rds"))
    cpflux=readRDS(paste0("output/",ciudad,"_flcn_",delito,".rds"))
    centflux=readRDS(paste0("output/",ciudad,"_ctfl_",delito,".rds"))
    
    # cvm=data.frame(
    

    lagcent_cvm <-  cvm_eval(model = lagcent, coords = coords)
    lagflux_cvm <-  cvm_eval(model = lagflux, coords = coords)
    cpcent_cvm <-  cvm_eval(model = cpcent, coords = coords)
    cpflux_cvm <-  cvm_eval(model = cpflux, coords = coords)
    centflux_cvm <-  cvm_eval(model = centflux, coords = coords)
    
    resultados <- list(lagcent_cvm, lagflux_cvm, cpcent_cvm, cpflux_cvm,
                        centflux_cvm)
    nombres_res <-  c("lagcent_cvm", "lagflux_cvm", "cpcent_cvm", "cpflux_cvm",
                      "centflux_cvm")
    
    names(resultados) <- nombres_res
    
    delitos_list[delito] <- resultados
    
    # 
    # 
    # names(cvm) <- names_type
    # registro <- data.frame(Ciudad=ciudad, Delito =delito, Test = "CVMtest", cvm)
    # 
    # cvm_registros <-  cvm_registros %>%
    #   rbind(registro)
    
  }
  ciudades_list[ciudad] <- delitos_list
}



rownames(cvm_registros) <- NULL
promedios <- cvm_registros %>% 
  group_by(Delito) %>% 
  summarise(Center_m = mean(Center),
            Flux_m = mean(Flux),
            Center.Controls_m = mean(Center.Controls),
            Flux.Controls_m = mean(Flux.Controls),
            Center.Flux_m = mean(Center.Flux))

saveRDS(promedios, "output/mean_city_cvm.rds")
saveRDS(cvm_registros, "output/registros_cvm.rds")
saveRDS(ciudades_list, "output/ciudades_cvm.rds")

# 
# data <- data.frame(x1 = c(0, 1, 0, 7, 1),             # Creating example data
#                    x2 = c(5, 0, 0, 0, 0),
#                    x3 = 3)
# data    
# colnames(data)[max.col(data, ties.method = "first")] 
# 
# print(colnames(reg)[max.col(data)])
# reg
# 
# res <- reg[(nrow(reg)-3):nrow(reg), ] %>% 
#   t() %>% 
#   `colnames<-`(.[1, ]) %>%
#   .[-1, ] %>% 
#   as.data.frame() %>% 
#   mutate_if(is.character, as.numeric)
# 
# 
# colnames(res)[max.col(res, ties.method = "first")] 
# 
# library(dplyr)
# 
# res[which.max(res),]
# 
# cvm.test_2 <- function(x,y){
#   r <- 1000 #permutation samples
#   reps <- numeric(r) #score of the two-sample Cramer-von Mises test
#   n <- length(x)
#   m <- length(y)
#   v.n <- numeric(n)
#   v1.n <- numeric(n)
#   v.m <- numeric(m)
#   v1.m <- numeric(m)
#   z <- c(x,y) #the combined sample
#   N <- length(z)
#   Ix <- seq(1:n)
#   Iy <- seq(1:m)
#   v.n <- (x-Ix)**2
#   v.m <- (y-Iy)**2
#   #test statistic
#   reps_0 <- ((n * sum(v.n)+m * sum(v.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
#   for (k in 1:r){#permutation samples
#     w <- sample(N,size=n,replace=FALSE)
#     x1 <- sort(z[w])
#     y1 <- sort(z[-w])
#     v1.n <- (x1-Ix)**2
#     v1.m <- (y1-Iy)**2
#     reps[k] <- ((n * sum(v1.n)+m * sum(v1.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
#   }
#   p <- mean(c(reps_0,reps) >= reps_0)
#   res <- list(reps=reps,reps0=reps_0,pvalue=p)
#   return(res)
# }
