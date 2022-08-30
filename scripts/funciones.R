

# Funciones ---------------------------------------------------------------

## spatial weights matrix
spatial_weights <-  function(city_df, nvec = 12,
                             style = "W", longlat = F){
  nb <- spdep::nb2listw(neighbours = spdep::knn2nb(
    knn = spdep::knearneigh( x = city_df, k = nvec, longlat = F)), 
    style = style)
  return(nb)
}

# transformaciones
t_cub <-  function(x) x**(1/3)
t_log1p <-  function(x) log1p(x)
lag_spatial <-  function(x, nb){
  as.numeric(spdep::lag.listw(x=nb, var=x))
}



lm.betam <-  function (MOD) # extract standarized coefficients for multilevel models
{
  b <- summary(MOD)$coefficients[-1, 1]
  sx <- apply(MOD@frame[2:(ncol(MOD@frame)-1)],2,FUN=sd)
  sy <- sd(unlist(MOD@frame[1]))
  beta <- b * sx/sy
  return(beta)
}

tabula_coef = function(cityp,lagcent,lagflux,cpcent,cpflux,centflux)   ## receives data and 5 models
{
  ### tabulate the results
  tab1 = data.frame(summary(lagcent)$coefficients)
  tab2 = data.frame(summary(lagflux)$coefficients)
  tab3 = data.frame(summary(cpcent)$coefficients)
  tab4 = data.frame(summary(cpflux)$coefficients)
  tab5 = data.frame(summary(centflux)$coefficients)
  
  ## standarize coefficients
  tab1[-1,1] = lm.betam(lagcent)
  tab2[-1,1] = lm.betam(lagflux)
  tab3[-1,1] = lm.betam(cpcent)
  tab4[-1,1] = lm.betam(cpflux)
  tab5[-1,1] = lm.betam(centflux)
  
  ## record statistical significance  Note:  *p<0.1; **p<0.05; ***p<0.01
  tab1$significancia = ifelse(tab1[,4]<0.01,"***",ifelse(tab1[,4]<0.05,"**",ifelse(tab1[,4]<0.1,"*","")))
  tab2$significancia = ifelse(tab2[,4]<0.01,"***",ifelse(tab2[,4]<0.05,"**",ifelse(tab2[,4]<0.1,"*","")))
  tab3$significancia = ifelse(tab3[,4]<0.01,"***",ifelse(tab3[,4]<0.05,"**",ifelse(tab3[,4]<0.1,"*","")))
  tab4$significancia = ifelse(tab4[,4]<0.01,"***",ifelse(tab4[,4]<0.05,"**",ifelse(tab4[,4]<0.1,"*","")))
  tab5$significancia = ifelse(tab5[,4]<0.01,"***",ifelse(tab5[,4]<0.05,"**",ifelse(tab5[,4]<0.1,"*","")))
  
  ## rename
  tab1$var = row.names(tab1)
  tab2$var = row.names(tab2)
  tab3$var = row.names(tab3)
  tab4$var = row.names(tab4)
  tab5$var = row.names(tab5)
  
  ## synthesize info for summary regression estimate ** (std dev)
  tab1$resumen = paste0(round(tab1$Estimate,2),tab1$significancia,"_(",round(tab1$Std..Error,2),")")
  tab2$resumen = paste0(round(tab2$Estimate,2),tab2$significancia,"_(",round(tab2$Std..Error,2),")")
  tab3$resumen = paste0(round(tab3$Estimate,2),tab3$significancia,"_(",round(tab3$Std..Error,2),")")
  tab4$resumen = paste0(round(tab4$Estimate,2),tab4$significancia,"_(",round(tab4$Std..Error,2),")")
  tab5$resumen = paste0(round(tab5$Estimate,2),tab5$significancia,"_(",round(tab5$Std..Error,2),")")
  
  ## merge all tables
  coeff = merge(merge(tab1[,c("var","resumen")],tab2[,c("var","resumen")],by="var",all=TRUE),merge(tab3[,c("var","resumen")],merge(tab4[,c("var","resumen")],tab5[,c("var","resumen")],by="var",all=TRUE),by="var",all=TRUE),by="var",all=TRUE)
  colnames(coeff) = c("Vars","Center","Flux","Center Controls","Flux Controls","Center Flux")
  
  ## order rows
  orden=names(summary(centflux)$coefficients[, 1])
  coeff=coeff[match(orden,coeff$Vars),]
  
  ## remove NAs
  coeff[is.na(coeff)] = "."
  
  ## compile error stats for all models
  errores = data.frame(t(rbind(
    summary(lagcent)$AICtab,
    summary(lagflux)$AICtab,
    summary(cpcent)$AICtab,
    summary(cpflux)$AICtab,
    summary(centflux)$AICtab
  )))
  
  ## order names
  errores$indicador = row.names(errores)
  colnames(errores) = c("Center","Flux","Center Controls","Flux Controls","Center Flux","Vars")
  errores = errores[,c("Vars","Center","Flux","Center Controls","Flux Controls","Center Flux")]
  
  ##  Kolmogoronov - Smirnoff statistic
  ks=c(
    "KStest",
    ks.test(cityp@data[,delito],round(exp(predict(lagcent)),0))$statistic,
    ks.test(cityp@data[,delito],round(exp(predict(lagflux)),0))$statistic,
    ks.test(cityp@data[,delito],round(exp(predict(cpcent)),0))$statistic,
    ks.test(cityp@data[,delito],round(exp(predict(cpflux)),0))$statistic,
    ks.test(cityp@data[,delito],round(exp(predict(centflux)),0))$statistic
  )
  
  ##  Cramer-von Mises statistic
  cvm=c(
    "CVMtest",
    goftest::cvm.test(x = round(exp(predict(lagcent)),0), estimated = T)$statistic,
    goftest::cvm.test(x = round(exp(predict(lagflux)),0), estimated = T)$statistic,
    goftest::cvm.test(x = round(exp(predict(cpcent)),0), estimated = T)$statistic,
    goftest::cvm.test(x = round(exp(predict(cpflux)),0), estimated = T)$statistic,
    goftest::cvm.test(x = round(exp(predict(centflux)),0), estimated = T)$statistic
  )
  
  
  ##  Spatial pseudo-KS
  dato=cityp@data[,delito]
  # 1
  pred=round(exp(predict(lagcent)))
  fits=data.frame(dato,pred)
  fits$catdato[fits$dato==0]=0
  fits$catdato[fits$dato>0]=as.numeric(cut(fits$dato[fits$dato>0],breaks=4))
  fits$catpred[fits$pred==0]=0
  if (max(fits$pred)>0) {fits$catpred[fits$pred>0]=as.numeric(cut(fits$dato[fits$pred>0],breaks=4))}
  fits$catfit=abs(fits$catdato-fits$catpred)
  fits$catfit[fits$catfit>1]=1
  lagcentfit0=1-sum(fits$catfit)/length(fits$catfit)
  lagcentfit1=1-sum(fits$catfit[fits$dato>0])/length(fits$catfit[fits$dato>0])
  # 2
  pred=round(exp(predict(lagflux)))
  fits=data.frame(dato,pred)
  fits$catdato[fits$dato==0]=0
  fits$catdato[fits$dato>0]=as.numeric(cut(fits$dato[fits$dato>0],breaks=4))
  fits$catpred[fits$pred==0]=0
  if (max(fits$pred)>0) {fits$catpred[fits$pred>0]=as.numeric(cut(fits$dato[fits$pred>0],breaks=4))}
  fits$catfit=abs(fits$catdato-fits$catpred)
  fits$catfit[fits$catfit>1]=1
  lagfluxfit0=1-sum(fits$catfit)/length(fits$catfit)
  lagfluxfit1=1-sum(fits$catfit[fits$dato>0])/length(fits$catfit[fits$dato>0])
  # 3
  pred=round(exp(predict(cpcent)))
  fits=data.frame(dato,pred)
  fits$catdato[fits$dato==0]=0
  fits$catdato[fits$dato>0]=as.numeric(cut(fits$dato[fits$dato>0],breaks=4))
  fits$catpred[fits$pred==0]=0
  if (max(fits$pred)>0) {fits$catpred[fits$pred>0]=as.numeric(cut(fits$dato[fits$pred>0],breaks=4))}
  fits$catfit=abs(fits$catdato-fits$catpred)
  fits$catfit[fits$catfit>1]=1
  cpcentfit0=1-sum(fits$catfit)/length(fits$catfit)
  cpcentfit1=1-sum(fits$catfit[fits$dato>0])/length(fits$catfit[fits$dato>0])
  # 4
  pred=round(exp(predict(cpflux)))
  fits=data.frame(dato,pred)
  fits$catdato[fits$dato==0]=0
  fits$catdato[fits$dato>0]=as.numeric(cut(fits$dato[fits$dato>0],breaks=4))
  fits$catpred[fits$pred==0]=0
  if (max(fits$pred)>0) {fits$catpred[fits$pred>0]=as.numeric(cut(fits$dato[fits$pred>0],breaks=4))}
  fits$catfit=abs(fits$catdato-fits$catpred)
  fits$catfit[fits$catfit>1]=1
  cpfluxfit0=1-sum(fits$catfit)/length(fits$catfit)
  cpfluxfit1=1-sum(fits$catfit[fits$dato>0])/length(fits$catfit[fits$dato>0])
  # 5
  pred=round(exp(predict(centflux)))
  fits=data.frame(dato,pred)
  fits$catdato[fits$dato==0]=0
  fits$catdato[fits$dato>0]=as.numeric(cut(fits$dato[fits$dato>0],breaks=4))
  fits$catpred[fits$pred==0]=0
  if (max(fits$pred)>0) {fits$catpred[fits$pred>0]=as.numeric(cut(fits$dato[fits$pred>0],breaks=4))}
  fits$catfit=abs(fits$catdato-fits$catpred)
  fits$catfit[fits$catfit>1]=1
  centfluxfit0=1-sum(fits$catfit)/length(fits$catfit)
  centfluxfit1=1-sum(fits$catfit[fits$dato>0])/length(fits$catfit[fits$dato>0])
  # output
  spksall=c("SPKSall",lagcentfit0,lagfluxfit0,cpcentfit0,cpfluxfit0,centfluxfit0)
  spkspos=c("SPKSpos",lagcentfit1,lagfluxfit1,cpcentfit1,cpfluxfit1,centfluxfit1)
  
  ## return a compiled table with model coefficients and error stats
  return(rbind(coeff,errores,ks,spksall,spkspos))
  
}



# Función de evaluación de crames von mises
cvm_eval <- function(model, coords, nperm = 999){
  predic_model <- predict(model, model@frame[, -1]) %>% 
    abs()%>% 
    as.vector()
  model_cvm <- ecespa::syrjala0(coords = coords,R = F,
                        var1 = model@frame[, 1],
                        var2 = predic_model, nsim = nperm)
  return(model_cvm)
}



vals_jenk_breaks <-  function(x, deciles = 10, seed=1 ){
  set.seed(seed)
  
  y <- cut( x,c(unique(getJenksBreaks(x,deciles)), Inf),
            labels = F, include.lowest = T, right = F)
  return(y)
}


ppp_maker <- function(model, city, delito, ppp_predict = T, w, nb){
  deciles <- 10
  city_sf <- city %>% st_as_sf() 

    
  if(ppp_predict== T){
    if (class(model)[1] == "Sarlm"){
      ypred=as.numeric(predict(model, listw = nb, pred.type = "TS")**3)
    }
    if (class(model)[1] == "glmerMod"){
      ypred=exp(predict(model))
    }
    city <- city_sf%>% 
      mutate(ypred= ypred)%>% 
      mutate(ypred_adj = vals_jenk_breaks(ypred)) %>% 
      filter(ypred_adj > 1) %>% 
      as("Spatial")
    # marks = city@data[,"ypred_adj"]
  }else{
    city <- city_sf %>%
      mutate(delito_adj = vals_jenk_breaks(!!rlang::sym(delito))) %>%
      filter(delito_adj>1) %>% 
      as("Spatial")
    # marks = city@data[,"delito_adj"]
  }
  
  pts <- coordinates(city)
  p <- ppp(pts[,1], pts[,2], window = w)
  # ds_al <- density(p, adjust=.25)
  # plot(ds_al, main='Densidad de Delito')
  return(p)
  
}






predict_df  <- function(model, city, longlat = F){
  
  if (class(model)[1] == "Sarlm"){
    nb <- spatial_weights(city_df = city, nvec= 12,longlat = longlat)
    dato <- data.frame(ypred=as.numeric(
      predict(model, listw = nb, pred.type = "TS"))**3) %>% 
      mutate(delitos = model$y,
             ypred_adj = vals_jenk_breaks(ypred),
             y_adj = vals_jenk_breaks(delitos)) %>% 
      tibble::rownames_to_column(var = "rows") 
    
    
  }
  if (class(model)[1] == "glmerMod"){
    doble_max <-  (max(model@frame[, 1])*1.5)
    dato = data.frame(ypred = exp(predict(model))) %>% 
      mutate(ypred = ifelse(ypred > doble_max, doble_max, ypred)) %>% 
      mutate(delitos = model@frame[, 1],
             ypred_adj = vals_jenk_breaks(ypred),
             y_adj = vals_jenk_breaks(delitos)) %>% 
      tibble::rownames_to_column(var = "rows") 
    
  }
  
  base <- city@coords %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "rows") %>% 
    left_join(dato, by = "rows") %>% 
    rename(x=2, y=3)
  
  
  return(base)
}


cvm_eval_spatial <- function(coords_x, coords_y, 
                             var1, var2, nperm=99){
  
  coords_df <- data.frame(x= coords_x, y =coords_y)
  cvm <- ecespa::syrjala(coords = coords_df, 
                         var1 = var1,
                         var2 = var2, 
                         nperm = nperm)
  return(cvm$cvm.sim)
}


rasterizar_density <- function(x, s, crs_used){
  x <- raster(x)
  proj4string(x)=crs_used
  sample <-raster(extent(s), nrow = nrow(s@data), 
                  ncol = ncol(s@data), crs = crs_used)
  res_x <- resample(x,sample)
  return(res_x)
}


rescale01 <- function(x) {
  val <- values(x)
  values(x) <- (val - min(val, na.rm = TRUE)) / (max(val, na.rm = TRUE) - min(val, na.rm = TRUE))
  x
}

