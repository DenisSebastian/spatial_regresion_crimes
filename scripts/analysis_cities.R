####################################
###### Crime models comparison #####
####################################

####### 1. Setup #######

options(scipen = 999, stringsAsFactors=F, encoding = "UTF-8",max.print=1000000)

#### load libraries
library(DescTools)
library(spatialreg)
library(MASS)
library(lme4)
library(doSNOW)
library(snow)
library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(xlsx)
library(clusterSim)
library(spdep)

# setwd("C:/Users/matia/Documents/Papers/crime_cellphones")

#### List of crimes (delitos) and cities (ciudades)
crimes=c("Minor property offenses","Robberies","Burglaries","Drunkennes, damages and disorders",
         "Family violence and aggression","Injuries, drugs and weapons","High violence and murder")
delitos=c("hurto","robo","asalto","amen","abus","crim","viol")

ciudades=gsub(".rds","",list.files("data",pattern="urb"))
ciudades=ciudades[ciudades=="urb_4_18"]
cities=ciudades

# cities=c("Antofagasta","Coquimbo-Serena","Temuco")
# ciudades=c("antofa","serena","temuco")
# ciudades="serena"

#### Results cleaning functions ######

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
  return(rbind(coeff,errores,ks,cvm, spksall,spkspos))

}

###################################

####### 2. Crime Analysis and Models by city #######

### Perform regressions for each city

####  Parallelization parameters
mcoptions=list(preschedule=FALSE, set.seed=FALSE)
ncores=2
cl = makeCluster(ncores,type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
iter=length(ciudades)

modelos <- foreach(obs=1:iter,.options.multicore=mcoptions,.packages=c("spatialreg","MASS","lme4","spdep","clusterSim")) %dopar% {
  # obs = 1
  city=ciudades[obs]
  ## Read city's spatialpointsdataframe
  cityp= readRDS(paste0("data/",city,".rds"))
  ## create spatial weigths matrix with 12 neighbors
  nvec = 12
  nb = nb2listw(neighbours=knn2nb(knn=knearneigh(cityp,nvec)),style = "W")

  ## generate spatially lagged variables
  norm = "n4"
  cityp$lagviol   <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$viol),type=norm))
  cityp$lagcrim   <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$crim),type=norm))
  cityp$lagabus   <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$abus),type=norm))
  cityp$lagamen   <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$amen),type=norm))
  cityp$lagasalto <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$asalto),type=norm))
  cityp$lagrobo   <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$robo),type=norm))
  cityp$laghurto  <- as.numeric(data.Normalization(lag.listw(x=nb, var=cityp$hurto),type=norm))

  ## Normalize independent variables
  cityp$tamhog = as.numeric(data.Normalization(cityp$tamhog,norm))
  cityp$neduc  = as.numeric(data.Normalization(cityp$educ,norm))
  cityp$vulzon = as.numeric(data.Normalization(cityp$vulzon,norm))
  cityp$mobil  = as.numeric(data.Normalization(cityp$mobil,norm))

for (delito in delitos){

    # ## Lagged crime and Centrality
    # lagcent=glmer.nb( formula(paste0(delito,"~lag",delito,"+lofcom+(1|zona)")),data=cityp@data,verbose=F)
    # saveRDS(lagcent,paste0("output/",city,"_lgct_",delito,".rds"))
    # 
    # ## Lagged crime and Flux
    # lagflux=glmer.nb( formula(paste0(delito,"~lag",delito,"+lpflot+(1|zona)")),data=cityp@data,verbose=F)
    # saveRDS(lagflux,paste0("output/",city,"_lgfl_",delito,".rds"))
    # 
    # ## Centrality with controls
    # cpcent=glmer.nb( formula(paste0(delito,"~lag",delito,"+lofcom+ldenspob+tamhog+neduc+vulzon+mobil+(1|zona)")),data=cityp@data,verbose=F)
    # saveRDS(cpcent,paste0("output/",city,"_ctcn_",delito,".rds"))
    # 
    # ## Flux with controls
    # cpflux=glmer.nb( formula(paste0(delito,"~lag",delito,"+lpflot+ldenspob+tamhog+neduc+vulzon+mobil+(1|zona)")),data=cityp@data,verbose=F)
    # saveRDS(cpflux,paste0("output/",city,"_flcn_",delito,".rds"))
    # 
    # ## Flux and Centrality
    # centflux=glmer.nb( formula(paste0(delito,"~lag",delito,"+lofcom+lpflot+ldenspob+tamhog+neduc+vulzon+mobil+(1|zona)")),data=cityp@data,verbose=F)
    # saveRDS(centflux,paste0("output/",city,"_ctfl_",delito,".rds"))
    
    

# lagsarlm ----------------------------------------------------------------

    
    ## Lagged crime and Centrality
    lagsar_cent=lagsarlm( formula(paste0(delito,"~lag",delito,"+lofcom")),
                          data=cityp@data, listw = nb, type="lag",method="eigen")
    
    saveRDS(lagsar_cent,paste0("output/",city,"_lagsarct_",delito,".rds"))
    
    ## Lagged crime and Flux
    lagsar_lux=lagsarlm( formula(paste0(delito,"~lag",delito,"+lpflot+zona")),data=cityp@data, listw = nb)
    saveRDS(lagsar_lux,paste0("output/",city,"_lagsarfl_",delito,".rds"))
    
    ## Centrality with controls
    lagsar_cpcent=lagsarlm( formula(paste0(delito,"~lag",delito,"+lofcom+ldenspob+tamhog+neduc+vulzon+mobil+zona")),data=cityp@data,listw = nb)
    saveRDS(lagsar_cpcent,paste0("output/",city,"_lagsar_ctcn_",delito,".rds"))
    
    ## Flux with controls
    lagsar_cpflux=lagsarlm( formula(paste0(delito,"~lag",delito,"+lpflot+ldenspob+tamhog+neduc+vulzon+mobil+zona")),data=cityp@data,listw = nb)
    saveRDS(lagsar_cpflux,paste0("output/",city,"lagsar_flcn_",delito,".rds"))
    
    ## Flux and Centrality
    lagsar_centflux=lagsarlm( formula(paste0(delito,"~lag",delito,"+lofcom+lpflot+ldenspob+tamhog+neduc+vulzon+mobil+zona")),data=cityp@data,listw = nb)
    saveRDS(lagsar_centflux,paste0("output/",city,"lagsar_ctfl_",delito,".rds"))

  }

  ## End parallel
  stopCluster(cl)
}

###################################

####### 3. Results compilation #######

#### Creation of Data structures for results gathering
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

resultados_test <- data.frame()
## for every city
for (ciud in 1:length(ciudades)) {
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

  ## for every crime
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

    ## Tabulate regression
    reg=tabula_coef(cityp,lagcent,lagflux,cpcent,cpflux,centflux)
    reg=rbind(rep(city,5),rep(crime,5),names(reg),reg)
    reg=reg[!substr(reg$Vars, 1, 4)%in%c("fitt","devi","df.r"),]
    reg$Vars=c("City","Crime Category","Model","(Intercept)","Lagged Crime",
               "Log(Flux Pop)","Log(Dens Pop)","Log(Dens Facilities)",
               "Household Size","Education","Vulnerability (Zone)","New residents (Zone)",
               "AIC","BIC","logLik","KS statistic","CVM statistic",  "PSpatial KS all","PSpatial Ks pos")
    regs=rbind(regs,reg)

    ## Select best models
    bestpos=as.numeric(reg[nrow(reg),2])
    bestall=as.numeric(reg[(nrow(reg)-1),2])
    bestcvm=as.numeric(reg[nrow(reg)-2,2])
    
#tablas de todos los test
    res <- reg[(nrow(reg)-3):nrow(reg), ] %>% 
      t() %>% 
      `colnames<-`(.[1, ]) %>%
      .[-1, ] %>% 
      as.data.frame() %>% 
      mutate_if(is.character, as.numeric) %>% 
      mutate(Ciudad = ciudad, Delito = delito) %>% 
      tibble::rownames_to_column("Tipo")
    
    resultados_test <-  resultados_test %>%
      rbind(res)
    
    
    modbest=modsimp="Center"
    if(as.numeric(reg[(nrow(reg)-1),3])>bestall) {
      bestpos=as.numeric(reg[nrow(reg),3])
      bestall=as.numeric(reg[(nrow(reg)-1),3])
      modbest=modsimp="Flux"
    }
    modcont="Center cont"
    if(as.numeric(reg[(nrow(reg)-1),4])>bestall) {
      bestpos=as.numeric(reg[nrow(reg),4])
      bestall=as.numeric(reg[(nrow(reg)-1),4])
      modbest="Center cont"
    }
    if(as.numeric(reg[(nrow(reg)-1),5])>bestall) {
      bestpos=as.numeric(reg[nrow(reg),5])
      bestall=as.numeric(reg[(nrow(reg)-1),5])
      modbest=modcont="Flux cont"
    }
    if(as.numeric(reg[(nrow(reg)-1),6])>bestall) {
      bestpos=as.numeric(reg[nrow(reg),6])
      bestall=as.numeric(reg[(nrow(reg)-1),6])
      modbest="Cent Flux"
    }

    citystmod=c(cityst,crime,as.integer(citystat$n[10-del]),modsimp,modcont,modbest,round(bestall,3),round(bestpos,3))
    syntstat=rbind(syntstat,citystmod)

  }
}

saveRDS(resultados_test, "output/test_tbl.rds")

stats=stats[match(c("1","2","hurto","robo","asalto","amen","abus","crim","viol",
                    "flux","pobdens","ofcom","tamhog","educ","vuln","mobil")
                  ,stats$Variables),]
stats$Variables=c("City","Statistic","Minor property offenses","Robberies","Burglaries","Drunkennes, damages and disorders",
                  "Family violence and aggression","Injuries, drugs and weapons","High violence and murder",
                  "Flux Population","Density Population","Density Facilities","Household size","Education Average",
                  "Vulnerability (Zones)","New residents (Zones)")
names(syntstat)=c("City","Inhabitants","Commerce_Offices","Flux","Crime_type","N_Crimes",
                  "Simple_model","Controls_model","Best_model","PSpatial KS all","PSpatial Ks pos")


improve=data.frame()
for (crime in crimes){
  imp=length(syntstat$City[syntstat$Crime_type==crime & syntstat$Best_model %in% c("Flux","Flux cont","Cent Flux")])/
  length(syntstat$City[syntstat$Crime_type==crime])  
  imp=c(crime,imp)
  improve=rbind(improve,imp)
}
crime="All"
imp=length(syntstat$City[syntstat$Best_model %in% c("Flux","Flux cont","Cent Flux")])/
  length(syntstat$City)  
imp=c(crime,imp)
improve=rbind(improve,imp)
names(improve)=c("Crime_Type","Improvement")

synt=merge(as.data.frame(table(syntstat$Simple_model)),as.data.frame(table(syntstat$Controls_model)),by="Var1",all.x=T,all.y=T)
synt=merge(synt,table(syntstat$Best_model),by="Var1",all.x=T,all.y=T)

names(synt)=c("Models","Simple","Controls","Best")

###################################

####### 4. Export Tables #######
write.xlsx(regs, file="output/results.xlsx", sheetName="Regressions", row.names=F, col.names=T)
write.xlsx(stats, file="output/results.xlsx", sheetName="Statistics", append=T, row.names=F, col.names=T)
write.xlsx(syntstat, file="output/results.xlsx", sheetName="Summary", append=T, row.names=F, col.names=T)
write.xlsx(synt, file="output/results.xlsx", sheetName="Best_models", append=T, row.names=F, col.names=T)
write.xlsx(improve, file="output/results.xlsx", sheetName="Improvements", append=T, row.names=F, col.names=T)
###################################

saveRDS(resultados_test, "output/test_tbl.rds")


