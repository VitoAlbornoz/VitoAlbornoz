rm(list=ls())

library(RStoolbox)
library(stringr)
library(raster)
library(sp)
library(rgdal)
library(snow)
library(landsat)
library(gtools)
library(beepr)

'%ni%'<-Negate('%in%')
DIRECCION <- "E:/FONDECYT/benjamin"
ade18s<-readOGR("E:/FONDECYT/benjamin/Landsat8/adecuad18s.shp")
ade19s<-readOGR("E:/FONDECYT/benjamin/Landsat8/adecuad19s.shp")
aderas<-raster("E:/FONDECYT/benjamin/aderas.tif")

dir<-file.path(DIRECCION,"escenas/sr")
dirfinal<-file.path(DIRECCION,"escenas/sr/TC")
setwd(dir)

carpetas<-list.files(dir)
carpetas<-carpetas[str_detect(carpetas,".tar") & str_detect(carpetas,"20010328")]# & str_detect(carpetas,"091")]# |
                   #str_detect(carpetas,".tar") & str_detect(carpetas,"2000") ]

azi<-vector()
ele<-vector()

mtl<-list()

sombra<-c("sombra_solo091_modificado.shp","sombra_solo092_modificado.shp")

#ciclo para hacer la corrección topográfica de las escenas, el primero (l) es para ir de escena en escena, el segundo (b) es para la correccion banda a banda
for(l in 1:length(carpetas)){
  # lectura del metadata y demas 
  
  {
    sce<-paste(
      #substr(carpetas[l],1,4),
      substr(carpetas[l],5,10),
      paste(substr(carpetas[l],11,14),
            substr(carpetas[l],15,16),
            substr(carpetas[l],17,18),sep="-"),
      sep="_")
      setwd(dir)
    esc<-file.path(getwd(),sce)
    
    if(dir.exists(esc)==T){
      sensor  <- substr(carpetas[l],1,4)
      pathrow <- substr(carpetas[l],5,10)
      setwd(esc)                     #definir el directorio de trabajo como la carpeta recién creada
      
    } else {
      untar(carpetas[l],exdir=esc)
      sensor  <- substr(carpetas[l],1,4)
      pathrow <- substr(carpetas[l],5,10)
      setwd(esc)                     #definir el directorio de trabajo como la carpeta recién creada
    }
    
    mtl[[l]]<-readMeta(list.files()[str_detect(list.files(),"MTL")==TRUE])
    
    if(file.exists(file.path(esc, "topocor",paste("TC_",sensor,"_",sce,".tif",sep="")))){
      print(paste("escena",sce,"ya procesada, continúa el análisis / hora =",Sys.time()))
      next()} else {print(paste("comienza escena",sce,"a las",Sys.time()))}
    
    azi[l]<- mtl[[l]]$SOLAR_PARAMETERS[1]*pi/180
    ele[l]<- mtl[[l]]$SOLAR_PARAMETERS[2]*pi/180
  }

  #corrección del extent de las escenas landsat
  {
    bandas<-list.files(pattern="_band")
    ex_lt<-extent(raster(bandas[1]))
    exok<-extent(ex_lt[1],ex_lt[2],ex_lt[3]+10000000,ex_lt[4]+10000000)
    ban<-stack(bandas)
    ban<-setExtent(ban,exok)
    
    if(str_detect(crs(ban),"zone=18")==T){
      print(paste("comienza proceso de arreglo de extent para path 232 a las",Sys.time()))
           proj4string(ban)<-"+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
           ban<-projectRaster(ban,aderas)*0.0001
           beep(sound=2)
           print("termina reproyección")
           #ban<-extend(ban,ade19s)
           #ban<-crop(ban,ade19s)*0.0001
           
           writeRaster(ban,"pretopo.tif",format="GTiff",datatype="FLT4S",overwrite=T)
           print(paste("termina proceso de arreglo de extent para path 232 a las",Sys.time()))
           beep(sound=2)
           
           } else {
      proj4string(ban)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
      ban<-crop(ban,ade19s)*0.0001
      writeRaster(ban,"pretopo.tif",format="GTiff",datatype="FLT4S",overwrite=T)
      }
  }
  
  # correccion radiometrica 
  {
    #print("comienza correccion radiometrica")
    #haze<-estimateHaze(ban,hazeBands=1:4,plot=F)
    #bsr<-radCor(ban,metaData=mtl[[l]],method="sdos",hazeValues=haze,hazeBands = 1:4)
  }
  
  # enmascarado de nubes
  {
  #  if(file.exists(paste(sce,"mascara.tif",sep="_"))){
  #    mascara<-raster(paste(sce,"mascara.tif",sep="_"))} else {
  #      qa     <- raster(list.files(esc,pattern="pixel_qa"))
  #      radsat <- raster(list.files(esc,pattern="radsat_qa"))
  #      qa<-setExtent(qa,exok)
  #      radsat<-setExtent(radsat,exok)
  #      if(sensor%in%c("LT04","LT05","LE07")){
  #        qa[qa%ni%c(72,136,80,112,144,176,96,112,160,176,224,130,132,136,144,160,176,224)]<-1
  #        qa[qa%in%c(72,136,80,112,144,176,96,112,160,176,224,130,132,136,144,160,176,224)]<-NA
  #        
  #        radsat[radsat%ni%c(8,16)]<-1
  #        radsat[radsat%in%c(8,16)]<-NA
  #        
  #        st<-stack(qa,radsat)
  #        
  #        if(str_detect(crs(st),"zone=18")==T){
  #           print(paste("comienza proceso de arreglo de extent para path 232 a las",Sys.time()))
  #           proj4string(st)<-"+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #           st<-projectRaster(st,aderas)
  #           beep(sound=2)
  #           print("termina reproyección")
  #      } else {
  #           proj4string(st)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #           st<-crop(st,ade19s)
  #      }
  #        
  #        mascara<-st[[1]]*st[[2]]
  #        
  #        writeRaster(st[[1]],paste(sce,"_qa",".tif",sep=""),
  #                   format="GTiff",datatype="FLT4S",overwrite=T)
  #        writeRaster(st[[2]],paste(sce,"_radsat",".tif",sep=""),
  #                   format="GTiff",datatype="FLT4S",overwrite=T)
  #        writeRaster(mascara,paste(sce,"_mascara",".tif",sep=""),
  #                   format="GTiff",datatype="FLT4S",overwrite=T)} else 
  #                     if (sensor=="LC08") {
  #          aerosol<- raster(list.files(esc,pattern="sr_aerosol"))
  #          proj4string(aerosol)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #          aerosol<- setExtent(aerosol,exok)
  #          
  #          qa[qa%ni%c(328, 392, 840, 904, 1350,336, 368, 400, 432, 848, 880, 912, 944, 1352,480, 992,834, 836, 840, 848, 864, 880, 898, 900, 904, 912, 928, 944, 992)]<-1
  #          qa[qa%in%c(328, 392, 840, 904, 1350,336, 368, 400, 432, 848, 880, 912, 944, 1352,480, 992,834, 836, 840, 848, 864, 880, 898, 900, 904, 912, 928, 944, 992)]<-NA
  #          
  #          radsat[radsat%ni%c(16,32,48)]<-1
  #          radsat[radsat%in%c(16,32,48)]<-NA
  #          
  #          aerosol[aerosol%ni%c(194, 196, 200, 208, 224, 228)]<-1
  #          aerosol[aerosol%in%c(194, 196, 200, 208, 224, 228)]<-NA
  #          
  #          st<-stack(qa,radsat,aerosol)
  #          
  #          if(str_detect(crs(st),"zone=18")==T){
  #            print(paste("comienza proceso de arreglo de extent para path 232 a las",Sys.time()))
  #            proj4string(st)<-"+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #            st<-projectRaster(st,aderas)
  #            beep(sound=2)
  #            print("termina reproyección")
  #            } else {
  #              proj4string(st)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #              st<-crop(st,ade19s)
  #          }
  #          
  #          mascara<-st[[1]]*st[[2]]*st[[3]]
  #          
  #          writeRaster(st[[1]],paste(esc,sce,"_qa",".tif",sep=""),
  #                      format="GTiff",datatype="FLT4S",overwrite=T)
  #          writeRaster(st[[2]],paste(esc,sce,"_radsat",".tif",sep=""),
  #                      format="GTiff",datatype="FLT4S",overwrite=T)
  #          writeRaster(st[[3]],paste(esc,sce,"_aerosol",".tif",sep=""),
  #                      format="GTiff",datatype="FLT4S",overwrite=T)
  #          
  #          writeRaster(mascara,paste(esc,sce,"_mascara",".tif",sep=""),
  #                      format="GTiff",datatype="FLT4S",overwrite=T)
  #        }
  #  aerosol<- NULL
  #  radsat <- NULL
  #  qa     <- NULL
  #  st     <- NULL
  #  #ban<-ban*mascara
  #  #mascara<-NULL
  #  }
  }
      
  # correccion topografica
  {
    print("inicia correccion topográfica")
    #ilu<-raster("illu.tif")
    if(any(str_detect(list.files(),"illu.tif")==T)) {
      ilu<-raster("illu.tif")
      print("raster de iluminacion ya existe, continúa el análisis")
      } else {
      print("se creará nuevo raster de iluminación")
      dem<-crop(raster("E:/FONDECYT/benjamin/Landsat8/demangela_19s.tif"),extent(ban))
      #dem<-crop(raster("E:/FONDECYT/benjamin/dem/dem_prefill_19s.tif"),ade19s)
      ilu<-topCor(ban,dem,metaData=mtl[[l]],method="illu")
      writeRaster(ilu,filename="illu.tif",format="GTiff",overwrite=T)
      rm(dem)
      }
    
    if(str_detect(esc,"091")==T){
      sombras<-readOGR(file.path(DIRECCION,"escenas",sombra[1]))} else {sombras<-readOGR(file.path(DIRECCION,"escenas",sombra[2]))}
      
    datosreg<-raster::extract(stack(ban,ilu),sombras,df=T)[,-1]
    colnames(datosreg)<-c(paste("ref_B",c(1:5,7),sep=""),"ilu")
    reg<-list()
    
    pdf(file=paste("ref_vs_ilu_",sce,".pdf",sep=""))  
    par(mfrow=c(3,2))
    for (h in 1:dim(ban)[3]){
      reg[[h]]<-lm(datosreg[,h]~ilu,data=datosreg)
      plot(datosreg[,h]~datosreg[,dim(ban)[3]+1],xlab="ilu",ylab=paste("banda",h),main=paste("banda",h),cex.main=0.75,cex.lab=0.75,pch=16)
      abline(reg=reg[[h]],lwd=3,col="red")
      legend("topleft", legend=c(paste("R2 adj =",round(as.numeric(summary(reg[[h]])[9]),4),sep=""),
             paste("y = ",round(reg[[h]]$coefficients[2],4),"x + ",round(reg[[h]]$coefficients[1],4),sep="")),
             cex=0.7,bty="n",fill="white")
    }
    dev.off()
    
    topocor<-list()
    c<-vector()
    pdf(paste("pretopo_vs_postopo",sce,".pdf",sep=""))
    par(mfrow=c(3,3))
    for(b in 1:length(bandas)){
      print(paste("comienza banda",b,"a las",Sys.time()))
      c[b]<-reg[[b]]$coefficients[1]/reg[[b]]$coefficients[2]
      
      topocor[[b]]<-ban[[b]]*( ( cos(pi/2-ele[l]) + c[b]) / (ilu+c[b]))
      # enmascarado de la escena en zonas vegetadas y no vegetadas #####
      #b1 <- mask(ban,m1)#,updatevalue=NA)
      #b2 <- mask(ban,m2)#,updatevalue=)
      
      #tccomp<-topCor(ban,dem=dem,metaData=mtl[[l]],method="C")#solarAngles = c(azi[l],ele[l]),method="C")#,illu=ilu)
      
      #b1.tp <- topCor(b1,dem=dem,metaData=mtl [[l]] ,method="C", illu = raster("illu.tif"))
      #b2.tp <- topCor(b2,dem=dem,metaData=mtl [[l]] ,method="C", illu = raster("illu.tif"))
      
      #b1.tp <- topocorr(b1,slope=sloasp[[1]],aspect=sloasp[[2]],sunelev=ele[l],sunazimuth = azi[l],method="ccorrection")
      #b2.tp <- topocorr(b2,slope=sloasp[[1]],aspect=sloasp[[2]],sunelev=ele[l],sunazimuth = azi[l],method="ccorrection")
      # unión de ambas escenas ya corregidas
      #b1.tp[is.na(b1.tp)] <- 0
      #b2.tp[is.na(b2.tp)] <- 0
      
      #tc[[b]] <-b1.tp+b2.tp
      #####
      
      #writeRaster(topocor[[b]],filename=paste("topocor/TC_B",b,sep=""),format="GTiff",overwrite=T,datatype = 'FLT4S')
      plot(ban[[b]],topocor[[b]],xlab="pretopo",ylab="correc topo",main=paste("banda",b),ylim=c(0,1),xlim=c(0,1),cex.main=0.75,cex.lab=0.75)
      print(paste("termina banda",b))#,"a las",Sys.time()))
    }
    
    dev.off()
    print("termina corrección topográfica")
  }
  
  # escritura de los datos corregidos
  {
    if(!dir.exists(file.path(esc, "topocor"))) {dir.create(file.path(esc, "topocor"))} else {setwd(file.path(esc, "topocor"))}
     
    setwd(file.path(esc, "topocor"))
    topocor<-stack(topocor)
    writeRaster(topocor,filename=paste("TC_",sensor,"_",sce,".tif",sep=""),
                format="GTiff",overwrite=T,datatype = 'FLT4S',bylayer=F)
  }
  
  rm(topocor)
  rm(ban)
  rm(ilu)
  rm(reg)
  rm(sombras)
  rm(datosreg)
   
  print(paste("termina escena",paste(substr(carpetas[l],5,10),substr(carpetas[l],11,14)),"a las",Sys.time()))
}




# pedazo parcialmente inutilde la correccion topografica ###### 
#calculo del NDVI y asignación del extent correcto
#ndvi<-(stack(bandas[5])-stack(bandas[4]))/(stack(bandas[5])+stack(bandas[4]))
#ndvi<-setExtent(ndvi,exok)
#proj4string(ndvi)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#ndvi<-crop(ndvi,ade19s)
#m1 <- reclassify(ndvi,rcl=matrix(c(-Inf,0.5,0.5,Inf,1,NA),ncol=3),include.lowest=T,right=F)
#m2 <- reclassify(ndvi,rcl=matrix(c(-Inf,0.5,0.5,Inf,NA,1),ncol=3),include.lowest=T,right=F)
#rm(ndvi)

# ciclo de corrección topográfica discriminando por cobertura

#dem <- crop(raster("E:/benjamin/landsat8/demangela_19s.tif"),ade19s)
#dem <- resample(dem,m1,method="bilinear")
#sloasp<-terrain(dem,opt=c("slope","aspect"),neighbors = 8,unit="degrees")
#ilu<-topCor(img=dem,dem=dem,method="illu",metaData=mtl[[l]])
#writeRaster(ilu,filename="illu.tif",format="GTiff",datatype="FLT4S",overwrite=T)
#####
  
