###==== reclasificacion por transiciones ====###
rm(list=ls())

library(raster)
library(rgdal)
library(beepr)

ruta<-"E:/FONDECYT/benjamin/clasificaciones_final"
setwd(ruta)

c1984<- raster("c1984_pretrans_maj.tif")
c2000<- raster("c2000_pretrans_maj.tif")
c2018<- raster("c2018_pretrans_maj.tif")

clasificaciones<-list(
  c1984<-raster("E:/FONDECYT/benjamin/clasificaciones_final/clas1984_feb2019_pretransicion.tif"),
  c2000<-raster("E:/FONDECYT/benjamin/clasificaciones_final/clas2000_pretransicion.tif"),
  c2018<-raster("E:/FONDECYT/benjamin/clasificaciones_final/clas2018_feb2019_pretransicion.tif"))

dem  <- crop(raster("E:/FONDECYT/benjamin/Landsat8/demangela_19s.tif"),c1984)
ade  <- raster("E:/FONDECYT/benjamin/aderas.tif")


'%ni%'<-Negate('%in%')


# transicion de praderas, pastoreo y matorral abierto a primario 
{
  h1000<-reclassify(dem,matrix(c(-Inf,1000,1000,Inf,0,1),ncol=3))
  #extension_primario84<-rasterize(readOGR("mascaras_transiciones/extensiones/extension_zonas_primario84.shp"),c1984,field="Id",background=0,
  #                                  filename="mascaras_transiciones/extensiones/extension_zonas_primario84.tif",format="GTiff",datatype="INT1U",overwrite=T)
  extension_primario84   <-raster("mascaras_transiciones/extensiones/extension_zonas_primario84.tif")
  
  #extension_secundario18<-rasterize(readOGR("mascaras_transiciones/extensiones/extension_zonas_secundario18.shp"),c1984,field="Id",background=0,
  #                                    filename="mascaras_transiciones/extensiones/extension_zonas_secundario18.tif",format="GTiff",datatype="INT1U",overwrite=T)
  extension_secundario18 <-raster("mascaras_transiciones/extensiones/extension_zonas_secundario18.tif")
  
  pr_pa_ma_84<-c1984
  pr_pa_ma_84<-reclassify(pr_pa_ma_84,matrix(c(7,8,11,12,15,16,17,1,2,3,4,5,6,9,10,13,14,18,19,20,1,1,1,1,1,1,1,rep(0,13)),ncol=2))
  #pr_pa_ma_84[pr_pa_ma_84%in%c(7,8,11,12)]<-1
  
  primario18<-c2018
  primario18<-reclassify(primario18,matrix(c(1:12,13,14:20,rep(0,12),1,rep(0,7)),ncol=2))
  #primario18[primario18==13]<-1
  
  correccion_prads_primario<-pr_pa_ma_84*primario18
  
  limite_cprim84<-mask(h1000+extension_primario84,extension_secundario18,maskvalue=1,updatevalue=0);limite_cprim84[limite_cprim84==2]<-1
  limite_csec18<-mask(reclassify(h1000,matrix(c(1,0,0,1),ncol=2))+extension_secundario18,extension_primario84,maskvalue=1,updatevalue=0);limite_csec18[limite_csec18==2]<-1
  
  c_prim84<-correccion_prads_primario*limite_cprim84
  c_sec18<-correccion_prads_primario*limite_csec18
  
  
  plot(c_prim84)
  plot(c_sec18)
  
  writeRaster(c_prim84*ade,"mascaras_transiciones/correccion_praderas84_primario84.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  writeRaster(c_sec18*ade,"mascaras_transiciones/correccion_primario18_secundario18.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  #rm(pr_pa_ma_84)
  #rm(primario18)
  #rm(correccion_prads_primario)
  beep()
}
# transicion de secundario a primario dentro de la máscara de Nahuel
{
  nubes<-raster("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zona_nubes_nahuel.tif")
  secundario84<-c1984*nubes
  primario84<-c1984*nubes
  pastoreo84<-c1984*raster("E:/FONDECYT/benjamin/clasificacion1984/febrero2019/zonas/zonas_secundario_pastoreo_nahuel.tif")
  
  primario00<-c2000*nubes
  sd00<-c2000*nubes
  veg_altura00<-c2000*nubes
  
  
  
  pastoreo84[pastoreo84!=14 | is.na(pastoreo84)]<-0
  pastoreo84[pastoreo84==14]<-1
  
  
  secundario84<-reclassify(secundario84,
                           matrix(c(1:13,14,15:20,rep(0,13),1,rep(0,6)),ncol=2))
  #secundario84[secundario84==14]<-1
  primario84<-reclassify(primario84,matrix(c(1,2:12,13,14:20,1,rep(0,11),1,rep(0,7)),ncol=2))
  #primario84[primario84==13]<-1
  
  primario00<-reclassify(primario00,matrix(c(1,2:12,13,14:20,1,rep(0,11),1,rep(0,7)),ncol=2))
  #primario00[primario00==13]<-1
  sd00<-reclassify(sd00,matrix(c(1:14,15,16,17:20,rep(0,14),1,1,1,rep(0,3)),ncol=2))
  #sd00[sd00%in%c(15,16)]<-1
  veg_altura00<-reclassify(veg_altura00,matrix(c(1:16,17,18:20,rep(0,16),1,rep(0,3)),ncol=2))
  
  correccion_secundario_primario<-secundario84*primario00
  correccion_primario_sd<-primario84*sd00
  plot(correccion_secundario_primario)
  plot(correccion_primario_sd)
  plot(veg_altura00)
  
  writeRaster(correccion_secundario_primario,"mascaras_transiciones/correccion_secundario_primario_nahuel.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  writeRaster(correccion_primario_sd,"mascaras_transiciones/correccion_primario_sd_nahuel.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  writeRaster(veg_altura00,"mascaras_transiciones/correccion_veg_altura_nahuel.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  writeRaster(pastoreo84,"mascaras_transiciones/correccion_pastoreo_nahuel.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  
  rm(nubes)
  rm(secundario84)
  rm(primario18)
  rm(correccion_secundario_primario)
  rm(correccion_primario_sd)
  rm(sd18)
  rm(primario84)

}

# zona de secundarios a pastoreo en la máscara de nahuel


# transicion de matorral arborescente a bosque primario en la zona de Tapera
{
  #zona_tapera<-rasterize(readOGR("zona_correccion_tapera.shp"),c1984,field="Name",background=0,
  #                       filename="zona_correccion_tapera.tif",format="GTiff",datatype="INT1U",overwrite=T)
  zona_correccion<-raster("zona_correccion_tapera.tif")+raster("E:/FONDECYT/benjamin/congresoCIEP/limite_pampa_bosque.tif")
  
  ma84<-c1984*zona_correccion
  ma84[ma84!=9]<-0
  ma84[ma84==9]<-1
  
  prim18<-c2018*zona_correccion
  prim18[prim18!=13]<-0
  prim18[prim18==13]<-1
  
  prim_sec18<-ma84*prim18
  plot(prim_sec18)
  
  writeRaster(prim_sec18,"mascaras_transiciones/primario18_secundario18_pampas.tif",
              format="GTiff",dfatatype="INT1U",overwrite=T)
  
}
beep(sound=3)

# ================================================================= # 
# ========================== ENMASCARADO ========================== # 
# ================================================================= #

mascaras<-file.path(getwd(),"mascaras_transiciones",list.files("mascaras_transiciones",pattern="nahuel.tif$"))

mascaras<-lapply(lapply(mascaras,raster),function(x){reclassify(x,matrix(c(1,0,0,1),ncol=2))})

#sd1984<-c1984
#sd1984[sd1984%ni%15]<-0
#sd1984[sd1984%in%15]<-1
#sd1984<-sd1984*reclassify(dem,matrix(c(-Inf,1200,1200,Inf,0,1),ncol=3))

#sd1984<-reclassify(sd1984,matrix(c(0,1,1,0),ncol=2))
#c1984<-c1984*sd1984+reclassify(sd1984,matrix(c(1,0,0,16),ncol=2))

# correccion de zonas de praderas en altura en 1984 a primario en 1984 por efecto de nieve y demás
#c1984<-c1984*mascaras[[1]]
#c1984<-c1984+reclassify(mascaras[[1]],matrix(c(1,0,0,13),ncol=2))

# 1984= corrección de secundario a pastoreo en los bosquecitos de la pampa en coyhaique alto
c1984<-c1984*mascaras[[1]]
c1984<-c1984+reclassify(mascaras[[1]],matrix(c(1,0,0,12),ncol=2))

# 1984= correccion de primarios a suelo desnudo en la zona de nubes a partir de la clasificación de Nahuel
#{c1984<-mask(c1984,raster(mascaras[[2]],maskvalue=1,updatevalue=16))}
c1984<-c1984*mascaras[[2]]
c1984<-c1984+reclassify(mascaras[[2]],matrix(c(1,0,0,16),ncol=2))

# correccion de zonas de primario en 2018 a secundario en 2018 por efecto de rebrote mal clasificado
#c2018<-c2018*mascaras[[3]]
#c2018<-c2018+reclassify(mascaras[[3]],matrix(c(1,0,0,14),ncol=2))

# 1984= correccion de bosques secundarios a primarios en la zona de nubes a partir de la clasificación de Nahuel
#{c1984<-mask(c1984,raster(mascaras[[4]],maskvalue=1,updatevalue=13))}
c1984<-c1984*mascaras[[3]]
c1984<-c1984+reclassify(mascaras[[3]],matrix(c(1,0,0,13),ncol=2))

# 1984= corrección de vegetación en altura en zona de nubes a partir de clasificación de Nahuel
#{c1984<-mask(c1984,raster(mascaras[[5]],maskvalue=1,updatevalue=17))}
c1984<-c1984*mascaras[[4]]
c1984<-c1984+reclassify(mascaras[[4]],matrix(c(1,0,0,17),ncol=2))

# 2018= corrección de bosques primarios a bosques secundarios en la pampa 
#c2018<-c2018*mascaras[[6]]
#c2018<-c2018+reclassify(mascaras[[6]],matrix(c(1,0,0,14),ncol=2))



noclas84<-c1984
noclas84[noclas84!=20]<-0
noclas84[noclas84==20]<-1

noclas18<-c2018
noclas18[noclas18!=20]<-0
noclas18[noclas18==20]<-1

par(mfrow=c(1,2));plot(noclas84);plot(noclas18)

noclas<-noclas84+noclas18
noclas[noclas==2]<-1
plot(noclas)

c1984<-mask(c1984,noclas,maskvalue=1,updatevalue=20)
c2018<-mask(c2018,noclas,maskvalue=1,updatevalue=20)

writeRaster(noclas,"mascara_noclas.tif",
            format="GTiff",datatype="INT1U",overwrite=T)

writeRaster(c1984,"c1984_para_transicion.rst",format="IDRISI",datatype="INT1U",overwrite=T)
writeRaster(c1984,"clas1984_mascaranahuel.tif",format="GTiff",datatype="INT1U",overwrite=T)

writeRaster(c2018,"c2018_para_transicion.rst",format="IDRISI",datatype="INT1U",overwrite=T)
writeRaster(c2018,"c2018_para_transicion.tif",format="GTiff",datatype="INT1U",overwrite=T)
beep(sound=5)


















area.r((raster(file.path(m2000,"mascara_pastoreo_praderas.tif"))+raster(file.path(m2000,"mascara_cultivos_matorral.tif")))*reclassify(dem,matrix(c(-Inf,1000,1000,Inf,0,1),ncol=3)))
area.r((raster(file.path(m2000,"mascara_pastoreo_praderas.tif"))+raster(file.path(m2000,"mascara_cultivos_matorral.tif")))*raster(file.path(z2018,"zona_vegalta.tif")))





