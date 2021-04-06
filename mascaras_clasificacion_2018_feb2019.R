##Código para hacer las mascaras
library(raster)
library(rgdal)
library(beepr)
rm(list=ls())

dir<-"E:/FONDECYT/benjamin/clasificacion2018/clasificaciones"
setwd(dir)

clas  <- raster("E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/noviembre2018/correccion_plantaciones_26112018")

ade   <- raster("E:/FONDECYT/benjamin/aderas.tif")
dem   <- crop(raster("E:/FONDECYT/benjamin/Landsat8/demangela_19s.tif"),clas)
h1200 <- reclassify(dem,matrix(c(-Inf,1200,1200,Inf,0,1),ncol=3))
h1000 <- reclassify(dem,matrix(c(-Inf,1000,1000,Inf,0,1),ncol=3))
h800  <- reclassify(dem,matrix(c(-Inf,800,800,Inf,0,1),  ncol=3))

'%ni%'<-Negate('%in%')

#======================================================================================================#
#===================================== Generación de las máscaras =====================================#
#======================================================================================================#

# mascara de artefactos por LaSRC  === USAR PARA ENMASCARAR LA CLASIFICACION ESPECTRAL PRIMERO ===
{
  #clas<-raster("febrero2019/MLC_2018_26022019_maj")
  
  #fuente<-extend(raster("noviembre2018/clas_correccion_artefactos_lasrc_2018"),clas)
  #mascara<-rasterize(readOGR("mascaras/mascara_para_lasrc.shp"),clas,field="Name",background=0,
  #                   filename="mascara_para_lasrc.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #mascara<-raster("mascara_para_lasrc.tif")
   #codigo para generar la máscara de SR Aerosol
  #if(file.exists("mascaras/sr_aerosol.tif")) {sr_aerosol<-raster("mascaras/sr_aerosol.tif")} else
  {
  #  sr_aerosol<-raster("mascaras/mascara_artefactos_lasrc.tif")
  #  ex_lt<-extent(sr_aerosol)
  #  exok<-extent(ex_lt[1],ex_lt[2],ex_lt[3]+10000000,ex_lt[4]+10000000)
  #  sr_aerosol<-setExtent(sr_aerosol,exok)
  #  proj4string(sr_aerosol)<-"+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #  sr_aerosol<-mask(crop(sr_aerosol,clas),mascara)
    
    # extraer pixeles para enmascarar
  #  sr_aerosol[sr_aerosol==80]<-1
  #  sr_aerosol[sr_aerosol!=1] <-0
  #  sr_aerosol[is.na(sr_aerosol)]<-0
  #  sr_aerosol<-mask(sr_aerosol,mascara,maskvalue=1,updatevalue=0)
  #  writeRaster(sr_aerosol,"mascaras/sr_aerosol.tif",format="GTiff",datatype="INT1U",overwrite=T)
  }
  
  #clas<-mask(clas,sr_aerosol,maskvalue=1,updatevalue=0)
  #mascara<-fuente*sr_aerosol
  
  #clas_enmascarada<-clas+mascara
  #plot(clas_enmascarada)
  
  #writeRaster(clas_enmascarada,"febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif",
  #            format="GTiff",datatype="INT1U",overwrite=T)
}

# bosques secundarios a achaparrados altos === LISTA === 
{
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  #adefocal<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_secundario_achaparrado_1984.shp"),clas,background=0,field="CLASS_ID")
  #adefocal<-rasterize(readOGR("mascaras/zonas_secundario_achaparrado2018.shp"),clas,background=0,field="CLASS_ID")
  #writeRaster(adefocal,"mascaras/zonas_secundario_achaparrado2018.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #adefocal<-raster("mascaras/zonas_secundario_achaparrado2018.tif")
  adefocal<-reclassify(raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/zonas_reclas_aa-ab-sec.tif"),
                       matrix(c(1,2,3,0,0,1),ncol=2))
  plot(adefocal)
  
  
  #zona_especial<-rasterize(readOGR("mascaras/zona_especial_achaparrado_primario.shp"),clas,field="CLASS_ID",background=0)
  #writeRaster(zona_especial,"mascaras/zona_especial_achaparrado_primario.tif",format="GTiff",datatype="INT1U",overwrite=T)
  zona_especial<-extend(raster("mascaras/zona_especial_achaparrado_primario.tif"),clas)
  h1100<-reclassify(dem,matrix(c(-Inf,1100,1100,Inf,0,1),ncol=3))
  #zonasur_h<-rasterize(readOGR("mascaras/zonasur_altura.shp"),clas,background=0,field="Name")
  #h600<-reclassify(dem,matrix(c(-Inf,600,600,Inf,0,1),ncol=3))*zonasur_h
  
  #mascara_altura<-mask(h1000,zonasur_h,maskvalue=1,updatevalue=0)+h600;plot(mascara_altura)
  
  #dejar solo las zonas de bosque secundario
  
  clas[clas%ni%c(41,42,43)]<-0
  clas[clas%in%c(41,42,43)]<-1
  plot(clas)
  #clas[clas==256]<-1
  
  #adefocal[is.na(adefocal)]<-0
  #adefocal[adefocal!=0]<-1
  
  
  masklv<-adefocal*clas#+clas*zona_especial*h1100;plot(masklv)
  #writeRaster(masklv,"noviembre2018/mascaras/mascara_secundarios_achaparradoLV.tif",format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(masklv,"febrero2019/mascaras/mascara_secundarios_achaparradoLV.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(masklv)
  rm(adefocal)
  rm(zona_especial)
  rm(mascara_agua_noclas)
  rm(zonas_aguas)
  rm(h1100)
  
}

# mascara glaciares  === LISTA === 
{
  glac<-raster("mascaras/mascara_glaciar.tif")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  glac[glac!=0]<-1
  glac<-extend(glac,clas)
  plot(glac)
  writeRaster(glac,"febrero2019/mascaras/mascara_glaciar.tif",format="GTiff",datatype="INT1U", overwrite=T)
  rm(glac)
  rm(clas)
}

## bosques achaparrados altos con h<1000 msnm a Secundarios === LISTA ===
{
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  
  clas[clas%in%c(1,2,3)]<-1
  clas[clas!=1]<-0
  
  mAB_achapabajo      <-reclassify(raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/zonas_reclas_aa-ab-sec.tif"),
                                   matrix(c(1,2,3,1,0,0),ncol=2))
  mAB_achapasecundario<-reclassify(raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/zonas_reclas_aa-ab-sec.tif"),
                                   matrix(c(1,2,3,0,1,0),ncol=2))
  
  mascara_achapaAlto_achapaBajo<- clas*mAB_achapabajo
  mascara_achapaAlto_secundario<- clas*mAB_achapasecundario
  mascara_achapaAlto_achapaBajo[is.na(mascara_achapaAlto_achapaBajo)]<-0
  mascara_achapaAlto_secundario[is.na(mascara_achapaAlto_secundario)]<-0
  
  plot(mascara_achapaAlto_achapaBajo)
  plot(mascara_achapaAlto_secundario)
  
  writeRaster(mascara_achapaAlto_secundario,"febrero2019/mascaras/mascara_achaparradoAlto_a_secundario.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(mascara_achapaAlto_achapaBajo,"febrero2019/mascaras/mascara_achaparradoAlto_a_achaparradoBajo.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(mascara_achapaAlto_achapaBajo)
  rm(mascara_achapaAlto_secundario)
  rm(mascaraAB)
  rm(h.bajo800.400)
  rm(mAB_achapasecundario)
  rm(mAB_achapabajo)
  rm(achapa_secundario)
  rm(achapaSur)
}

# máscara de las zonas no clasificadas cerca del agua === LISTA === 
{
  print("máscara de las zonas no clasificadas cerca del agua")
  a<-Sys.time();print(a)
  beep_on_error(
    {
      clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
      #zonas_agua<-rasterize(readOGR("mascaras/zonas_agua.shp"),clas,field="Id",background=0,
      #                      filename="mascaras/zonas_agua.tif",format="GTiff",datatype="INT1U", overwrite=T)
      zonas_agua<-raster("mascaras/zonas_agua.tif")
      agua<-clas*ade
      agua[agua%ni%c(8:11,39,40)]<-0
      agua[agua%in%c(8:11,39,40)]<-1
      
      #agua_p<-agua
      #pdte<-reclassify(terrain(dem,"slope",unit="degrees"),matrix(c(-Inf,15,15,Inf,0,1),ncol=3))
      #mascagua<-agua_p*pdte*reclassify(dem,matrix(c(-Inf,600,600,Inf,0,1),ncol=3))
      #plot(mascagua)
      #writeRaster(mascagua,"noviembre2018/mascaras/mascara_agua_pdte.tif",format="GTiff",datatype="INT1U",overwrite=T)#
      #
      #agua<-agua*reclassify(mascagua,matrix(c(1,0,0,1),ncol=2))#mask(agua,mascagua,maskvalue=1,updatevalue=0)
      plot(agua)
      agua[agua==0]<-NA
      agua<-mask(agua,raster("mascaras/zonas_laderas_aguas.tif"),maskvalue=1,updatevalue=NA)
      
      agua_shp<-rasterToPolygons(agua,na.rm=T,dissolve=T);shapefile(agua_shp,"mascaras/zonas_agua_paraNoclas.shp",overwrite=T)
      agua_shp<-shapefile("mascaras/zonas_agua_paraNoclas.shp")
      beep()
      bufagua_shp<-buffer(agua_shp,width=30);shapefile(bufagua_shp,"mascaras/buffer_agua_paraNoclas.shp",overwrite=T)
      bufagua<-shapefile("mascaras/buffer_agua_paraNoclas.shp")
      beep()
      plot(bufagua_shp)
      bufagua<-rasterize(bufagua_shp,clas,background=0,
                         filename="mascaras/buffer_noclas_a_agua.tif",format="GTiff",datatype="INT1U",overwrite=T)
      bufagua<-raster("mascaras/buffer_noclas_a_agua.tif")
      beep()
      plot(bufagua)
      
      noclas<-clas
      noclas[noclas%ni%c(0,44)]<-NA
      noclas[noclas%in%c(0,44)]<-1
      noclas[is.na(noclas)]<-0
      plot(noclas)
      
      mascara_noclas_agua<-noclas*bufagua+zonas_agua
      mascara_noclas_agua[mascara_noclas_agua!=0]<-1
      plot(mascara_noclas_agua)
      mascara_noclas_agua<-mask(mascara_noclas_agua,raster("mascaras/zonas_laderas_aguas.tif"),maskvalue=1,updatevalue=0)
      writeRaster(mascara_noclas_agua,"febrero2019/mascaras/mascara_bordes_agua.tif",format="GTiff",datatype="INT1U",overwrite=T)
      beep(sound=2)
      rm(clas)
      rm(agua)
      rm(agua_shp)
      rm(agua_p)
      rm(bufagua)
      rm(bufagua_shp)
      rm(zonas_agua)
      rm(mascara_noclas_agua)
      }
    ,sound=7)
  b<-Sys.time();print(b-a)
}

# terrenos sin vegetación sobre el límite arbóreo === LISTA ===
{
  ndvi<-stack("E:/FONDECYT/benjamin/clasificacion2018/bandas_sr_2018.tif")[[7]]
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")

  #pampa<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_praderas_densas_1984.shp"),clas,field="codigo",background=0)
  #writeRaster(pampa,"E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_praderas_densas_1984.tif",format="GTiff",datatype="INT1U",overwrite=T)
  pampa<-raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/limite_pampa_bosque.tif")
  
  #máscara para zonas no clasificadas
  noclas<-clas
  noclas[noclas==0]<-256
  noclas[noclas!=256]<-0;noclas[noclas==256]<-1
  mascara_noclas<-noclas*h1000
  #mascar para zonas de agua
  #clas[clas%ni%c(10:13)]<-0
  #clas[clas%in%c(10:13)]<-1
  
  ndvi<-reclassify(ndvi,matrix(c(-Inf,0.5,0.5,Inf,1,0),ncol=3))*ade #mascara de zonas con NDVI < 0.25
  
  # máscara con todas las zonas clasificadas como sin vegetación sobre el límite de la vegetación
  clas[clas%ni%c(12:18,21:23,45:49)]<-0
  clas[clas%in%c(12:18,21:23,45:49)]<-1
  sv_altura<-clas*h1000
  
  mascara_sv_altura<-mask(ndvi*mascara_noclas+sv_altura,pampa,maskvalue=1,updatevalue=0)
  mascara_sv_altura<-mask(mascara_sv_altura,raster("febrero2019/mascaras/mascara_glaciar.tif"),maskvalue=1,updatevalue=0)
  mascara_vegetacion_altura_noclas<-mask(mask(mascara_noclas,ndvi,maskvalue=1,updatevalue=0),pampa,maskvalue=1,updatevalue=0)
  
  plot(mascara_sv_altura)
  plot(mascara_vegetacion_altura_noclas)
  sv_altura[sv_altura==2]<-1
  
  plot(mascara_sv_altura)
  
  writeRaster(mascara_sv_altura,"febrero2019/mascaras/mascara_sv_altura.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(noclas)
  rm(mascara_sv_altura)
  rm(mascara_noclas)
  rm(mascara_vegetacion_altura_noclas)
  rm(sv_altura)
  rm(zonas_nubes)
  rm(ndvi)
  
  
  
}

# mascara urbano === LISTA === 
{
  print("mascara urbano")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")*ade
  m.urbano<-mask(extend(raster("febrero2019/ISODATA_2018_urbanoCoy_27022019"),clas),readOGR("febrero2019/shapes/zona_urbano.shp"))
  m.urbano[m.urbano%ni%1:3]<-0
  m.urbano[m.urbano%in%1:3]<-1
  #zonas_urbano<-rasterize(readOGR("mascaras/zonas_urbano.shp"),clas,field="Name",background=0,
  #                        filename="mascaras/zonas_urbano.tif",format="GTiff",datatype="INT1U",overwrite=T)
  zonas_urbano<-raster("mascaras/zonas_urbano.tif")
  
  m.urbano[is.na(m.urbano)]<-0
  urbano<-m.urbano+zonas_urbano
  urbano[urbano%in%1:2]<-1
  
  plot(urbano)
  writeRaster(urbano,"febrero2019/mascaras/mascara_urbano.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(m.urbano)
  rm(clas)
  rm(zonas_urbano)
  rm(urbano)
}

# máscara de zonas de suelo desnudo a pampa === LISTA === 
{
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  #pampa<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/mascaras/zonas_sv_pampaabierta_septiembre.shp"),clas,field="Name",background=0,
  #                 filename="mascaras/zonas_sv_pampa.tif",format="GTiff",datatype="INT1U",overwrite=T)
  pampa<-raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/limite_pampa_bosque.tif")
  #plot(pampa)
  #zonas_agua<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_lechos_pampa_1984.shp"),clas,field="CLASS_ID",background=0,
  #                      filename="mascaras/zonas_agua_pampa.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #zonas_agua<-raster("mascaras/zonas_agua_pampa.tif")
  
  sv<-clas*ade
  sv[sv%ni%c(0,12:18,45:49)]<-256
  sv[sv%in%c(0,12:18,45:49)]<-1
  sv[sv==256]<-0
  estepa<-clas*ade
  estepa[estepa%ni%c(21:23,34,35)]<-0
  estepa[estepa%in%c(21:23,34,35)]<-1
  
  sv_pampa<-mask(sv*pampa,raster("febrero2019/mascaras/mascara_urbano.tif"),maskvalue=1,updatevalue=0)
  estepa_ok<-estepa*pampa
  
  plot(sv_pampa)
  plot(estepa_ok)
  sv_pampa<-mask(sv_pampa,raster("febrero2019/mascaras/mascara_bordes_agua.tif"),maskvalue=1,updatevalue=0)
  
  writeRaster(sv_pampa,"febrero2019/mascaras/mascara_sv_a_pampaAbierta.tif",format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(estepa_ok,"febrero2019/mascaras/mascara_estepa.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(sv)
  rm(sv_pampa)
  rm(pampa)
  rm(estepa)
  rm(estepa_ok)
}

# Máscara de Vegetación en altura === LISTA === 
{
  print("Máscara de Vegetación en altura")
  ##  trabajar con clases praderas y con matorrales abiertos y sin estepa  (id= [27,28,35,36])
  #zona_h1000     <- raster("febrero2019/zonas/zona_h1000.tif")
  #zona_novegalta <- raster("febrero2019/zonas/zona_exclusion_sv_altura.tif")
  #zona_vegalta<-mask(h1200+h1000*zona_h1000,zona_novegalta,maskvalue=1,updatevalue=0)
  #zona_vegalta[zona_vegalta==2]<-1
  #writeRaster(zona_vegalta,"febrero2019/zonas/zona_vegalta.tif",
  #            format="GTiff",datatype="INT1U",overwrite=T)
  
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  clas[clas%in%c(25,26,27,28,29,34,35)]<-256
  clas[clas!=256]<-0
  clas[clas==256]<-1
  plot(clas)
  
  #mascara_pradera_altura<-clas*h1200
  mascara_pradera_altura<-mask(clas*raster("febrero2019/zonas/zona_vegalta.tif"),
                               raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/limite_pampa_bosque.tif"),maskvalue=1,updatevalue=0)
  
  plot(mascara_pradera_altura)
  writeRaster(mascara_pradera_altura,"febrero2019/mascaras/mascara_vegetacion_en_altura.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(mascara_pradera_altura)
}  

# reclasificacion de praderas y agricola a otras cosas === LISTA ===
{
  print("reclasificacion zonas agricolas y praderas")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  zonas_cultivos_noreclas<-raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/zona_agro_noreclas.tif")
  # ======== trabajar con la clase Cultivos solamente (id=5 y 6) ======== #
  {
    cultivo<-clas*ade
    cultivo[cultivo%ni%6:7]<-0
    cultivo[cultivo%in%6:7]<-1
    
    # eliminar de las mascaras de cultivos todas las zonas que están dentro de lo que se considera zona de cultivos
    #zonas_cultivos_noreclas<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/mascaras/zona_agro_noreclas.shp"),clas,field="Id",background=0,
    #                                   filename="E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/mascaras/zona_agro_noreclas.tif",format="GTiff",datatype="INT1U",overwrite=T)
    
    #zona pa reclasificar cultivos a achaparrado
    #zonas_cultivos_achaparrado<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_agricola_achaparrado_1984.shp"),clas,field="CLASS_ID", background=0,
    #                                      filename="mascaras/zonas_agricola_achaparrado_2018.tif",format="GTiff",datatype="INT1U",overwrite=T)
    zonas_cultivos_achaparrado<-raster("mascaras/zonas_agricola_achaparrado_2018.tif")
    zonas_cultivos_achaparrado[zonas_cultivos_achaparrado==2]<-1
    
    cultivo<-mask(cultivo,zonas_cultivos_noreclas,maskvalue=1,updatevalue=0)
    # hacer máscara en función de la pendiente. 
    # Cultivos con pdte > 7.5° pasan a matorrales
    # cultivos con pdte < 7.5° pasan a praderas
    
    slo12.5<-reclassify(terrain(dem,"slope",unit="degrees")*ade,matrix(c(-Inf,12.5,12.5,Inf,0,1),ncol=3))
    plot(slo12.5)
    
    mascara_cultivos_matorral <- mask(cultivo*slo12.5,zonas_cultivos_achaparrado,maskvalue=1,updatevalue=0)
    mascara_cultivos_praderas <- mask(mask(cultivo,slo12.5,maskvalue=1,updatevalue=0),zonas_cultivos_achaparrado,maskvalue=1,updatevalue=0)
    mascara_cultivos_achaparrado<-cultivo*zonas_cultivos_achaparrado
    
    plot(mascara_cultivos_matorral)
    plot(mascara_cultivos_praderas)
    plot(mascara_cultivos_achaparrado)
    
    writeRaster(mascara_cultivos_achaparrado,"febrero2019/mascaras/mascara_cultivos_achaparrado.tif",format="GTiff",datatype="INT1U",overwrite=T)
    writeRaster(mascara_cultivos_matorral,"febrero2019/mascaras/mascara_cultivos_matorral.tif",format="GTiff",datatype="INT1U",overwrite=T)
    writeRaster(mascara_cultivos_praderas,"febrero2019/mascaras/mascara_cultivos_praderas.tif",format="GTiff",datatype="INT1U",overwrite=T)
  }
  # ======= fin de la reclasificacion de las zonas de cultivos ======= #
  
  # ======= reclasificacion de barbecho fuera del área de cultivo a otras cosas ======= #
  {
    zona_pampa <- raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/limite_pampa_bosque.tif")
    
    barbecho<-clas*ade
    
    barbecho[barbecho%ni%45:49]<-0
    barbecho[barbecho%in%45:49]<-1
    
    barbecho.pradera <- mask(barbecho,zonas_cultivos_noreclas, maskvalue=1,updatevalue=0) # sacar zonas de cultivos
    barbecho.pradera <- mask(barbecho.pradera,zona_pampa,      maskvalue=1,updatevalue=0) # sacar zonas de pampa
    barbecho.pradera <- mask(barbecho.pradera,h1000,           maskvalue=1,updatevalue=0) # sacar zonas que serán sd de altura
    barbecho.pradera <- mask(barbecho.pradera,raster("febrero2019/mascaras/mascara_urbano.tif"),maskvalue=1,updatevalue=0)
    
    plot(barbecho.pradera)
    
    writeRaster(barbecho.pradera,"febrero2019/mascaras/mascara_barbecho_pradera.tif",
                format="GTiff",datatype="INT1U",overwrite=T)
    rm(barbecho)
    rm(barbecho.pradera)
  }
  
  # ======= trabajar con la clase de Pastoreo  ======= #
  {
    zonas_agricolas<-raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/zonas_pastoreo_noreclas_2000.tif")
    plot(zonas_agricolas)
    agro<-clas
    agro[agro%ni%c(4:5)]<-0
    agro[agro%in%c(4:5)]<-1
    agro<-mask(agro,zonas_agricolas,maskvalue=1,updatevalue=0)
    
    mascara_agro_a_praderas<-mask(agro,zonas_cultivos_achaparrado,maskvalue=1,updatevalue=0)
    mascara_agro_achaparrado<-agro*zonas_cultivos_achaparrado
    
    plot(mascara_agro_a_praderas)
    plot(mascara_agro_achaparrado)
    
    pastoreo<-raster("E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/MLC_2018_23102018_grande")
    
    pastoreo[pastoreo%ni%c(6,7,8)]<-0
    pastoreo[pastoreo%in%c(6,7,8)]<-1
    
    pdte<-reclassify(terrain(dem,"slope",unit="degrees"),matrix(c(-Inf,5.71,5.71,Inf,1,0),ncol=3))
    
    mascara_pastoreo<-pastoreo*pdte*zonas_agricolas
    plot(mascara_pastoreo)
    
    writeRaster(mascara_agro_a_praderas,"febrero2019/mascaras/mascara_pastoreo_praderas.tif",format="GTiff",datatype="INT1U",overwrite=T)
    writeRaster(mascara_pastoreo,"febrero2019/mascaras/mascara_pastoreo_plano.tif",format="GTiff",datatype="INT1U",overwrite=T)
    writeRaster(mascara_agro_achaparrado,"febrero2019/mascaras/mascara_pastoreo_achaparrado.tif",format="GTiff",datatype="INT1U",overwrite=T)
    
    rm(clas)
    rm(slo12.5)
    rm(slo7.5)
    rm(mascara_cultivos_matorral)
    rm(mascara_cultivos_praderas)
    rm(mascara_agro_a_praderas)
    rm(mascara_pastoreo)
    rm(agro)
    rm(cultivo)
    rm(zonas_agricolas)
    rm(zonas_cultivos)
    rm(zonas_cultivos_noreclas)
    rm(pdte)
    rm(pastoreo)
    rm(zonas_cultivos_achaparrado)
    rm(mascara_cultivos_achaparrado)
    rm(mascara_agro_achaparrado)
  }
}

# reclas de zonas agricolas de pastoreo en ñirehuao a cultivos === LISTA ===
{
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  print("reclasificacion de zonas de pastoreo en ñirehuao")
  adecultivos<-rasterize(readOGR("mascaras/zonas_pastoreo_cultivos.shp"),clas,field="Id",background=0)
  clas[clas%in%c(4,5)]<-256
  clas[clas!=256]<-0
  clas[clas==256]<-1
  
  mascara_cultivos<-adecultivos*clas
  plot(mascara_cultivos)
  writeRaster(mascara_cultivos,"febrero2019/mascaras/mascara_pastoreo_cultivos.tif",format="GTiff",datatype="INT1U",overwrite=T)
}

# mascara de zonas de praderas pampas densas === LISTA ===
{
  clas<-raster("MLC_2018_11092018")
  zonas_praderas_pampas_densas<-rasterize(readOGR("E:/FONDECYT/benjamin/clasificacion1984/mascaras/zonas_pampasdensas_1984.shp"),clas,field="codigo",background=0)
  plot(zonas_praderas_pampas_densas)
  
  pradens<-clas*ade
  pradens[pradens!=9]<-0
  pradens[pradens==9]<-1
  plot(pradens)
  
  pradens<-pradens*zonas_praderas_pampas_densas
  plot(pradens)
  
  writeRaster(pradens,"febrero2019/mascaras/mascara_praderas_densas.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(pradens)
}

# máscara de zonas quemadas === LISTA ===
{clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")

fuego<-rasterize(readOGR("mascaras/zonas_fuego.shp"),clas,field="Name",background=0)
writeRaster(fuego,"febrero2019/mascaras/mascara_fuego.tif",format="GTiff",datatype="INT1U",overwrite=T)
rm(clas)
}

# completar máscara de cuncunas === LISTA ===
{
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  #cun<-clas
  #cun[cun%ni%c(19,20)]<-0
  #cun[cun%in%c(19,20)]<-1
  #cunfoc<-focal(cun,matrix(rep(1,25),ncol=5),fun=function(x){return(modal(x,na.rm=T))}) #hacer filtro modal con matrix de 3x3 para eliminar pixeles solos#
  #cunfoc[cunfoc==0]<-NA
  #writeRaster(cunfoc,"febrero2019/mascaras/cuncunas_focal_2018.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #cunfoc<-raster("febrero2019/mascaras/cuncunas_focal_2018.tif")
  #cunpol<-rasterToPolygons(cunfoc, na.rm=T,dissolve=T);shapefile(cunpol,"febrero2019/shapes/cuncunas_focal_2018.shp",overwrite=T)
  #buffer_cuncunas<-buffer(cunpol,width=30) #hacer un buffer de 1 pixel al rededor de las zonas clasificadas como cuncunas
  #shapefile(buffer_cuncunas,"mascaras/shapes/buffer_rocas_cuncunas2018.shp",overwrite=T)
  buffer_cuncunas<-shapefile("mascaras/shapes/buffer_rocas_cuncunas2018.shp")
  beepr::beep()
  
  # extraer solo las rocas 
  clas[clas%ni%c(12:18)]<-0
  clas[clas%in%c(12:18)]<-1
  rocas_cuncunas<-mask(clas,buffer_cuncunas,updatevalue=0)
  plot(rocas_cuncunas)
  
  rocas_cuncunas<-mask(rocas_cuncunas,raster("febrero2019/mascaras/mascara_sv_altura.tif"),maskvalue=1,updatevalue=0)
  rocas_cuncunas<-mask(rocas_cuncunas,raster("mascaras/zonas_cuncunas_2018.tif"))
  writeRaster(rocas_cuncunas,"febrero2019/mascaras/mascara_rocas_a_cuncunas.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(cun)
  rm(cunfoc)
  rm(cunpol)
  rm(buffer_cuncunas)
  rm(rocas_cuncunas)
  beep(sound=2)
}

# mascara de cuncunas a otras cosas === LISTA ===
{
  print("mascara de cuncunas a otras cosas")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  b_deg<-raster("noviembre2018/bosques_degradados_2018")
  
  zonas_cuncunas          <-rasterize(readOGR("mascaras/zonas_cuncunas_2018.shp"),clas, field="Id",background=0,
                                 filename="mascaras/zonas_cuncunas_2018.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #zonas_cuncunas_total     <- rasterize(readOGR("noviembre2018/zonas_cuncunas_correccion_total.shp"),clas,field="Id",background=0,
  #                                filename="mascaras/zonas_cuncunas_2018_correccion_total.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #zonas_cuncunas_espectral <- rasterize(readOGR("noviembre2018/zonas_cuncunas_correccion_espectral.shp"),clas,field="Id",background=0,
  #                                filename="mascaras/zonas_cuncunas_2018_correccion_espectral.tif",format="GTiff",datatype="INT1U",overwrite=T)
  
  zonas_cuncunas           <- raster("mascaras/zonas_cuncunas_2018.tif")
  zonas_cuncunas_total     <- raster("mascaras/zonas_cuncunas_2018_correccion_total.tif")
  zonas_cuncunas_espectral <- raster("mascaras/zonas_cuncunas_2018_correccion_espectral.tif")
  
  zonas_cuncunas[zonas_cuncunas==2]<-1
  zonas_cuncunas_total[zonas_cuncunas_total==2]<-1 # con esto puedo jugar a ver si utilizo todas las zonas definidas como brotes de cuncunas o las zonas optimizadas en noviembre. 
  zonas_cuncunas_espectral[zonas_cuncunas_espectral==2]<-1
  
  # zonas
  zonas_primario<-raster("mascaras/zonas_primarios.tif")
  zonas_matorral_denso<-raster("mascaras/zonas_matorral_denso.tif")
  zonas_matorral_denso<-mask(zonas_matorral_denso,zonas_primario,maskvalue=1,updatevalue=0)
  
  # correccion espectral
  {
    c_es<-b_deg*zonas_cuncunas_espectral
    c_es[c_es%ni%c(19,20,48)]<-0
    c_es[c_es%in%c(19,20,48)]<-1
  }
  
  # correccion total
  {
    c_t<-b_deg*zonas_cuncunas_total
    c_t[c_t%ni%c(19,20,28,29,48)]<-0
    c_t[c_t%in%c(19,20,28,29,48)]<-1
  }
  
  cun_mask<-c_es+c_t
  cun_mask[cun_mask==2]<-1
  
  cun<-clas
  
  cun[cun%ni%c(19,20)]<-0
  cun[cun%in%c(19,20)]<-1
  
  cun<-cun+cun_mask#mask(cun,zonas_cuncunas,maskvalue=1,updatevalue=0)+cun_mask
  cun[cun==2]<-1
  
  # máscara de bosques degradados dentro de las zonas de cuncunas, se reclasificarán estas zonas a bosque primario
  bosques_deg<-mask(b_deg*zonas_cuncunas,zonas_matorral_denso,maskvalue=1,updatevalue=0)
  bosques_deg[bosques_deg%ni%c(47)]<-0
  bosques_deg[bosques_deg%in%c(47)]<-1
  writeRaster(bosques_deg,"febrero2019/mascaras/mascara_bosquesdeg_primario.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  
  cun<-cun+raster("febrero2019/mascaras/mascara_rocas_a_cuncunas.tif")#+matorrales*zonas_cuncunas
  
  mascara_cuncunas_achaparrado    <- cun*h1200
  mascara_cuncunas_primario       <- cun*zonas_primario-mascara_cuncunas_achaparrado+bosques_deg
  mascara_cuncunas_matorral_denso <- cun*zonas_matorral_denso
  mascara_cuncunas_secundario     <- cun-mascara_cuncunas_primario-mascara_cuncunas_matorral_denso-mascara_cuncunas_achaparrado-bosques_deg
  
  mascara_cuncunas_achaparrado[mascara_cuncunas_achaparrado==2]<-1
  mascara_cuncunas_secundario[mascara_cuncunas_secundario==2]<-1
  mascara_cuncunas_primario[mascara_cuncunas_primario==2]<-1
  mascara_cuncunas_matorral_denso[mascara_cuncunas_matorral_denso==2]<-1
  
  writeRaster(mascara_cuncunas_primario,
              "febrero2019/mascaras/mascara_cuncunas_primario.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(mascara_cuncunas_secundario,
              "febrero2019/mascaras/mascara_cuncunas_secundario.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(mascara_cuncunas_achaparrado,
              "febrero2019/mascaras/mascara_cuncunas_achaparrado.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(mascara_cuncunas_matorral_denso,
              "febrero2019/mascaras/mascara_cuncunas_matorral_denso.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  
  rm(mascara_cuncunas_primario)
  rm(mascara_cuncunas_secundario)
  rm(mascara_cuncunas_achaparrado)
  rm(mascara_cuncunas_matorral_denso)
  rm(clas)
  rm(cun)
  rm(zonas_cuncunas)
  rm(b_deg)
  rm(bosques_deg)
  rm(cun_mask)
  rm(c_es)
  rm(c_t)
}

# reclasificacion de aguas en pendientes y wea === LISTA ===
{
  print("reclasificacion aguas en pendientes")
  agua<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")*ade
  zonas_aguas<-raster("mascaras/zonas_laderas_aguas.tif")
  
  agua[agua%ni%c(8:11,24,39,40)]<-0
  agua[agua%in%c(8:14,24,39,40)]<-1
  
  mascara_agua_noclas<-agua*zonas_aguas
  plot(mascara_agua_noclas)
  
  writeRaster(mascara_agua_noclas,"febrero2019/mascaras/mascara_agua_noclas.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  rm(zonas_agua)
  rm(zonas_aguas)
  rm(mascara_agua_noclas)
  rm(agua)
  rm(clas)
}

# nueva mascara de humedales === LISTA ===
{
  print("nueva mascara de humedales")
  clas1<-raster("MLC_14082018_3")
  clas2<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  clas3<-raster("MLC_2018_12092018")
  
  zonas_agua<-raster("mascaras/zonas_agua.tif")
  clas1<-extend(clas1,clas2)
  zonas_sombras<-raster("mascaras/zonas_laderas_aguas.tif")
  zonas_sombras[is.na(zonas_sombras)]<-0
  
  aguabase<-clas1
  aguabase<-mask(aguabase,zonas_sombras,maskvalue=1,updatevalue=0)
  humedales<-clas2
  sombras<-clas2
  
  aguabase[aguabase%ni%10:13]<-0
  aguabase[aguabase%in%10:13]<-1
  aguabase<-aguabase+zonas_agua
  
  humedales[humedales%ni%24]<-0
  humedales[humedales%in%24]<-1
  
  mascara_humedales<-mask(humedales,aguabase,maskvalue=1,updatevalue=0)
  mascara_humedales<-mask(mascara_humedales,zonas_sombras,maskvalue=1,updatevalue=0)
  plot(mascara_humedales)
  
  sombras[sombras!=44]<-0
  sombras[sombras==44]<-1
  mascara_sombras<-humedales*zonas_sombras+sombras
  plot(mascara_sombras)
  
  mascara_sombras<-mask(mascara_sombras,
                        reclassify(raster("febrero2019/mascaras/mascara_bordes_agua.tif")+raster("mascaras/congresoCIEP/mascara_aguabase.tif"),matrix(c(2,1),ncol=2)),
                        maskvalue=1,updatevalue=0)
  
  writeRaster(mascara_humedales,"febrero2019/mascaras/mascara_humedales.tif",format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(mascara_sombras,"febrero2019/mascaras/mascara_sombras.tif",format="GTiff",datatype="INT1U",overwrite=T)
  writeRaster(aguabase,"febrero2019/mascaras/mascara_aguabase.tif",format="GTiff",datatype="INT1U",overwrite=T)
  rm(mascara_humedales)
  rm(aguabase)
  rm(mascara_sombras)
  rm(sombras)
  rm(humedales)
  rm(clas1)
  rm(clas2)
  rm(clas3)
  rm(zonas_sombras)
  rm(zonas_agua)
}

# Máscara de plantaciones y Bosques primarios === LISTA ===
{
  print("mascara de plantaciones")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  #plantaciones<-rasterize(readOGR("noviembre2018/zonas_plantaciones_base.shp"),clas,field="ID",background=0,
  #                        filename="mascaras/zonas_plantaciones_2018.tif",format="GTiff",datatype="INT1U",overwrite=T)
  plantaciones<-raster("mascaras/zonas_plantaciones_2018.tif")
  
  prim<-clas
  prim[prim%ni%c(36:38,41:43)]<-0
  prim[prim%in%c(36:38,41:43)]<-1
  
  mascara_prim_plant<-prim*plantaciones
  plot(mascara_prim_plant)
  
  
  writeRaster(mascara_prim_plant,"febrero2019/mascaras/mascara_primario_plantaciones.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(prim)
  rm(mascara_prim_plant)
  rm(plantaciones)
  
}

# Máscara de bosques alterados en cercanías a la pampa === LISTA ===
{
  print("mascara de bosques alterados cerca de la pampa")
  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  b_deg<-raster("E:/FONDECYT/benjamin/clasificacion2018/clasificaciones/noviembre2018/bosques_degradados_2018_22112018")
  #h1100<-reclassify(dem,matrix(c(-Inf,1100,1100,Inf,0,1),ncol=3))
  #buf_pampa<-rasterize(buffer(readOGR("E:/FONDECYT/benjamin/clasificacion2002/mascaras/shapes/limite_pampa_bosque.shp"),5000),clas,background=0)-raster("E:/FONDECYT/benjamin/clasificacion2002/mascaras/zonas/limite_pampa_bosque.tif")
  #buf_pampa<-buf_pampa*ade
  #plot(buf_pampa)
  #writeRaster(buf_pampa,"febrero2019/zonas/zonas_bufer_pampa_bosques_degradados.tif",
  #           format="GTiff",datatype="INT1U",overwrite=T)
  buf_pampa<-raster("febrero2019/zonas/zonas_bufer_pampa_bosques_degradados.tif")*ade
  
  md<-clas
  md[md%ni%c(29,30)]<-0
  md[md%in%c(29,30)]<-1
  
  b_deg[b_deg!=48]<-0
  b_deg[b_deg==48]<-1
  
  mascara_bosques_degradados<-md*b_deg*buf_pampa
  plot(mascara_bosques_degradados)
  
  writeRaster(mascara_bosques_degradados,"febrero2019/mascaras/mascara_bosques_degradados.tif",
              format="GTiff",datatype="INT1U",overwrite=T)
  rm(clas)
  rm(buf_pampa)
  rm(md)
  rm(b_deg)
  rm(mascara_bosques_degradados)
  
}

stop()

#======================================================================================================#
#======================================== Máscaras POR REVISAR ========================================#
#======================================================================================================#

# mascara de achaparrados altos a primarios en zonas definidas con h entre 500 y 1000 msnm === OBSOLETA PORQUE SE FUSIONARON LAS CLASES ===
{
#  clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
#  zona_primario_norte<-rasterize(readOGR("mascaras/zonas_achaparrado_primario.shp"),
#                                 clas,field="CLASS_ID",background=0)
#  zona_primario_sur  <-rasterize(readOGR("mascaras/zonas_achaparrado_primario_sur.shp"),
#                                 clas,field="CLASS_ID",background=0)
#  zona_especial      <-rasterize(readOGR("mascaras/zona_especial_achaparrado_primario.shp"),
#                                 clas,field="CLASS_ID",background=0)
#  
#  h600.1000<-reclassify(dem,matrix(c(-Inf,600,1000,600,1000,Inf,0,1,0),ncol=3))
#  h400.800<-reclassify(dem,matrix(c(-Inf,400,800,400,800,Inf,0,1,0),ncol=3))
#  h1000.1100<-reclassify(dem,matrix(c(-Inf,1000,1100,1000,1100,Inf,0,1,0),ncol=3))
#  
#  clas[clas%ni%1:3]<-0
#  clas[clas%in%1:3]<-1
#  
#  mascara_achaparrado_primario<-h600.1000*zona_primario_norte*clas+h400.800*zona_primario_sur*clas+h1000.1100*zona_especial*clas
#  plot(mascara_achaparrado_primario)
#  mascara_achaparrado_primario<-mask(mascara_achaparrado_primario,raster("mascaras/congresoCIEP/mascara_achaparradoAlto_a_secundario.tif",maskvalue=1,updatevalue=0))
#  
#  writeRaster(mascara_achaparrado_primario,"mascaras/congresoCIEP/mascara_achaparrado_primario.shp",format="GTiff",datatype="INT1U",overwrite=T)
#  rm(clas)
#  rm(mascara_achaparrado_primario)
#  rm(h400.1000)
#  rm(h400.800)
#  rm(h1000.1100)
}

# bosques achaparrados bajos hacia la pampa    === NO USAR PORQUE NO ESTÁ ESA CLASE ===
# la idea es generar una máscara que sirva para reclasificar las zonas de bosques achaparrados bajos fuera de la pampa 
# y a menos de 800msnm a secundario, y una que sirva para reclasificar las zonas sobre 800msnm a bosque achaparrado
{
  #clas<-raster("febrero2019/clasificacion2018_enmascarada_artefactos_27022019.tif")
  #h.bajo400<-raster("mascaras/h_bajo400.tif")
  #mascara_AB<-raster("mascaras/zonas_pampa_oriental_ABfinal.tif")
  #zonas_oriental_pampa<-rasterize(readOGR("mascaras/zona_oriental_achaparradosBajos.shp"),clas,field="Name",background=0)*ade  #generacion de la máscara de zonas orientales pampa nueva
  
  #alturas<-reclassify(dem,matrix(c(-Inf,800,800,Inf,1,0),ncol=3))*ade
  #mascara_AB<-(zonas_oriental_pampa*alturas)+alturas+1
  #mascara_AB<-reclassify(mascara_AB,matrix(c(1,2,3,3,2,1),ncol=2))
  #plot(mascara_AB)
  #writeRaster(mascara_AB,"mascaras/zonas_pampa_oriental_ABfinal.tif",format="GTiff",datatype="INT1U",overwrite=T)
  
  #mascaraAB<-raster("mascaras/zonasAB_cortada.tif")*ade
  #clas[clas!=4]<-0
  #clas[clas==4]<-1
  # zonas que se reclasificaran de bosque achaparrado bajo a secundario
  #maskAB.secundario<-mask(clas,reclassify(mascara_AB,matrix(c(1,2,3,1,0,1),ncol=2)),
  #                        maskvalue=1,updatevalue=0)*h.bajo400 #seleccionar solo achaparrados bajos fuera de las pampas y bajo 1000 msnm, y bajo los 400msnm en los interiores montañosos del norte y sur#
  #plot(maskAB.secundario)
  # zonas que se reclasificarán de bosques achaparrados bajos a bosque achaparrado alto
  #maskAB.achaparrado<-mask(clas,reclassify(mascara_AB,matrix(c(1,2,3,1,1,0),ncol=2)),
  #                         maskvalue=1,updatevalue=0)*ade
  #maskAB.achaparrado[maskAB.achaparrado!=1]<-0
  #plot(maskAB.achaparrado)
  
  #writeRaster(maskAB.secundario,"noviembre2018/mascaras/mascara_achaparradobajo_a_secundario.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #writeRaster(maskAB.achaparrado,"noviembre2018/mascaras/mascara_achaparradobajo_a_achaparraado.tif",format="GTiff",datatype="INT1U",overwrite=T)
  #rm(clas)
}
