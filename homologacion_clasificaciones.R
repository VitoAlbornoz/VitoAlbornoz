# homologación de clases entre aysen y biobio


# Aysen
clasificaciones<-stack(
  c1984<-raster("productos_finales/clasificaciones/aysen_1984_sintetizada.tif"),
  c2000<-raster("productos_finales/clasificaciones/aysen_2000_sintetizada.tif"),
  c2018<-raster("productos_finales/clasificaciones/aysen_2018_sintetizada.tif"))

clases<-{c("Bosque Adulto",       #1
           "Renoval",             #2
           "Matorrales",          #3
           "Praderas",            #4
           "Agropecuario",        #5
           "Agua",                #6
           "Humedales",           #7 
           "Suelo Desnudo",       #8
           "Plantación Forestal", #9
           "Incendio",            #10
           "Urbano",              #11
           "Estepa",              #12
           "Estepa abierta",      #13
           "Glaciares-Nieve",     #14
           "Vegetación de altura",#15
           "No clasificado")}     #16
           


clases<-{c("agua", #2
           "agropecuario", #3
           "estepa", #4
           "fuego", #5
           "glaciares_nieve", #6
           "humedales", #7
           "matorrales", #8 
           "plantacion forestal", #10
           "praderas", #11
           "Bosque adulto", #13
           "renoval", #14
           "suelo desnudo", #15
           "vegetacion en altura", #17
           "urbano", #18
           "estepa abierta", #19
           "no clasificado")} #20

clasificaciones <- reclassify(clasificaciones,matrix(c(2,3,4 ,5 ,6 ,7,8,10,11,13,14,15,17,18,19,20,
                     6,5,12,10,14,7,3,9 ,4 ,1 ,2 ,8 ,15,11,13,16),ncol=2))
writeRaster(clasificaciones,"productos_finales/clasificaciones/aysen_sintetizada_homologada.tif",format="GTiff",datatype="INT1U",bylayer=T,overwrite=T)


# biobio
ruta<-"E:/FONDECYT/biobio"
setwd(ruta)
clasificaciones<-stack(
  c1986<-raster("analisis_paisaje/Reclass_1986.tif"),
  c2001<-raster("analisis_paisaje/Reclass_2001.tif"),
  c2011<-raster("analisis_paisaje/Reclass_2011.tif"),
  c2017<-raster("analisis_paisaje/Reclass_2017.tif"))

dem.bb <- raster("dem_biobio.tif")

clases<-{c("Bosque Adulto",        #1
           "Renoval",              #2
           "Matorrales",           #3
           "Praderas",             #4
           "Agropecuario",         #5
           "Agua",                 #6
           "Humedales",            #7 
           "Suelo Desnudo",        #8
           "Plantación Forestal",  #9
           "Incendio",             #10
           "Urbano",               #11
           "Estepa",               #12
           "Estepa abierta",       #13
           "Glaciares-Nieve",      #14
           "Vegetación de altura", #15
           "No clasificado")}      #16 
 
clases<-{c("bosque adulto",        #1
           "renoval",              #2
           "plantaciones",         #3
           "matorrales",           #4
           "agricola",             #5
           "suelo desnudo",        #6
           "agua",                 #7 
           "humedales",            #8
           "nieve",                #9
           "urbano",               #10
           "praderas",             #11
           "no clasificado",       #12
           "incendio")}            #13

clasificaciones <- reclassify(clasificaciones,matrix(c(1,2,3,4,5,6,7,8,9 ,10,11,12,13,
                                                       1,2,9,3,5,8,6,7,14,11,4, 16,10),ncol=2))

for(i in 1:4) {clasificaciones[[i]][clasificaciones[[i]] %in% c(3,4,7) & dem.bb >= 1900] <- 15}
unique(clasificaciones[[1]])
writeRaster(clasificaciones,"clasificaciones/biobio_homologada.tif",format="GTiff",datatype="INT1U",bylayer=T,overwrite=T)
