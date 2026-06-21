## Cargo todas las librerias a usar
rm(list=ls())
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ppcor)
library(pheatmap)
library(igraph)
library(FOCI)
library(VineCopula)
library(RColorBrewer)
library(XICOR)
library(dHSIC)
library(RColorBrewer)
library(CondCopulas)
library(visNetwork)
library(geomnet)
library(igraph)
library(corrplot)
library(ggnet)
library(dplyr)
library(corrplot)
library(readxl)

## Cargo los datos que voy a usar

AnotacionesE <- read_excel("C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/scripts/DatosSinTrasncripcion.xlsx")
transcD <- read.csv2("C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/scripts/Developmental_transcriptome_Drosophila.csv")
names(transcD)<- c("Identificador", names(transcD)[-1])

## Creo y cargo las funciones que uso despues

# Normalizacion 

est1<-function(u){u/max(u)}
est2<-function(u){ qnorm((rank(u) - 0.5)/length(u))}

# Funciones de correlacion 

corp<-function(u,v,metodo=mmi){ 
  u<-as.numeric(u)
  v<-as.numeric(v)  
  fg<-c()
  for(jj in 1: (24-larg)){ 
    fg[jj] <- cor(u[jj:(jj+larg)],v[jj:(jj+larg)],
                  method = metodo)
  }
  max(fg,na.rm = T)
}

corM<-function(M){
  MM<-matrix(0, ncol=ncol(M),nrow=ncol(M))
  for(ii in 1:ncol(MM)){ for(jj in 1:ncol(MM))
    
  {MM[ii,jj]<-corp(M[,ii],M[,jj])
  }}
  colnames(MM)<-colnames(M)
  rownames(MM)<-colnames(M)
  MM
}

# Cargo la paleta de colores

paleta_divergente <- colorRampPalette(c("#B94E63", "#FFF0F5", "#B94E63"))(200)

