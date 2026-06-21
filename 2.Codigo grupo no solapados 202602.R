set.seed(pi)
#############
## Elijo aquellas funciones que tienen Ng genes asociados
Ng <- 30
NN <- 4
nivel_elegido <- 5 # <- CAMBIÁ ESTE NIVEL cuando quieras

# 1) Filtrar por nivel
Anot_nivel <- AnotacionesE %>%
  filter(depth_min == nivel_elegido)

xx <- table(Anot_nivel$GOasociado)
yy <- names(xx[xx %in% Ng])
yy <- sample(yy, length(yy))

yy <- sample(yy, NN)

PP <- list()
for (jj in 1:NN) {
  GO <- AnotacionesE %>%
    filter(GOasociado == yy[jj])
  
  AA <- transcD %>%
    filter(Identificador %in% GO$Identificador)
  
  AA$Identificador2 <- jj
  AA$GOasociado <- yy[jj]
  PP[[jj]] <- AA
}

datos <- do.call(rbind, PP)

asas<-datos$Identificador[duplicated(datos$Identificador)]
datos<-subset(datos,  ! Identificador %in% asas)
table(datos$GOasociado)


colores_grupo <- c(
  "1" = "#5DA5DA",  # azul
  "2" = "#F2C300",  # amarillo
  "3" = "#F15854",  # rojo
  "4" = "#60BD68"   # verde
)

tl_colors <- colores_grupo[as.character(datos$Identificador2)]

### Calcula la maxima correlacion de ventanas para todos los pares de tiempos, 
### moviendo el metodo de calculo

largg<- 24 ##########largo de la cadena de corr
larg<-largg-1
mmi<-"kendall"

dat <-t(datos[,2:25])
colnames(dat)<-datos[,26]
date1<-apply(dat, 2,est1)

# Usar en corrplot()
corrplot(as.matrix(corM(date1)),
         method = "shade",
         type = "full",
         diag = TRUE,
         tl.pos = "lt",   # importante: mostrar labels
         tl.cex = 0.8,
         tl.col = tl_colors,
         bg = "white",
         col = paleta_divergente,
         cl.pos = "n")  # Oculta la leyenda

colnames(date1)<-datos[,1]
############################################################3
#########grafos

umbral<-0.6
cor_mat <- corM(date1)
cor_mat<-cor_mat - diag(nrow(cor_mat))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
new_g <- graph_from_data_frame(only_sig, F)

# Data Preparation --------------------------------------------------------

#Load dataset

nodes<-data.frame(cbind(datos$Identificador,datos$Identificador))
colnames(nodes) <- c("id", "label")

nodes$group<-datos$Identificador2
#Edges
edges <- only_sig
colnames(edges) <- c("from", "to", "width")
edges$width<-4*(edges$width-umbral)
#Create graph for Louvain
graph <- graph_from_data_frame(edges, directed = FALSE)

#Louvain Comunity Detection
cluster <- cluster_louvain(graph)

cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- rownames(cluster_df)

#Create group column

nodes2 <- left_join(nodes, cluster_df, by = "label")
colnames(nodes2)[4] <- "group2"

#visNetwork(nodes, edges)
visNetwork(nodes2, edges, width = "100%") %>%
  visIgraphLayout() %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#FF8A9E",
      border = "#B94E63",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = TRUE, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#FF8A9E", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group2") %>% 
  visLayout(randomSeed = 11)

visNetwork(nodes, edges, width = "100%") %>%
  visIgraphLayout(layout = "layout_with_kk") %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#FF8A9E",
      border = "#B94E63",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = TRUE, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#FF8A9E", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group") %>% 
  visLayout(randomSeed = 11)


####################################################
###################################################

## Correlaciones parciales

 asas <- pcor(date1)
 a <- as.matrix(asas$estimate)
 a <- a %>% round(7)
 colnames(a) <- colnames(date1)
 rownames(a) <- colnames(date1)
 
 corrplot(a,
          method = "shade",
          type = "full",
          diag = TRUE,
          tl.pos = "lt",   # importante: mostrar labels
          tl.cex = 0.8,
          tl.col = tl_colors,
          bg = "white",
          col = paleta_divergente,
          cl.pos = "n")  # Oculta la leyenda

 
 umbral<-0.5
 cor_mat <- a
 cor_mat<- cor_mat - diag(nrow(cor_mat))
 cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
 cor_edge_list <- as_data_frame(cor_g, 'edges')
 only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
 new_g <- graph_from_data_frame(only_sig, F)
 
 # Data Preparation --------------------------------------------------------
 
 #Load dataset
 
 nodes<-data.frame(cbind(datos$Identificador,datos$Identificador))
 colnames(nodes) <- c("id", "label")
 
 nodes$group<-datos$Identificador2
 #Edges
 edges <- only_sig
 colnames(edges) <- c("from", "to", "width")
 edges$width<-4*(edges$width-umbral)
 #Create graph for Louvain
 graph <- graph_from_data_frame(edges, directed = FALSE)
 
 #Louvain Comunity Detection
 cluster <- cluster_louvain(graph)
 
 cluster_df <- data.frame(as.list(membership(cluster)))
 cluster_df <- as.data.frame(t(cluster_df))
 cluster_df$label <- rownames(cluster_df)
 
 #Create group column
 
 nodes2 <- left_join(nodes, cluster_df, by = "label")
 colnames(nodes2)[4] <- "group2"
 
 #visNetwork(nodes, edges)
 visNetwork(nodes2, edges, width = "100%") %>%
   visIgraphLayout(layout = "layout_with_kk") %>%
   visNodes(
     shape = "dot",
     color = list(
       background = "#FF8A9E",
       border = "#B94E63",
       highlight = "#FF8000"
     ),
     shadow = list(enabled = TRUE, size = 10)
   ) %>%
   visEdges(
     shadow = FALSE,
     color = list(color = "#FF8A9E", highlight = "#C62F4B")
   ) %>%
   visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
              selectedBy = "group2") %>% 
   visLayout(randomSeed = 11)
 
 visNetwork(nodes, edges, width = "100%") %>%
   visIgraphLayout(layout = "layout_with_kk") %>%
   visNodes(
     shape = "dot",
     color = list(
       background = "#FF8A9E",
       border = "#B94E63",
       highlight = "#FF8000"
     ),
     shadow = list(enabled = TRUE, size = 10)
   ) %>%
   visEdges(
     shadow = FALSE,
     color = list(color = "#FF8A9E", highlight = "#C62F4B")
   ) %>%
   visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
              selectedBy = "group") %>% 
   visLayout(randomSeed = 11)
 
###############################
##FOCI JOE 2022
combin1<-data.frame(t(combn(ncol(date1),2)))
rg<-list()
for( jj in 1:nrow(combin1) ){
  
  a1<-data.frame(setdiff(1:ncol(date1), combin1[jj,]))
  aa<-t(replicate(nrow(a1),as.numeric(combin1[jj,])))
  rg[[jj]]<-cbind(aa,a1)
}
combin<-do.call(rbind,rg)
names(combin)<-c("X1","X2","X3")


mm<-c()

for( jj in 1:nrow(combin) ){
  
  
  i<-date1[,combin[jj,1]]
  j<-date1[,combin[jj,2]]
  k<-date1[,combin[jj,3]]
  num<-abs(corp(i,j))-abs(corp(i,k)*corp(j,k))
  den<-sqrt(1 - corp(i,k)^2) *sqrt(1 - corp(j,k)^2) 
  mm[jj]<-num/den
}

combin$partial<-mm

library(reshape)
names(combin)

###correlaciones parciales
asas<-cast(combin, X1+X2~X3)

asas[is.na(asas)]<-1

correlacionespar<-asas[,! names(datos) %in% c("X1","X2")]
umbral2<-0.2
correlacionespar1<-correlacionespar<umbral2
aristas<-as.numeric(apply(correlacionespar1,1,sum)==0)

df<-data.frame(asas[,c(1,2)],aristas)

m1 <- matrix(0, nrow(datos), nrow(datos))
m1[as.matrix(df[1:2])] <- df$aristas


cor_mat <- m1+t(m1)
rownames(cor_mat )<-datos$Identificador2
colnames(cor_mat )<-datos$Identificador2
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list
new_g <- graph_from_data_frame(only_sig, F)

corrplot(cor_mat,
         method = "shade",
         type = "full",
         diag = TRUE,
         tl.pos = "lt",   # importante: mostrar labels
         tl.cex = 0.8,
         tl.col = tl_colors,
         bg = "white",
         col = paleta_divergente,
         cl.pos = "n")  # Oculta la leyenda

rownames(cor_mat )<-datos$Identificador
colnames(cor_mat )<-datos$Identificador

umbral<-0.5
cor_mat<- cor_mat - diag(nrow(cor_mat))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
new_g <- graph_from_data_frame(only_sig, F)

# Data Preparation --------------------------------------------------------

#Load dataset

nodes<-data.frame(cbind(datos$Identificador,datos$Identificador))
colnames(nodes) <- c("id", "label")

nodes$group<-datos$Identificador2
#Edges
edges <- only_sig
colnames(edges) <- c("from", "to", "width")
edges$width<-4*(edges$width-umbral)
#Create graph for Louvain
graph <- graph_from_data_frame(edges, directed = FALSE)

#Louvain Comunity Detection
cluster <- cluster_louvain(graph)

cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- rownames(cluster_df)

#Create group column

nodes2 <- left_join(nodes, cluster_df, by = "label")
colnames(nodes2)[4] <- "group2"

#visNetwork(nodes, edges)
visNetwork(nodes2, edges, width = "100%") %>%
  visIgraphLayout(layout = "layout_with_kk") %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#FF8A9E",
      border = "#B94E63",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = TRUE, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#FF8A9E", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group2") %>% 
  visLayout(randomSeed = 11)

visNetwork(nodes, edges, width = "100%") %>%
  visIgraphLayout(layout = "layout_with_kk") %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#FF8A9E",
      border = "#B94E63",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = TRUE, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#FF8A9E", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group") %>% 
  visLayout(randomSeed = 11)

############################################################################

cols_no_tiempo <- c("Identificador","Identificador2","GOasociado")
cols_tiempo <- setdiff(names(datos), cols_no_tiempo)

# escalar expresión entre 0 y 1 por gen
datos_scaled <- datos %>%
  rowwise() %>%
  mutate(across(all_of(cols_tiempo),
                ~ .x / max(c_across(all_of(cols_tiempo)), na.rm = TRUE))) %>%
  ungroup()

# formato largo
datos_long <- datos_scaled %>%
  pivot_longer(cols = all_of(cols_tiempo),
               names_to = "Tiempo",
               values_to = "Expresion")

datos_long$Tiempo <- factor(datos_long$Tiempo, levels = cols_tiempo)

# definir colores de los grupos
colores_grupo <- c(
  "1" = "#5DA5DA",
  "2" = "#F2C300",
  "3" = "#F15854",
  "4" = "#60BD68"
)

# gráfico
ggplot(datos_long,
       aes(x = Tiempo,
           y = Expresion,
           group = Identificador,
           color = factor(Identificador2))) +
  
  geom_line(alpha = 0.5, linewidth = 1) +
  
  facet_wrap(~ GOasociado, nrow = 2, ncol = 2) +
  
  scale_color_manual(values = colores_grupo) +
  
  labs(
    x = "Tiempo",
    y = "Expresión estandarizada",
    color = "Grupo"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
