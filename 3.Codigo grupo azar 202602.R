set.seed(pi)

Ng <- 1
NN <- 95
nivel_elegido <- 5 # <- CAMBIÁ ESTE NIVEL cuando quieras

# 1) Filtrar por nivel
Anot_nivel <- AnotacionesE %>%
  filter(depth_min == nivel_elegido)

xx_azar <- table(Anot_nivel$GOasociado)
yy_azar <- names(xx_azar[xx_azar %in% Ng])   

yy_azar <- sample(yy_azar, length(yy_azar))

yy_azar <- sample(yy_azar, NN)

PP_azar <- list()
for (jj in 1:NN) {
  GO <- Anot_nivel %>%
    filter(GOasociado == yy_azar[jj])
  
  AA_azar <- transcD %>%
    filter(Identificador %in% GO$Identificador)
  
  AA_azar$Identificador2 <- jj
  AA_azar$GOasociado <- yy_azar[jj]
  PP_azar[[jj]] <- AA_azar
}

azar <- do.call(rbind, PP_azar)

asas_az<-azar$Identificador[duplicated(azar$Identificador)]
azar<-subset(azar,  ! Identificador %in% asas_az)
table(azar$GOasociado)


az <-t(azar[,2:25])
colnames(az)<-azar[,1]

azar1 <- apply(az, 2,est1)

largg<- 24 ##########largo de la cadena de corr
larg<-largg-1
mmi<-"kendall"

# Usar en corrplot()
corrplot(as.matrix(corM(azar1)),
         method = "shade",
         type = "full",
         diag = TRUE,
         tl.pos = "n",
         bg = "white",
         col = paleta_divergente,
         cl.pos = "n")  # Oculta la leyenda

umbral<-0.6
cor_mat <- corM(azar1)
cor_mat<-cor_mat - diag(nrow(cor_mat))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
new_g <- graph_from_data_frame(only_sig, F)

# Data Preparation --------------------------------------------------------

#Load dataset

nodes<-data.frame(cbind(azar$Identificador,azar$Identificador))
colnames(nodes) <- c("id", "label")

nodes$group<-azar$Identificador2
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

####################################################
###################################################

## Correlaciones parciales

asas <- pcor(azar1)
a <- as.matrix(asas$estimate)
a <- a %>% round(7)
colnames(a) <- colnames(azar1)
rownames(a) <- colnames(azar1)

corrplot(a,
         method = "shade",
         type = "full",
         diag = TRUE,
         tl.pos = "n",
         bg = "white",
         col = paleta_divergente,
         cl.pos = "n")  # Oculta la leyenda


umbral<-0.6
cor_mat <- a
cor_mat<- cor_mat - diag(nrow(cor_mat))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
new_g <- graph_from_data_frame(only_sig, F)

# Data Preparation --------------------------------------------------------

#Load dataset

nodes<-data.frame(cbind(azar$Identificador,azar$Identificador))
colnames(nodes) <- c("id", "label")

nodes$group<-azar$Identificador2
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

###############################
##FOCI JOE 2022
combin1<-data.frame(t(combn(ncol(azar1),2)))
rg<-list()
for( jj in 1:nrow(combin1) ){
  
  a1<-data.frame(setdiff(1:ncol(azar1), combin1[jj,]))
  aa<-t(replicate(nrow(a1),as.numeric(combin1[jj,])))
  rg[[jj]]<-cbind(aa,a1)
}
combin<-do.call(rbind,rg)
names(combin)<-c("X1","X2","X3")


mm<-c()

for( jj in 1:nrow(combin) ){
  
  
  i<-azar1[,combin[jj,1]]
  j<-azar1[,combin[jj,2]]
  k<-azar1[,combin[jj,3]]
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

correlacionespar<-asas[,! names(asas) %in% c("X1","X2")]
umbral2<-0.2
correlacionespar1<-correlacionespar<umbral2
aristas<-as.numeric(apply(correlacionespar1,1,sum)==0)

df<-data.frame(asas[,c(1,2)],aristas)

m1 <- matrix(0, ncol(azar1), ncol(azar1))
m1[as.matrix(df[1:2])] <- df$aristas

cor_mat <- m1+t(m1)
rownames(cor_mat )<-azar$Identificador
colnames(cor_mat )<-azar$Identificador
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list
new_g <- graph_from_data_frame(only_sig, F)

corrplot(cor_mat,
         method = "shade",
         type = "full",
         diag = TRUE,
         tl.pos = "n",
         bg = "white",
         col = paleta_divergente,
         cl.pos = "n")


umbral<-0.6
cor_mat<- cor_mat - diag(nrow(cor_mat))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = 'correlation')
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[cor_edge_list$correlation > umbral, ]
new_g <- graph_from_data_frame(only_sig, F)

# Data Preparation --------------------------------------------------------

#Load dataset

nodes<-data.frame(cbind(azar$Identificador,azar$Identificador))
colnames(nodes) <- c("id", "label")

nodes$group<-azar$Identificador2
#Edges
edges <- only_sig
colnames(edges) <- c("from", "to", "width")
edges$title <- paste0("Correlación: ", round(edges$width, 3))
#edges$width <- 4*(edges$correlation - umbral)
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

############################################################################
cols_no_tiempo <- c("Identificador","Identificador2","GOasociado")
cols_tiempo <- setdiff(names(azar), cols_no_tiempo)

# escalar expresión entre 0 y 1 por gen
azar_scaled <- azar %>%
  rowwise() %>%
  mutate(across(all_of(cols_tiempo),
                ~ .x / max(c_across(all_of(cols_tiempo)), na.rm = TRUE))) %>%
  ungroup()

# formato largo
azar_long <- azar_scaled %>%
  pivot_longer(cols = all_of(cols_tiempo),
               names_to = "Tiempo",
               values_to = "Expresion")

azar_long$Tiempo <- factor(azar_long$Tiempo, levels = cols_tiempo)

# gráfico final
ggplot(azar_long,
       aes(x = Tiempo,
           y = Expresion,
           group = Identificador)) +
  
  geom_line(color = "#F4A6B8", alpha = 0.5, linewidth = 1) +
  
  # línea promedio (opcional pero recomendado)
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               color = "black",
               linewidth = 1.2) +
  
  labs(
    x = "Tiempo",
    y = "Expresión estandarizada"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

