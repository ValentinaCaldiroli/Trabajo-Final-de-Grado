library(ontologyIndex)
library(tidyverse)
library(writexl)
set.seed(1998)
############descargo bases
transcD <- read.csv2("C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/TFG/Developmental_transcriptome_Drosophila.csv")

Anotaciones <- read.delim("C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/TFG/Anotaciones_BP_Drosophila", header=FALSE)%>% 
  filter(V4 != "IEA") 
Anotaciones<-Anotaciones[, c(1,3)]
names(Anotaciones)<-c("Identificador", "GOasociado")

go <- get_ontology(
  "C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/scripts/go-basic.obo",
  extract_tags = "everything"
)


Anotaciones$Ontologia <- go$namespace[Anotaciones$GOasociado]
names(transcD)<- c("Identificador", names(transcD)[-1])
AnotacionesE<-subset(Anotaciones, Identificador %in% transcD$Identificador)
xx<-table(AnotacionesE$GOasociado) # GO:0006355 tiene 575 genes asociados
max(xx)
table(xx)  #######hay 1404 genes que solo expresaron en 1 funcion.
AnotacionesE$Ontologia <- as.character(AnotacionesE$Ontologia)
table(Anotaciones$Ontologia)

depth_min_to_root <- function(go, term, root) {
  # devuelve la distancia mínima desde term hasta root (0 si term==root)
  if (is.na(term) || !(term %in% names(go$parents))) return(NA_integer_)
  if (term == root) return(0L)
  
  visited <- character(0)
  frontier <- term
  depth <- 0L
  
  while (length(frontier) > 0) {
    depth <- depth + 1L
    # subir a padres (puede haber varios)
    parents <- unique(unlist(go$parents[frontier], use.names = FALSE))
    parents <- parents[!is.na(parents)]
    if (length(parents) == 0) return(NA_integer_)
    if (root %in% parents) return(depth)
    
    # evitar loops
    parents <- setdiff(parents, visited)
    visited <- union(visited, frontier)
    frontier <- parents
    
    # por seguridad
    if (depth > 2000) return(NA_integer_)
  }
  
  NA_integer_
}

get_depths <- function(go, terms, root) {
  terms <- as.character(terms)
  sapply(terms, function(t) depth_min_to_root(go, t, root))
}

# Ejemplo:
terms <- unique(AnotacionesE$GOasociado)
depths_bp <- get_depths(go, terms, root = "GO:0008150")

# agregar a tu tabla (por GO)
depth_df <- data.frame(GOasociado = names(depths_bp),
                       depth_min = as.integer(depths_bp),
                       stringsAsFactors = FALSE)

AnotacionesE <- merge(AnotacionesE, depth_df, by="GOasociado", all.x=TRUE)

write_xlsx(AnotacionesE, "C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/scripts/DatosSinTrasncripcion.xlsx") 
