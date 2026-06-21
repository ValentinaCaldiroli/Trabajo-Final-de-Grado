library(ontologyIndex)
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(kableExtra)

# =========================
# 1. Columnas de expresión
# =========================

cols_tiempo <- setdiff(names(transcD), "Identificador")

# =========================
# 2. Resumen de anotaciones GO por gen
# =========================

anot_por_gen <- AnotacionesE %>%
  distinct(Identificador, GOasociado) %>%
  group_by(Identificador) %>%
  summarise(
    n_GO = n_distinct(GOasociado),
    GO_asociados = paste(unique(GOasociado), collapse = ", "),
    .groups = "drop"
  )

# =========================
# 3. Unir expresión + anotaciones
# =========================

datos_total <- transcD %>%
  inner_join(anot_por_gen, by = "Identificador")

# =========================
# 4. Calcular descriptivas por gen
# =========================

descriptiva_genes <- datos_total %>%
  rowwise() %>%
  mutate(
    media_exp = mean(c_across(all_of(cols_tiempo)), na.rm = TRUE),
    sd_exp = sd(c_across(all_of(cols_tiempo)), na.rm = TRUE),
    min_exp = min(c_across(all_of(cols_tiempo)), na.rm = TRUE),
    max_exp = max(c_across(all_of(cols_tiempo)), na.rm = TRUE),
    rango_exp = max_exp - min_exp,
    cv_exp = sd_exp / media_exp
  ) %>%
  ungroup()

# =========================
# 5. Top 10 genes con mayor expresión promedio
# =========================

top10_expresion <- descriptiva_genes %>%
  arrange(desc(media_exp)) %>%
  slice(1:20) %>%
  dplyr::select(
    Identificador,
    n_GO,
    media_exp,
    sd_exp,
    min_exp,
    max_exp,
    rango_exp
  ) %>%
  dplyr::mutate(
    media_exp = round(media_exp, 3),
    sd_exp = round(sd_exp, 3),
    min_exp = round(min_exp, 3),
    max_exp = round(max_exp, 3),
    rango_exp = round(rango_exp, 3)
  )

##################################################################################
go <- get_ontology("C:/Users/mrocc/OneDrive/Desktop/Facultad/TFG/scripts/go-basic.obo", extract_tags = "everything")

terms0 <- c("GO:0042981", "GO:0007619", "GO:0035293", "GO:0007169")

nombres_GO <- c(
  "GO:0042981" = "Regulación del proceso apoptótico",
  "GO:0007619" = "Comportamiento de cortejo",
  "GO:0035293" = "Regulación positiva de la transcripción",
  "GO:0007169" = "Señalización por receptor tirosina quinasa"
)

# Obtener ancestros de un término GO
get_ancestors_all <- function(term, go) {
  visited <- term
  frontier <- term
  
  repeat {
    parents <- unlist(go$parents[frontier], use.names = FALSE)
    parents <- unique(parents[!is.na(parents)])
    parents <- setdiff(parents, visited)
    
    if (length(parents) == 0) break
    
    visited <- c(visited, parents)
    frontier <- parents
  }
  
  unique(visited)
}

# Distancia desde un término hasta un ancestro
dist_to_ancestor <- function(term, ancestor, go) {
  if (term == ancestor) return(0)
  
  frontier <- term
  visited <- term
  dist <- 0
  
  repeat {
    dist <- dist + 1
    parents <- unlist(go$parents[frontier], use.names = FALSE)
    parents <- unique(parents[!is.na(parents)])
    
    if (ancestor %in% parents) return(dist)
    
    parents <- setdiff(parents, visited)
    if (length(parents) == 0) return(Inf)
    
    visited <- c(visited, parents)
    frontier <- parents
  }
}

# Profundidad mínima desde raíz aproximada
get_depth <- function(term, go) {
  ancestors <- get_ancestors_all(term, go)
  length(ancestors)
}

pares_GO <- expand.grid(GO1 = terms0, GO2 = terms0, stringsAsFactors = FALSE)

resultado_dist <- pares_GO %>%
  rowwise() %>%
  mutate(
    anc_GO1 = list(get_ancestors_all(GO1, go)),
    anc_GO2 = list(get_ancestors_all(GO2, go)),
    anc_comunes = list(intersect(anc_GO1, anc_GO2)),
    
    ancestro_comun = {
      comunes <- anc_comunes
      
      if (length(comunes) == 0) {
        NA_character_
      } else {
        profundidades <- sapply(comunes, get_depth, go = go)
        comunes[which.max(profundidades)]
      }
    },
    
    distancia_GO1 = dist_to_ancestor(GO1, ancestro_comun, go),
    distancia_GO2 = dist_to_ancestor(GO2, ancestro_comun, go),
    distancia_total = distancia_GO1 + distancia_GO2,
    nombre_ancestro = go$name[ancestro_comun]
  ) %>%
  ungroup() %>%
  dplyr::select(
    GO1, GO2, ancestro_comun, nombre_ancestro,
    distancia_GO1, distancia_GO2, distancia_total
  )


mat_dist <- resultado_dist %>%
  dplyr::select(GO1, GO2, distancia_total) %>%
  pivot_wider(names_from = GO2, values_from = distancia_total)

mat_dist

df_plot <- resultado_dist %>%
  dplyr::mutate(
    GO1_label = nombres_GO[GO1],
    GO2_label = nombres_GO[GO2]
  )

ggplot(df_plot, aes(x = GO1_label, y = GO2_label, fill = distancia_total)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = distancia_total), size = 5) +
  scale_fill_gradient(low = "#FADADD", high = "#C2185B") +
  labs(
    x = "",
    y = "",
    fill = "Distancia",
    title = "Distancia ontológica entre términos GO"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

tabla_GO_distancias <- resultado_dist %>%
  filter(GO1 < GO2) %>%
  dplyr::mutate(
    Termino_1 = nombres_GO[GO1],
    Termino_2 = nombres_GO[GO2]
  ) %>%
  dplyr::select(
    GO1, Termino_1,
    GO2, Termino_2,
    ancestro_comun,
    nombre_ancestro,
    distancia_GO1,
    distancia_GO2,
    distancia_total
  )

tabla_GO_distancias


library(kableExtra)
library(dplyr)

tabla_latex <- tabla_GO_distancias %>%
  mutate(
    GO1 = as.character(GO1),
    GO2 = as.character(GO2),
    ancestro_comun = as.character(ancestro_comun)
  )

kable(tabla_latex,
      format = "latex",
      booktabs = TRUE,
      align = "c",
      caption = "Distancias ontológicas entre pares de términos GO, calculadas a partir del ancestro común más específico en la jerarquía de Gene Ontology.") %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )


library(dplyr)
library(tidyr)
library(ggplot2)

df_long <- as.data.frame(date1) %>%
  mutate(Tiempo = rownames(.)) %>%
  pivot_longer(
    cols = -Tiempo,
    names_to = "Identificador",
    values_to = "Expresion"
  )

df_long <- df_long %>%
  left_join(
    datos[, c("Identificador", "GOasociado")],
    by = "Identificador"
  )

resumen_GO <- df_long %>%
  group_by(GOasociado, Tiempo) %>%
  summarise(
    Expresion = median(Expresion, na.rm = TRUE),
    .groups = "drop"
  )

nombres_GO <- c(
  "GO:0042981" = "Apoptosis",
  "GO:0007619" = "Cortejo",
  "GO:0035293" = "Transcripción",
  "GO:0007169" = "Señalización RTK"
)

resumen_GO$GO_label <- nombres_GO[resumen_GO$GOasociado]

resumen_GO$Tiempo <- factor(resumen_GO$Tiempo, levels = rownames(date1))

library(ggplot2)

ggplot(resumen_GO, aes(x = Tiempo, y = GO_label, fill = Expresion)) +
  
  geom_tile(color = "white") +
  
  scale_fill_gradient(
    low = "#FADADD",
    high = "#C2185B"
  ) +
  
  labs(
    x = "Tiempo",
    y = "Término GO",
    fill = "Expresión\n(mediana)",
    title = "Perfil temporal promedio por función biológica"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
