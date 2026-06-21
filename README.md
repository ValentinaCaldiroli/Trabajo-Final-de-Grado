

# Detección de comunidades para inferir relaciones funcionales entre genes a partir de perfiles temporales de expresión en *Drosophila melanogaster*

## Trabajo Final de Grado – Licenciatura en Estadística

**Autora:** Valentina Caldiroli
**Tutores:** Leonardo Moreno y Flavio Pazos

---

## Resumen

Este repositorio contiene el código, resultados y documentación asociados al Trabajo Final de Grado de la Licenciatura en Estadística.

El objetivo del trabajo es evaluar la capacidad de diferentes medidas de dependencia para recuperar estructura funcional biológica a partir de perfiles temporales de expresión génica en *Drosophila melanogaster*.

Para ello se construyen redes de coexpresión génica utilizando distintas medidas de asociación entre genes y posteriormente se aplican algoritmos de detección de comunidades. Finalmente, se evalúa la coherencia biológica de las comunidades obtenidas mediante análisis de enriquecimiento funcional basado en Gene Ontology (GO).

---

## Objetivos

* Analizar perfiles temporales de expresión génica durante el desarrollo de *Drosophila melanogaster*.
* Construir redes de asociación entre genes utilizando distintas medidas de dependencia.
* Detectar comunidades mediante el algoritmo de Louvain.
* Evaluar la correspondencia entre las comunidades detectadas y la funcionalidad biológica de los genes.
* Comparar el desempeño de distintos enfoques estadísticos para la inferencia funcional.

---

## Metodología

El análisis se basa en la comparación de distintas medidas de dependencia:

* Correlación de Pearson
* Correlación de Spearman
* Correlación de Kendall
* Correlaciones parciales
* FOCI (Feature Ordering by Conditional Independence)

A partir de las matrices de asociación obtenidas:

1. Se construyen grafos de coexpresión génica.
2. Se aplican distintos umbrales de asociación.
3. Se detectan comunidades mediante el algoritmo de Louvain.
4. Se evalúa el enriquecimiento funcional de las comunidades utilizando términos Gene Ontology.

---

## Datos

Se utilizan datos de expresión génica temporal de *Drosophila melanogaster* obtenidos durante distintas etapas del desarrollo:

* Embrión
* Larva
* Pupa
* Adulto

Los genes fueron asociados a términos Gene Ontology (GO) pertenecientes a la ontología de Procesos Biológicos (Biological Process).

---

## Principales herramientas utilizadas

### Lenguajes

* R
* LaTeX

### Librerías de R

* tidyverse
* dplyr
* ggplot2
* igraph
* corrplot
* ppcor
* FOCI
* visNetwork
* knitr
* kableExtra

---

## Resultados principales

Los resultados obtenidos muestran que:

* Los perfiles temporales de expresión contienen información relevante sobre la funcionalidad génica.
* Las correlaciones marginales permiten recuperar comunidades con mayor coherencia funcional que las medidas de dependencia condicional consideradas.
* La correlación de Pearson fue el método que produjo, en general, las comunidades con mayor enriquecimiento funcional.
* Los conjuntos de genes seleccionados aleatoriamente no presentaron patrones funcionales comparables a los observados en los grupos definidos por Gene Ontology.

---

## Reproducibilidad

Todo el código utilizado para generar las figuras, tablas y resultados presentados en la memoria se encuentra disponible en este repositorio.

---

## Contacto

**Valentina Caldiroli**

Licenciatura en Estadística
Universidad de la República (Uruguay)
