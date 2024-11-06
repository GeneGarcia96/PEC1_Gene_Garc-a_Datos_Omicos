# Cargar las bibliotecas necesarias
library(SummarizedExperiment)

# Cargar los datos
data_values <- read.csv("DataValues_S013.csv", row.names = 1)
sample_info <- read.csv("DataInfo_S013.csv", row.names = 1)
metabolite_info <- read.csv("AAInformation_S006.csv", row.names = 1)

# Asegurarse de que las columnas de sample_info coincidan con las columnas de data_values
sample_info <- sample_info[match(colnames(data_values), rownames(sample_info)), ]
# Asegurarse de que las filas de metabolite_info coincidan con las filas de data_values
metabolite_info <- metabolite_info[match(rownames(data_values), rownames(metabolite_info)), ]

# Verificar dimensiones de los datos
print(dim(data_values))         # Dimensiones de la matriz de abundancia
print(dim(sample_info))         # Dimensiones de la información de muestras
print(dim(metabolite_info))     # Dimensiones de la información de metabolitos


# Crear el objeto SummarizedExperiment si las dimensiones coinciden
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(data_values)),
  colData = sample_info,
  rowData = metabolite_info
)

# Revisar el objeto
se
colData(se)  
assays(se)[[1]][1:5, 1:5]

rowData(se[1:5])      

metadata(se)

assay(se)
# Seleccionar solo las columnas de metabolitos que contienen valores numéricos
metabolites_data <- assay(se)[, 6:ncol(assay(se))]

# Verificar si todas las columnas son numéricas y convertir si es necesario
metabolites_data <- as.data.frame(metabolites_data)  # Convertir a data.frame si es necesario
metabolites_data <- data.matrix(metabolites_data)    # Convertir todas las columnas a matriz numérica


##Primeros analisis 

# Graficar el histograma de abundancia de metabolitos
hist(
  as.vector(metabolites_data),
  main = "Distribución de Abundancias de Metabolitos",
  xlab = "Abundancia",
  breaks = 30
)

# Boxplot de las abundancias de los metabolitos por muestra
boxplot(
  metabolites_data,
  main = "Boxplot de Abundancias por Muestra",
  xlab = "Muestras",
  ylab = "Abundancia",
  las = 2
)

cor_matrix <- cor(t(metabolites_data), use = "complete.obs")
heatmap(
  cor_matrix,
  main = "Matriz de Correlación entre Metabolitos",
  xlab = "Metabolitos",
  ylab = "Metabolitos"
)


## Filtrado de datos para su analisis 

## Mapa de calor 

# Reemplazar NA/NaN/Inf por la mediana de cada metabolito (columna)
# Convertir metabolites_data a una matriz numérica limpia

metabolites_data_clean <- as.matrix(metabolites_data)

# Reemplazar todos los valores NA, NaN, e Inf con la mediana de cada columna
metabolites_data_clean <- apply(metabolites_data_clean, 2, function(col) {
  col[is.na(col) | is.nan(col) | is.infinite(col)] <- median(col, na.rm = TRUE)
  return(col)
})

metabolite_clean <- metabolites_data_clean[, colSums(is.na(metabolites_data_clean)) == 0]

#Heatmap 
heatmap(metabolite_clean, main = "Mapa de Calor de Clustering", 
        xlab = "Muestras", ylab = "Metabolitos")

# Realizar PCA
pca <- prcomp(t(metabolite_clean), scale. = TRUE)

# Graficar los primeros dos componentes
plot(pca$x[, 1], pca$x[, 2], 
     xlab = "PC1", ylab = "PC2", 
     main = "Análisis de Componentes Principales (PCA)")

#Cluster 

dist_matrix <- dist(metabolite_clean, method = "euclidean")

# 2. Realizar el clustering jerárquico
hc <- hclust(dist_matrix, method = "ward.D2")

# 3. Graficar el dendrograma
plot(hc, main = "Clustering Jerárquico de Muestras", xlab = "", sub = "", ylab = "Distancia")


## Analisis de coeficiente 

cv <- apply(metabolite_clean, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
barplot(cv, main = "Coeficiente de Variación de Metabolitos", ylab = "CV", xlab = "Metabolitos")




