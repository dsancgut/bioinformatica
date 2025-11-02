# Daniela Sánchez Gutiérrez_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)
install.packages("readr")
library(readr)
datos_biomed <- read_csv("datos_biomed.csv")
# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
head(datos_biomed)
summary(datos_biomed)
dim(datos_biomed)
str(datos_biomed)
unique(datos_biomed$Tratamiento)
# Hay 5 variables(cada una de las columnas) y hay 3 tipos de tratamientos:FarmacoA,FarmacoB,Placebo.
library(ggplot2)

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
library(tidyr)
library(ggplot2)

# Pasar de formato ancho a largo (long)
datos_largos <- pivot_longer(datos_biomed,
                             cols = c(Glucosa, Presion, Colesterol),
                             names_to = "Variable",
                             values_to = "Valor")

ggplot(datos_largos, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Distribución de variables biomédicas por tratamiento",
       x = "Tratamiento",
       y = "Valor") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Realiza un violin plot (investiga qué es). (1 pt)

library(ggplot2)
library(tidyr)

# Leer los datos
datos <- read.csv("datos_biomed.csv", header = TRUE, sep = ",")

# Pasar los datos a formato largo
datos_largo <- datos |>
  pivot_longer(cols = c(Glucosa, Presion, Colesterol),
               names_to = "Variable",
               values_to = "Valor")

ggplot(datos_largo, aes(x = Tratamiento, y = Valor)) +
  geom_violin(fill = "#40E0D0", color = "white") +       
  geom_boxplot(width = 0.1, fill = "#00008B", color = "white") +  
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Distribución de variables biomédicas por tratamiento",
       x = "Tratamiento",
       y = "Valor") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )



# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)

datos <- read.csv("datos_biomed.csv", header = TRUE, sep = ",")

# Crear el gráfico de dispersión
plot(datos$Glucosa, datos$Presion,
     col = c("red", "blue", "green")[as.factor(datos$Tratamiento)],
     pch = 19,                         # tipo de punto sólido
     xlab = "Glucosa (mg/dL)",
     ylab = "Presión (mmHg)",
     main = "Relación entre Glucosa y Presión por tratamiento")

# Agregar la leyenda
legend("bottomright", 
       legend = levels(as.factor(datos$Tratamiento)),  # nombres de grupos
       col = c("red", "blue", "green"),                # mismos colores que en el gráfico
       pch = 19,                                       # tipo de punto
       title = "Tratamiento")

# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)

library(ggplot2)

# Leer los datos
datos <- read.csv("datos_biomed.csv", header = TRUE, sep = ",")

# Crear el gráfico con facet grid
ggplot(datos, aes(x = Presion, y = Colesterol)) +
  geom_point(aes(color = Tratamiento), size = 3) +   # puntos de color por tratamiento
  facet_grid(~ Tratamiento) +                        # un panel por tratamiento (en columnas)
  labs(title = "Relación entre Colesterol y Presión por tratamiento",
       x = "Presión (mmHg)",
       y = "Colesterol (mg/dL)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5))



# 7. Realiza un histogramas para cada variable. (0.5 pts)

# Asignar colores específicos a cada variable
colores <- c("Glucosa" = "yellow", "Presion" = "purple", "Colesterol" = "orange")

# Crear los histogramas
ggplot(datos_largo, aes(x = Valor, fill = Variable)) +
  geom_histogram(color = "black", bins = 15) +
  facet_wrap(~ Variable, scales = "free", ncol = 1) +  # un panel por variable, en columnas
  scale_fill_manual(values = colores) +
  labs(title = "Histogramas de variables biomédicas",
       x = "Valor",
       y = "Frecuencia") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # no hace falta leyenda, el nombre ya está en cada panel
  )


# 8.Crea un factor a partir del tratamiento. Investiga factor(). (1 pt)
datos$Tratamiento <- factor(datos$Tratamiento)

# Verificar el cambio
str(datos$Tratamiento)
levels(datos$Tratamiento)

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)

#Verificar los nombres de las columnas
names(datos)

#Calcular la media de glucosa por tratamiento
media_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = mean)

#Calcular la desviación estándar de glucosa por tratamiento
sd_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = sd)

#Unir ambos resultados en una sola tabla
resumen_glucosa <- merge(media_glucosa, sd_glucosa, by = "Tratamiento", suffixes = c("_media", "_sd"))

#Mostrar el resultado final
print(resumen_glucosa)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
head(datos)
names(datos)                # para confirmar el nombre de la columna
unique(datos$Tratamiento)  # para ver qué tratamientos hay exactamente

#Extraer los datos por tratamiento
farmacoA <- subset(datos, Tratamiento == "FarmacoA")
farmacoA
farmacoB <- subset(datos, Tratamiento == "FarmacoB")
farmacoB
placebo  <- subset(datos, Tratamiento == "Placebo")
placebo

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
#Test de normalidad (Shapiro-Wilk) para cada tratamiento
shapiro.test(farmacoA$Glucosa)# p-value= 0.3008 > 0.05, datos normales
shapiro.test(farmacoB$Glucosa)# p-value= 0.411 > 0.05, datos normales
shapiro.test(placebo$Glucosa) # p-value= 0,402 > 0.05, datos normales

#Hacer comparativa de medias
anova_resultado <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_resultado)
# p-value= 0.1 > 0.05, el resultado no es significativo, es decir no se puede rechazar la hipotesis nula a igualdad de medias.

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)



