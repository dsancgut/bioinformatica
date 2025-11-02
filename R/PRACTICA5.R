#############################################################################
#
# PRACTICA R
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")
a
# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) # significan los datos de nuestra tabla, tiene 12488 filas que son genes en 6 columnas
head(data)# 6 primeras filas de la tabla
tail(data) # 6 ultimas filas de la tabla

# Hacemos un primer histograma para explorar los datos
hist(data, col = "gray", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo para que se coloquen mejor 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data2 = log2(data)
# para organizar y que tus datos queden más bonitos
hist(data2, col = "gray", main="GSE5583 (log2) - Histogram")


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# las cajas son nuestras seis muestras, las azules wildtipe y naranjas knockout, en col hemos puesto los colores en el orden de las columnas de la tabla.
#main se refiere a l título de nuestro boxplot y las=2 sirve para poner la orientación de las etiquetas en peroendicular al eje
# la linea del medio es la mediana, cuando ponemos main el titulo, las 2 para que pongan los nombres en vertical
# ¿Qué es un boxplot?
# Es una representación gráfica de la distribución de nuestros datos, con la mediana, el Q1,Q2 y Q3.
boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots",las=2)
	



# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
# de los valores de expresión. ¿Es correcta la separación?
# es correcta porque nos separa las muestras de knockout y wildtype, deberían estar separados.
# sirve pra confirmar que el experimento está bien, lo piden en publicaciones
hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################
head(data)
# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado? una matriz
wt <- data[,1:3] #delante de la coma van las filas y detrás las columnas, seleccionamos todas las filas y solo las tres primeras cplumnas es decir los genes wildtype
ko <- data[,4:6] #hay que quitarle el dos poqrue si no nos sale sobre los datos logarítmicos
class(wt)
head(wt)
head(ko)
# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)#para calcular la media con apply, objeto, 1 para todas las filas y un 3 para las columnas y mean porque es la media
ko.mean = apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)

# ¿Cuál es la media más alta? 37460.5
limit = max(wt.mean, ko.mean) 
limit

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Añadir una línea diagonal con abline nos indica una correlación
abline(0, 1, col = "red")

# ¿Eres capaz de añadirle un grid? si, nos permite generar una cuadricula, para que me ejecute tengo que quitar los comentarios
grid()
#abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
#abline(1, 2, col = "red")     # línea y = 2x + 1
#abline(h = 2, col = "green")  # línea y = 2
#abline(v = 3, col = "violet") # línea x = 3

# Calculamos la diferencia entre las medias de las condiciones, si tenemos sobreexpreión en ko el resultdo será positivo y si es mayor el wt será positivo
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias hay poco gens que tengas expresión diferencial
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? porque para hacer los análisis siempre tienes que utilizar los datos brutos
# ¿Cuántas valores tiene cada muestra? 12488, un p-value para cada gen
pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)#hemos creado pvalue vacio y hemos creado un bucle para cada fila de mis datos
#dentro de cada fila me  va a sacar un vector x(3valors de wt) y un vector y( 3valors de ko)
# tendríamos que comprobar si es una distribución normal, si lo es utilizamos ´-student y si no el wilcox-test
# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)
# un pvalue para cada gen 12488
# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10? coloca y organiza nuestros datos de forma bonita.
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico? en las x la diferencia de medias, el corte en 0 no hay cambio, cada punto es un gen
# los sobreexpresados en ko estarán en ngativo por lo que estarán hacia la izquierda, la línea verde es el corte de p-value, 0,05 en el 2, los significativos por encima de dos
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios? 426 genes
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? porque los sobreexpresados en wt estarán reprimidos en ko y los ko estarán reprimidos en wt 
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")
#rojo genes ko y azul wt, y los negros son los que no pasan el p-value los que no son significativo


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# filtered: es la matriz de datos que se va a reperesentar, rowv:controla cómo se ordenan las filas, colv:igual que rowv pero para las columnas,cexCol: ajusta el tamaño del texto de las etiquetas de las columnas y labRow:FALSE, oculta los nombres de las filas 
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)

heatmap(filtered)
# deberían salir wildtipes en un lado y ko en otro, en rojo oscuro tiene sobreexpresión, y amarillo reprimido

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")
# a la derecha los nombres de los genes

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE) 

