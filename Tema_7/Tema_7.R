# ============================================================================================================================================

# Ejercicio 1.

# ============================================================================================================================================

# Cargamos los datos

data(iris)
head(iris)
summary(iris)
str(iris)

### Utilización del método JERARQUICO

# cargamos los paquetes que vamos a necesitar

library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)

# Paso 1
datos <- iris[,-5]

# Normalizamos los datos
datos.norm <- decostand(datos, "normalize") 
summary(datos.norm)

# PASO 2. Elección de la técnica cluster. ===============================
datos.ch <- vegdist(datos.norm, "euc") #calculamos la matriz de distancias euclídea

# método de conglomeración. Usaremos distintos métodos para ver luego cual va mejor.
par(mfrow=c(2,3)) #para disponer juntos los gráficos
datos.ch.single <- hclust(datos.ch, method="single") # Método "single" (método del vecino más cercano o mínima distancia).  
plot(datos.ch.single) # graficamos el dendrograma
datos.ch.complete <- hclust(datos.ch, method="complete") # Método "complete" (método del vecino más lejano o máxima distancia).
plot(datos.ch.complete) 
datos.ch.UPGMA <- hclust(datos.ch, method="average") # Método "UPGMA" o "Average" (método del promedio)
plot(datos.ch.UPGMA)
datos.ch.centroid <- hclust(datos.ch, method="centroid") # Método "centroid" (método del centroide).  
plot(datos.ch.centroid)
datos.ch.ward <- hclust(datos.ch, method="ward.D") # Método de mínima varianza de Ward
plot(datos.ch.ward)
datos.ch.median<- hclust(datos.ch, method="median") # Método "Median" (método del promedio)
plot(datos.ch.median)

# Determinar qué método de conglomeración es el mejor.  -----------------------#            
#   Correlación cofenética. mayor coeficiente = mejor modelo (menor distorsión)
datos.ch.single.coph <- cophenetic(datos.ch.single) # Single
cor(datos.ch, datos.ch.single.coph)
datos.ch.comp.coph <- cophenetic(datos.ch.complete) # Complete
cor(datos.ch, datos.ch.comp.coph)
datos.ch.UPGMA.coph <- cophenetic(datos.ch.UPGMA) # Average
cor(datos.ch, datos.ch.UPGMA.coph)
datos.ch.ward.coph <- cophenetic(datos.ch.ward) # Ward 
cor(datos.ch, datos.ch.ward.coph)
datos.ch.median.coph <- cophenetic(datos.ch.median) # Median 
cor(datos.ch, datos.ch.median.coph)
# La de Promedio es la que mayor valor arroja. Sin embargo, el que mejor representa el dendrograma es Ward
# El resto parece no afinar en las terminaciones del dicha visualización. Además es el segundo con valor más alto.

# Representación gráfica de las correlaciones cofenéticas: Diagrama de Shepard ----
# relaciones entre la matriz de distancias original y la matriz de distancias cofenética 
par(mfrow=c(2,2))
plot(datos.ch, datos.ch.single.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Single linkage", paste("Cophenetic correlation =",
                                    round(cor(datos.ch, datos.ch.single.coph),3))))
abline(0,1);  lines(lowess(datos.ch, datos.ch.single.coph), col="red")

plot(datos.ch, datos.ch.comp.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Complete linkage", paste("Cophenetic correlation =",
                                      round(cor(datos.ch, datos.ch.comp.coph),3))))
abline(0,1);  lines(lowess(datos.ch, datos.ch.comp.coph), col="red")

plot(datos.ch, datos.ch.UPGMA.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("UPGMA", paste("Cophenetic correlation =",
                           round(cor(datos.ch, datos.ch.UPGMA.coph),3))))
abline(0,1);  lines(lowess(datos.ch, datos.ch.UPGMA.coph), col="red")

plot(datos.ch, datos.ch.ward.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), 
     ylim=c(0,max(datos.ch.ward$height)),
     main=c("Ward clustering", paste("Cophenetic correlation =",
                                     round(cor(datos.ch, datos.ch.ward.coph),3))))
abline(0,1);  lines(lowess(datos.ch, datos.ch.ward.coph), col="red")

plot(datos.ch, datos.ch.median.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Median linkage", paste("Cophenetic correlation =",
                                      round(cor(datos.ch, datos.ch.median.coph),3))))
abline(0,1);  lines(lowess(datos.ch, datos.ch.median.coph), col="red")

# El más centrado corresponde con la de Promedio UPGMA.

# Distancia de Gower (1983) ------
(gow.dist.single <- sum((datos.ch-datos.ch.single.coph)^2))
(gow.dist.comp <- sum((datos.ch-datos.ch.comp.coph)^2))
(gow.dist.UPGMA <- sum((datos.ch-datos.ch.UPGMA.coph)^2))
(gow.dist.ward <- sum((datos.ch-datos.ch.ward.coph)^2))
# De nuevo el que menor valor arroja es la de Promedio.

# PASO 3: Elección del número de conglomerados finales y validación. =======================

# Gráfico de nivel de fusión ------------
plot(datos.ch.UPGMA$height, nrow(datos):2, type="S", 
     main="Fusion levels - Chord - Average", 
     ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(datos.ch.UPGMA$height, nrow(datos):2, nrow(datos):2, col="red", cex=0.8)

# Gráfico del ancho de silueta. ---------
asw <- numeric(nrow(datos))
#creamos un loop para los cálculos
for(k in 2:(nrow(datos)-1)){
    sil <- silhouette(cutree(datos.ch.UPGMA, k=k), datos.ch)
    asw[k] <- summary(sil)$avg.width}
k.best <- which.max(asw)

plot(1:nrow(datos), asw, type="h", 
     main="Silhouette-optimal number of clusters, Average", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an median silhouette width of", max(asw), "\n")
# Salto entre 2 y 3. Seleccionar 2.

# Estadístico de Mantel (Pearson) -----------------------
# Creamos una función para calcular la matriz de distancia binaria con los grupos
grpdist <- function(X){
    require(cluster)
    gr <- as.data.frame(as.factor(X))
    distgr <- daisy(gr, "gower")
    distgr}

# calculamos los estadísticos para el método de Average
kt <- data.frame(k=1:nrow(datos), r=0)
for(i in 2:(nrow(datos)-1)){
    gr <- cutree(datos.ch.UPGMA, i)
    distgr <- grpdist(gr)
    mt <- cor(datos.ch, distgr, method="pearson")
    kt[i,2] <- mt}
kt
k.best <- which.max(kt$r)

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Average", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# En este caso también nos dice que son 2 grupos los optimos.

# Todos los métodos están marcando como número óptimo de clusters = 2.

# Validación final.  --------------
#   Gráfico de silueta para la partición final 
#   observamos si los miembros de los grupos están bien clasificados.
k <- 2 # para 2 clusters
cutg <- cutree(datos.ch.UPGMA, k=k)
sil <- silhouette(cutg, datos.ch)
rownames(sil) <- row.names(datos)

plot(sil, main="Silhouette plot - Chord - Average", 
     cex.names=0.8, col=2:(k+1), nmax=100)

# PASO 4. Representación e interpretación de los resultados.   =============
# Clusters con colores por grupo
source("hcoplot.R")
hcoplot(datos.ch.UPGMA, datos.ch, k=2)
# Aunque realiza correctamente tanto el dendrograma como el gráfico de silueta
# para la partición final. El número de clusters respecto al original no es correcto.
# En el original son 3 clases, mientras que sólo ha detectado 2 por el método jerárquico.
# Está agrupando 2 tipologías dentro de un cluster.


### Utilización del método NO JERARQUICO

# Paso 1.

datos <- iris[,-5]

# Normalizamos los datos
datos.norm <- decostand(datos, "normalize") 
summary(datos.norm)

# PASO 2. Elección de la técnica de conglomeración. =============
# Si fijamos de antemano el número de grupos a formar:
library("factoextra")
fviz_nbclust(datos.norm, kmeans, method = "gap_stat")

set.seed(1) # para hacer repetible los resultados
(datos.kmeans <- kmeans(datos.norm, centers=3, nstart=100)) # le decimos el numero de grupos 3 grupos

# PASO 3. Elección del número de conglomerados finales y validación.  ========
# Si queremos determinar el número óptimo de grupos a formar en vez de seleccionarlos como en el paso anterior (3):  
# Elección del número de conglomerados finales
# Regla de Calinski y Harabasz (calinski). mayor valor = mejor k 
# Índice de estructura simple (SSI). mayor valor = mejor k
# Suma de errores cuadrados (SSE). mayor k = menor SSE. Mirar cuando no hay una mejora significativa en SSE.
# Validación: Gráfico del ancho de silueta.  

# calculamos particiones desde 2 a 6 clusters
datos.KM.cascade <- cascadeKM(datos.norm, inf.gr=2, sup.gr=6, iter=100, criterion="ssi")
datos.KM.cascade$results # vemos que el criterio SSI se maximiza para 3 grupos.
# Miramos el indice de estructura simple (SSI). mayor valor = mejor k, pues Calinski aún no lo hemos hallado (es 0.186 - 3 grupos)
head(datos.KM.cascade$partition) # vemos la partición creada y la graficamos
plot(datos.KM.cascade, sortg=TRUE) 

# Validación
(datos.kmeans <- kmeans(datos.norm, centers=3, nstart=100)) # con 3 grupos
dissE <- daisy(datos.norm) # para obtener la matriz de disimilaridades
sk <- silhouette(datos.kmeans$cl, dissE) # gráfico de silueta
plot(sk)
# Vemos que están bien clasificados, todos desde el 0 hacia la derecha.
# Pero por la tabla sabemos que 5 elementos del primer cluster no están correctamente clasificados. Pues hay 50 de cada modelo.

# PASO 4. Representación e interpretación de los resultados.  ========

# Graficamos los clusters
clusplot(datos.norm, datos.kmeans$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
# En dicho gráfico apreciamos que existe una intersección entre el cluster 1 y 2, debido a los elementos que están catalogados en
# el grupo 1 pero deberían realmente estar dentro del grupo 2.

# Se podría crear un dendrograma para ver cómo deberían catalogarse.
library("factoextra")
res <- hcut(datos.norm, k = 3, stand = TRUE)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))

# Respecto a los datos inciales, el método No Jerarquico ha trabajado mejor que el Jerárquico pues el ultimo metodo 
# no detectó correctamente los clusters, pasando de 3 que posee el original a 2 que detecto dicho metodo.

# ============================================================================================================================================

# Ejercicio 2.

# ============================================================================================================================================

# ANALISIS DISCRIMINANTE (LDA).
# Con el conjunto de datos "iris" realice el analisis discriminante lineal (LDA). Grafique e interprete los resultados.

# Cargamos los datos

data(iris)
head(iris)
summary(iris)
str(iris)

# cargamos los paquetes que vamos a necesitar
library(ade4)
library(vegan)
library(ellipse)
library(mvnormtest)
library(psych)
library(MASS)
library(klaR)

# Realizamos la agrupacion previa.

datos.iris <- iris[,-5]
datos.iris.hel <- decostand(datos.iris, "hellinger")
gr <- cutree(hclust(vegdist(datos.iris.hel, "euc"), "ward.D"),3)
table(gr)

# Evaluamos los supuestos de:

# 1) Homogeneidad. ------------------------

datos.iris.d1 <- dist(datos.iris)
datos.MHV <- betadisper(datos.iris.d1, gr)
anova(datos.MHV)
permutest(datos.MHV)
# Vemos que NO son homogéneas pues tenemos en el test de permutacion un
# p-valor menor al 0.05 (0.001 exactamente), por lo que tenemos que transformar.
# miramos la o las variables a transformar...

par(mfrow=c(1,ncol(datos.iris)))
for(i in 1:ncol(datos.iris)){ 
    boxplot(datos.iris[,i],main=colnames(datos.iris)[i])
}

# transformación usarmos sqrt
datos.iris.pars2 <- sqrt(datos.iris)

datos.iris.pars2.d1 <- dist(datos.iris.pars2)
datos.MHV <- betadisper(datos.iris.pars2.d1, gr)
permutest(datos.MHV)
# Ahora ya no existen problemas de homogeneidad (p-valor más alto del 0.05)

# 2) Normalidad multivariante. ----------------  

par(mfrow=c(1,ncol(datos.iris.pars2)))
for(j in 1:ncol(datos.iris.pars2)){
    hist(datos.iris.pars2[,j])}

#con el paquete mvnormtest
mshapiro.test(t(datos.iris.pars2))
# no rechazamos la normalidad multivariante (p-value mayor de 0.05)

# Multicolinealidad
as.dist(cor(datos.iris.pars2))
# Son ligeramente altas, alguna en el borde del límite con 0.97.

# Linealidad
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(datos.iris.pars2, lower.panel = panel.smooth, upper.panel = panel.cor)
# Siendo el más bajo 0.12 correspondiente a Sepal.Length con Sepal.Width, por lo general la relacion entre las
# variables puede asumirse bastante lineal.

# Exploración gráfica
par(mfrow=c(2,3))
partimat(factor(gr) ~ ., data = datos.iris.pars2 , method = "lda", nplots.vert=1)
# Por variables vemos cómo se dividiría el espacio entre los objetos (3 grupos que hemos seleccionado).
# Se observa que las particiones son aceptables a excepción de las creadas en las variables Sepal.Lenght con
# Sepal.Width, que, de hecho, arroja un ratio de error superior al resto con 0.19.

# Realizamos el LDA. =======================
datos.iris.pars2.df<- as.data.frame(datos.iris.pars2)
(datos.lda <- lda(gr ~., data=datos.iris.pars2.df))
# Coeficientes de las funciones discriminantes,
# LD1 = 4.1*Sepal.Length + 5.63*Sepal.Width - 9.62*Petal.Length - 5.78*Petal.Width
# LD2 = 1.82*Sepal.Length + 7.97*Sepal.Width - 2.26*Petal.Length + 3.57*Petal.Width

# Gráficos. ----------------------
#scores. posición de los objetos en el espacio de covariables
(Fp <- predict(datos.lda)$x)
# clasificación de los objetos
(datos.class <- predict(datos.lda)$class)

par(mfrow=c(1,1))
plot(Fp[, 1], Fp[, 2], type="n")
text(Fp[, 1], Fp[, 2], col=c(as.numeric(datos.class)+1))
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
# elipses del 95% alrededor de los grupos
for(i in 1:length(levels(as.factor(gr)))){
    cov <- cov(Fp[gr==i, ])
    centre <- apply(Fp[gr==i, ], 2, mean)
    lines(ellipse(cov, centre=centre, level=0.95))
}
# Hay solapamiento entre grupos azul y verde, mientras que el rojo está muy delimitado.


# Evaluamos la magnitud de la discriminación canónica. -------------------------

# 1) Valores propios.   
datos.lda$svd^2 #vemos que la primera función es la que tiene mayor poder discriminatorio. El primer eje.
100*datos.lda$svd^2/sum(datos.lda$svd^2) #en términos relativos ("proportion of trace"). El primer eje el 99% de lo que logran discriminar los ejes.

# 2) Correlaciones canónicas.   
#scores: valores de cada observación en cada eje canónico.
(scores <- predict(datos.lda)$x)
#calculamos la correlación canónica para cada eje
summary(lm(scores~gr))
# Vemos que el modelo es significativo.
# El R-cuadrado múltiple corresponde al valor del coeficiente de correlación canónica (cuadrado)
# R2 = 1.2% no discrimina bien entre grupos.

# Para evaluar la precisión de la clasificación, usaremos: --------------
# Vectores propios normalizados: coeficientes (estandarizados) de las funciones discriminantes.
(Cs <- datos.lda$scaling)
# Scores: posición de los objetos en el espacio de covariables
(Fp <- predict(datos.lda)$x)

# Clasificación de los objetos
(datos.class <- predict(datos.lda)$class)

# probabilidades posteriores de los objetos que pertenecen a los grupos
(datos.post <- predict(datos.lda)$posterior)

# Calculamos la tabla de contingencia de las clasificaciones a prori y las predichas
(datos.table <- table(gr, datos.class))
# las filas corresponden a los grupos (los grupos conocidos a priori) y 
# las columnas corresponden a los grupos predichos por el modelo LDA.
# La diagonal da las predicciones correctas (frecuencia) y lo demás son errores

# proporción de clasificaciones correctas para cada grupo
diag(prop.table(datos.table, 1)) # Sería 100%, 96% y 100% en el tercero.

# proporción de calsificaciones correctas (CCR) de manera global
sum(diag(datos.table))/sum(datos.table)
# CCR=98.6% -> proporción de clasificaciones correctas de manera global.
# El resultado es que NO discrimina del todo bien pero aún así las clasificaciones son un 98.6% correctas.
# Hay 4 tipos que no logra agrupar correctamente. El grupo 1 y 3 sí los detecta correctamente versus 
# el resto.

# Vamos a corregir los problemas del azar. ------------
# Cierto porcentaje de los datos estarán clasificados correctamente solo por azar. 
# Proporcional al tamaño del grupo. 
# Vamos a quitar dicha proporción para ver el poder real discriminativo de nuestro análisis.

# Coeficiente Kappa de Cohen.  
cohen.kappa(datos.table)
# el 98% del 98% (CCR) corresponden a clasificaciones correctas que no se deben al azar.

# Validación del MODELO.  -------------------
# Jackknife.
(datos.lda.jac <- lda(gr ~ ., data=datos.iris.pars2.df, CV=TRUE))

# proporción de clasificaciones correctas
datos.lda.jac <- datos.lda.jac$class
datos.lda.table <- table(gr, datos.lda.jac)
diag(prop.table(datos.lda.table, 1))
# Vemos que nos arroja los mismos datos que la anterior.

# Interpretación de las funciones canónicas. Usaremos:  --------------------

# vectores propios (o coeficienes canónicos): corresponden a los 
# datos originales (no estandarizadas, no normalizadas) 
# pero no son coeficientes del todo interpretables.  
datos.lda$scaling

# por lo que se suele usar la estructura canónica de correlaciones (las correlaciones de Pearson 
# entre cada variable discriminante y cada función canónica -scores-) 
# indica la fuerza y naturaleza (positiva o negativa) de la relación y 
# puede utilizarse para interpretar mejor o más adecuadamente cada eje.
# cargas al cuadrado indican el porcentaje de varianza de la variable 
# que está dada por esa función canónica.
# Podemos usar las originales o las estrucutas de correlaciones canónicas.
dim<-ncol(scores)
cutoff<-0

#calculamos la estructura de correlationes
z<-cor(datos.iris.pars2,scores[,1:dim])
z<-round(z,digits=2)
z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
(z<-as.data.frame(z))

# Menos la variable Sepal.Width, el resto aportan de manera negativa y alta al primer eje, 
# Todas aportan de manera positiva al segundo eje.
# En general vemos si aporta de forma positiva o negativa y por el valor, cual es más importante.
# juegan un papel importante para la discriminación a lo largo del primer eje.
# En esta caso Petal.Length y Petal.Width son las que juegan un papel más importante en el eje 1
# Sepal.Width es la que lo juega en el eje 2.

# ============================================================================================================================================

# Ejercicio 3.

# ============================================================================================================================================

library(vegan)
data(iris)
head(iris)
attach(iris)
datos_iris <- iris[,-5]

### MANOVA. Version Parametrica.

(m1 <- manova(as.matrix(datos_iris)~iris$Species))
# muestra residuales y otros valores, pero para buscar si existen diferencias significativas, realizamos lo siguiente:

summary(m1)
# EL metodo de pillai nos dice que hay diferecias significativas p es menor de 0.05

summary(m1, test = "Wilks")
# EL metodo de wilks nos dice también que hay diferecias significativas p es menor de 0.05
# POR LO QUE RECHAZAMOS LA IGUALDAD ENTRE GRUPOS.

# salida de Anova univariado. Muestra los valores para cada de las variables:
head(summary.aov(m1))

# comparaciones múltiples post hoc:
table(iris$Species)
# Recordamos que hay 3 tipos. Realizamos las pruebas 2 a 2...

# Vamos a corregir el error: se recomienda α = 0.15 y corregir por α/m donde m=nºtests realizados. 0.15/3=0.05
# si es menor a 0.025 cosideramos una diferencia significativa, en caso contrario no.
summary(manova(as.matrix(datos_iris) ~ Species, 
               subset = Species %in% c("setosa", "versicolor")))

summary(manova(as.matrix(datos_iris) ~ Species, 
               subset = Species %in% c("setosa", "virginica")))

summary(manova(as.matrix(datos_iris) ~ Species, 
               subset = Species %in% c("versicolor", "virginica")))

# En los tres casos los p-valores son menores a 0.05 por lo que son significativos y existen diferencias entre los 3 conjuntos o especies.

# Probamos otra opción:
library(RVAideMemoire)
pairwise.perm.manova(datos_iris,Species,nperm=49,test="Hotelling-Lawley")
# Nos muestra que sí hay diferencia entre los 3 (p-valor es 0.02 menor a 0.05)

# Probamos otra opción:
library(ICSNP)
HotellingsT2(as.matrix(datos_iris[Species %in% c("setosa", "versicolor"),]) ~ Species[Species %in% c("setosa", "versicolor")]) #p-value < 0.05 igual que el resultado de arriba.
HotellingsT2(as.matrix(datos_iris[Species %in% c("virginica", "versicolor"),]) ~ Species[Species %in% c("virginica", "versicolor")]) #p-value < 0.05 igual que el resultado de arriba.
HotellingsT2(as.matrix(datos_iris[Species %in% c("setosa", "virginica"),]) ~ Species[Species %in% c("setosa", "virginica")]) #p-value < 0.05 igual que el resultado de arriba.

# Como es una prueba paramétrica, vamos a probar los supuestos de normalidad y homogeneidad de varianza:
library(car)
par(mfcol=c(2,2))
hist(residuals(m1));
qqnorm(residuals(m1));qqline(residuals(m1)) # outliers
plot(m1$residuals~m1$fitted.values)

# Los residuales si parecen bastante normales (superior izda.), en el gráfico inferior se aprecia como van separándose de la línea ligeramente por los extremos,
# Se aprecia un patron de embudo en el gráfico superior derecho, lo que nos dice que existe un problema de homogeneidad de varianza.

#normalidad
library(MVN)
mardiaTest(residuals(m1)) # Data are not multivariate normal. No son bastante normales.
library(mvnormtest)
mshapiro.test(t(datos_iris)) # En este caso nos quedamos con este dato y los gráficos. p-value = 0.02 que nos dice que son normales.

#homogeneidad
library(biotools)
boxM(datos_iris, Species) #p-value < 2.2e-16. No podemos considerar que haya homogeneidad.


### Como no cumple con los requisitos para usar con el manova parametrico........     

# -------------------------------------------------------- #
# ADONIS, PERMANOVA o NPMANOVA.  
# -------------------------------------------------------- #

adonis(datos_iris ~ Species, data=iris, permutations=99) # existen diferencias significativas en cuanto a las especies.

# comparaciones múltiples post hoc. con dist utiliza adonis
pairwise.perm.manova(dist(datos_iris,"euclidean"),Species,nperm=49) 
# Vemos diferencias significativas entre las 3, las detecta como en el anterior metodo.

# También se puede usar...

# -------------------------------------------------------- #
### MRPP
# -------------------------------------------------------- #

(iris.mrpp <- mrpp(datos_iris, Species, distance = "bray")) #default=distance = "euclidean",
# como significance = 0.001, concluimos que las especies difieren significativamente en funcion de las 4 variables.
# A: 0.6

# graficamos
def.par <- par(no.readonly = TRUE)
layout(matrix(1:2,nr=1))
plot(iris.ord <- metaMDS(datos_iris), type="text", display="sites" )
ordihull(iris.ord, Species)

with(iris.mrpp, {
      fig.dist <- hist(boot.deltas, xlim=range(c(delta,boot.deltas)), 
                       main="Test of Differences Among Groups")
      abline(v=delta); 
      text(delta, 2*mean(fig.dist$counts), adj = -1.5,
           expression(bold(delta)), cex=1)  }
)
par(def.par)
# Se puede apreciar que los tres grupos estan diferenciados, pero dos de ellos se solapan en algunas observaciones.

## meandist
iris.md <- meandist(vegdist(datos_iris), Species)
iris.md
summary(iris.md)
plot(iris.md)
# Setosa es el que más difiere, aunque los tres grupos son diferentes entre sí. Entre Versicolor y Virginica no hay tan clara
# diferencia debido a las observaciones que quedan solapadas entre ambos grupos.

# También se puede usar...

# -------------------------------------------------------- #
# ANOSIM. 
# -------------------------------------------------------- #
# utiliza distancias
# utiliza pruebas de permutación Monte Carlo

iris.dist <- vegdist(datos_iris) #default method="bray"
iris.ano <- anosim(iris.dist, Species)
summary(iris.ano)
#R= 0.85 (muy alto), con Significance = 0.001 hay diferencias significativas, concluimos que los grupos difieren significativamente en funcion de las 4 variables.

#gráfico
plot(iris.ano)
# la diferencia entre los grupos es bastante grande (el primero de ellos)
# dentro de cada grupo la variacion es similiar entre los 3.