# Iniciamos con el método Jerárquico o HC
# Cargamos los datos y paquetes
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)

attach(iris)


# Paso 1 Selección de la muestra de datos.
iris_species <- iris[,5]
iris_data <- iris[,1:4]

iris_norm <- decostand(iris_data, "normalize")

# Paso 2. Calculamos la matriz de distancias
iris.ch <- vegdist(iris_norm, "euc") # cálculo de matriz de distancias euclideas

# Trabajamos varios métodos para posteriomente seleccionar el mejor.
par(mfrow=c(2,3)) #para disponer juntos los gráficos
iris.ch.single <- hclust(iris.ch, method="single") # Método "single" (método del vecino más cercano o mínima distancia).  
plot(iris.ch.single) # graficamos el dendrograma
iris.ch.complete <- hclust(iris.ch, method="complete") # Método "complete" (método del vecino más lejano o máxima distancia).
plot(iris.ch.complete) 
iris.ch.UPGMA <- hclust(iris.ch, method="average") # Método "UPGMA" o "Average" (método del promedio)
plot(iris.ch.UPGMA)
iris.ch.centroid <- hclust(iris.ch, method="centroid") # Método "centroid" (método del centroide).  
plot(iris.ch.centroid)
iris.ch.ward <- hclust(iris.ch, method="ward.D") # Método de mínima varianza de Ward
plot(iris.ch.ward)
iris.ch.ward$height <- sqrt(iris.ch.ward$height) # con transformación de la raíz cuadrada 
plot(iris.ch.ward)

# Correlaciones cofenéticas -------  
iris.ch.single.coph <- cophenetic(iris.ch.single) # Single
cor(iris.ch, iris.ch.single.coph)
iris.ch.comp.coph <- cophenetic(iris.ch.complete) # Complete
cor(iris.ch, iris.ch.comp.coph)
iris.ch.UPGMA.coph <- cophenetic(iris.ch.UPGMA) # Average
cor(iris.ch, iris.ch.UPGMA.coph)
iris.ch.ward.coph <- cophenetic(iris.ch.ward) # Ward 
cor(iris.ch, iris.ch.ward.coph)
cor(iris.ch, iris.ch.ward.coph, method="spearman")
# LA MAYOR CORRESPONDE AL MÉTODO DE LA MEDIA.

# Graficas las correlaciones:
par(mfrow=c(2,2))
plot(iris.ch, iris.ch.single.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Single linkage", paste("Cophenetic correlation =",
                                    round(cor(iris.ch, iris.ch.single.coph),3))))
abline(0,1);  lines(lowess(iris.ch, iris.ch.single.coph), col="red")

plot(iris.ch, iris.ch.comp.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Complete linkage", paste("Cophenetic correlation =",
                                      round(cor(iris.ch, iris.ch.comp.coph),3))))
abline(0,1);  lines(lowess(iris.ch, iris.ch.comp.coph), col="red")

plot(iris.ch, iris.ch.UPGMA.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("UPGMA", paste("Cophenetic correlation =",
                           round(cor(iris.ch, iris.ch.UPGMA.coph),3))))
abline(0,1);  lines(lowess(iris.ch, iris.ch.UPGMA.coph), col="red")

plot(iris.ch, iris.ch.ward.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), 
     ylim=c(0,max(iris.ch.ward$height)),
     main=c("Ward clustering", paste("Cophenetic correlation =",
                                     round(cor(iris.ch, iris.ch.ward.coph),3))))
abline(0,1);  lines(lowess(iris.ch, iris.ch.ward.coph), col="red")
# TAMBIEN SE APRECIA QUE EL QUE MÁS CENTRADO ESTÁ ES EL DE LA MEDIA. COINCIDE CON EL ANTERIOR.

# Distancia de Gower (1983) ------
(gow.dist.single <- sum((iris.ch-iris.ch.single.coph)^2))
(gow.dist.comp <- sum((iris.ch-iris.ch.comp.coph)^2))
(gow.dist.UPGMA <- sum((iris.ch-iris.ch.UPGMA.coph)^2))
(gow.dist.ward <- sum((iris.ch-iris.ch.ward.coph)^2))
# TAMBIÉN COINCIDE CON EL ANTERIOR. TODO PARECE CORRECTO.


# Paso 3. Elección del número de conglomerados.
# Elegimos el de la Media que es el que mejor se comporta.

# Gráfico de nivel de fusión ------------
plot(iris.ch.ward$height, nrow(iris_data):2, type="S", 
     main="Fusion levels - Chord - Ward", 
     ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(iris.ch.UPGMA$height, nrow(iris_data):2, nrow(iris_data):2, col="red", cex=0.8)
# Salto entre 2 y 3. Seleccionar 3.

# Gráfico del ancho de silueta. ---------
asw <- numeric(nrow(iris_data))
#creamos un loop para los cálculos
for(k in 2:(nrow(iris_data)-1)){
    sil <- silhouette(cutree(iris.ch.UPGMA, k=k), iris.ch)
    asw[k] <- summary(sil)$avg.width}
k.best <- which.max(asw)

plot(1:nrow(iris_data), asw, type="h", 
     main="Silhouette-optimal number of clusters, Ward", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# Salto entre 2 y 3. Seleccionar 2.

# Estadístico de Mantel (Pearson) -----------------------
# Creamos una función para calcular la matriz de distancia binaria con los grupos
grpdist <- function(X){
    require(cluster)
    gr <- as.data.frame(as.factor(X))
    distgr <- daisy(gr, "gower")
    distgr}

# calculamos los estadísticos para el método de Ward
kt <- data.frame(k=1:nrow(iris_data), r=0)
for(i in 2:(nrow(iris_data)-1)){
    gr <- cutree(iris.ch.ward, i)
    distgr <- grpdist(gr)
    mt <- cor(iris.ch, distgr, method="pearson")
    kt[i,2] <- mt}
kt
k.best <- which.max(kt$r)

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# En este caso nos dice que son 2 grupos los optimos.

