# Ejercicio 1

library(vcd)
attach(Arthritis)
head(Arthritis)
library(ggplot2)

# Comparar edades de pacientes sin mejoría con las de mejoría
mejoria <- subset(Arthritis, Improved == "Marked" | Improved == "None")
mejoria$Improved <- droplevels(mejoria$Improved)
detach(Arthritis)

shapiro.test(mejoria$Age)
ggplot(mejoria, aes(factor(Improved), Age)) + geom_boxplot() + labs(title = "Tratamiento")
# No sigue una distribución normal.
# Existen Outliers.

# pruebas paramétricas (suponemos Normalidad y homocedasticidad)
# Hipotesis Nula <- Son iguales.
# Hipotesis alternativa <- son diferentes.

t.test(Age ~ Improved, data = mejoria)

# p-valor < 0.05. Tenemos evidencia en los datos para afirmar con un 95% de confianza
# que las edades no son iguales en ambos grupos. La edad media en el 
# grupo de mejora evidente es superior (la media) que en el de ninguna mejora.

# pruebas no paramétricas.
# Hipotesis Nula <- Son iguales.
# Hipotesis alternativa <- son diferentes.

wilcox.test(Age ~ Improved, data = mejoria)

# p-valor < 0.05. Tenemos evidencia en los datos para afirmar con un 95% de confianza
# que las edades no son iguales en ambos grupos. La edad media en el 
# grupo de mejora evidente es superior (la media) que en el de ninguna mejora.

# pruebas Robustas
# Hipotesis Nula <- Son iguales.
# Hipotesis alternativa <- son diferentes.

library(DescTools)
YuenTTest(formula = Age ~ Improved, data = mejoria)

# p-valor < 0.05. Tenemos evidencia en los datos para afirmar con un 95% de confianza
# que las edades no son iguales en ambos grupos. La edad media en el 
# grupo de mejora evidente es superior (la media) que en el de ninguna mejora.



# Ejercicio 2
library(MASS)
attach(immer)
head(immer)

# En esta caso son muestras relacionadas.
cosecha <- immer[,c("Y1", "Y2")]
boxplot(cosecha)
# Se aprecia la existencia de outliers.

# hacemos la preba de normalidad
shapiro.test(cosecha$Y1)
# No sigue una distribución normal.
shapiro.test(cosecha$Y2)
# Sigue una distribución normal

# PRUEBAS PARAMÉTRICAS (suponemos Normalidad y homocedasticidad). Hay relación.
# Hipotesis Nula <- Poseen el mismo peso medio. En cada campo por año.
# Hipotesis alternativa <- No poseen el mismo peso medio (Y1 fue mayor que Y2).

(t1 <- t.test(cosecha$Y1, cosecha$Y2, paired = T, alternative = "greater"))
# p-valor < 0.05: Se acepta la hipótesis alternativa. El año Y2 obtuvo un peso menor a Y1

# PRUEBAS NO PARAMÉTRICAS
# Hipotesis Nula <- Poseen el mismo peso medio. En cada campo por año.
# Hipotesis alternativa <- No poseen el mismo peso medio (Y1 fue mayor que Y2).

(t1 <- wilcox.test(cosecha$Y1, cosecha$Y2, paired = T, alternative = "greater"))
# p-valor < 0.05: Se acepta la hipótesis alternativa. El año Y2 obtuvo un peso menor a Y1


# PRUEBAS ROBUSTAS. Existen outliers (en este caso).
# Hipotesis Nula <- Poseen el mismo peso medio.
# Hipotesis alternativa <- No poseen el mismo peso medio (Y1 fue mayor que Y2).

library(WRS2)

(tC <- yuend(x = cosecha$Y1, y = cosecha$Y2))
# p-valor < 0.05: Se acepta la hipótesis alternativa. El año Y2 obtuvo un peso menor a Y1
detach(immer)


# Ejercicio 3

attach(PlantGrowth)
head(PlantGrowth)
levels(PlantGrowth$group)
table(PlantGrowth$group)

ggplot(PlantGrowth, aes(factor(group), weight)) + geom_boxplot() + labs(title = "Peso seco Plantas")
# Comentar que son variables independientes, n < 30 (muestra pequeña)
# Existen outliers y se puede apreciar también en el gráfico que hay diferencia en los pesos entre los 3 grupos.

# PRUEBAS PARAMÉTRICAS
# Hipotesis Nula <- Todas las medias son iguales (pesos similares, pesos medios similares)
# Hipotesis alternativa <- Existe alguna media diferente.
m1 <- aov(weight ~ group, data = PlantGrowth)
summary(m1)
# Existen diferencias entre los pesos de los lotes, aunque no muy significativas. p-valor = 0.02
# Marca además una estrella.
# La moyor diferencia aparece entre tratamiento 2 y tratamiento 1 que fue el que aportó peor resultado en cuanto a peso (el de menor peso).

## Pruebas post hoc
plot(TukeyHSD(m1), las=2)
TukeyHSD(m1)
# La moyor diferencia aparece entre tratamiento 2 y tratamiento 1 que fue el que aportó peor resultado en cuanto a peso (el de menor peso).


# PRUEBAS NO PARAMÉTRICAS.
# Hipotesis Nula <- Todas las medias son iguales (pesos similares, pesos medios similares)
# Hipotesis alternativa <- Existe alguna media diferente.

kruskal.test(weight ~ group, data = PlantGrowth)
# Existen diferencias entre los pesos de los lotes, aunque no muy significativas. p-valor = 0.02

## Pruebas post hoc
require(PMCMR)
posthoc.kruskal.nemenyi.test(x = weight, g = group, data = PlantGrowth, dist = "Chisquare")
# De nuevo muestra que entre trt1 y trt2 existe mayor diferencia.


# PRUEBAS ROBUSTAS.
# Hipotesis Nula <- Todas las medias son iguales (pesos similares, pesos medios similares)
# Hipotesis alternativa <- Existe alguna media diferente.

boxplot(weight ~ group, data = PlantGrowth, las=2)

library(WRS2)
t1way(weight ~ group, data = PlantGrowth)
# Rechazamos H0. Existen diferencias entre los pesos de los lotes. p-valor = 0.008

## Pruebas post hoc
lincon(weight ~ group, data = PlantGrowth)
# La mayor diferencia se detecta entre trt1 y trt2.
detach(PlantGrowth)


# Ejercicio 4
library(car)
library(reshape)
attach(WeightLoss)
head(WeightLoss)
# Separamos peso y estado de animo y reindexamos
data <- subset(WeightLoss, group != "Control")
# Vamos a detectar la existencia de outliers en las diferentes columnas.
library(rapportools)
rp.outlier(data$wl1)
rp.outlier(data$wl2)
rp.outlier(data$wl3)
rp.outlier(data$se1)
rp.outlier(data$se2)
rp.outlier(data$se3)
# Existencia en outlier valor 11 en la columna de datos se3.

dataW <- cbind(ID = as.factor(1:22),data[,2:4])
rownames(dataW) <- NULL
dataS <- cbind(ID = as.factor(1:22),data[,5:7])
rownames(dataS)<- NULL
# pasamos las columnas de wl como una columna variable.
dataLS<-melt( dataS, id.vars = "ID", measure.vars = 2:4)
dataLS <- dataLS[order(dataLS$ID),]
dataLW<-melt( dataW, id.vars = "ID", measure.vars = 2:4)
dataLW<- dataLW[order(dataLW$ID),]
row.names(dataLS)<-NULL
row.names(dataLW)<-NULL

# PRUEBAS PARAMETRICAS
# Hipotesis Nula <- Todas las medias son iguales (no existe aumento)
# Hipotesis alternativa <- Existe alguna media diferente.

library(ez)
library(WRS2)
ezANOVA(data=dataLS, dv=value, wid=ID, within = variable)
# Primera línea p-valor < 0.05. Rechazamos igualdad de medias. Hay diferencia.
# La prueba de esfericidad dice que también deberíamos rechazarla (Hip Nula).

# Vamos a ver entre que medidas hay diferencia.
aggregate(value~variable, dataLS, mean)
with(dataLS, pairwise.t.test(value, variable, p.adjust.method = "bonferroni", paired = T))
# Hay diferencias entre los 3.
library(gplots)
attach(dataLS)
plotmeans(value ~ variable, xlab="Medicion", ylab = "valor")
# Los valores más significativos en cuanto aumento se dieron en SE3". SE3 posee los niveles de autoestima más altos.


# PRUEBAS NO PARAMETRICAS
# Hipotesis Nula <- Todas las medias son iguales (no existe aumento)
# Hipotesis alternativa <- Existe alguna media diferente.

# representamos los datos
x <- reshape(dataLS, v.names = "value", idvar="ID", timevar = "variable", direction = "wide")
boxplot(x[,-1])
# se aprecia que en los 3 casos hay diferencias, más significativas en SE3 con el resto.

friedman.test(as.matrix(x[,-1]))
# Viendo el p-valor se concluye que la prueba es significativa. Hay diferencias, se rechaza la hipotesis nula.
# Vemos las diferencias:

pairwise.wilcox.test(dataLS$value, dataLS$variable, p.adjust.method = "bonferroni", exact=F, paired = T)
# SE3 es el que mayor aumento obtuvo. Mayor diferencia entre SE3 y SE2 y luego SE3 con SE1. SE3 posee los niveles de autoestima más altos.


# PRUEBAS ROBUSTAS
# Hipotesis Nula <- Todas las medias son iguales (no existe aumento)
# Hipotesis alternativa <- Existe alguna media diferente.

# Mostramos de nuevo los valores
boxplot(value~variable, data= dataLS)

rmanova(y = value, groups = variable, blocks = ID)
# Se obtiene un p-valor menor a 0.05 por lo que rechazamos la hipotesis nula, hay diferencias entre ellos.

rmmcp(y = value, groups = variable, blocks = ID)
# En todos los casos nos dan valores significativos.Entre SE2 y SE3 es donde obtenemos valores más diferenciados. SE3 posee los niveles de autoestima más altos.
rm(list = ls())


# EJERCICIO 5

Convictions <-matrix(c(2, 10, 15, 3), nrow = 2, 
                     dimnames = list(c("Dizygotic", "Monozygotic"),c("Convicted", "Not convicted")))
Convictions

# Siendo p1 monocigóticos y p2 dicigóticos:
# Hipotesis Nula <- p1 menor o igual a p2
# Hipotesis Alternativa <- p1 mayor que p2

# Como la muestra es pequeña, utilizamos Fisher.

prop.test(Convictions)
# Como p-valor es menor a 0.05 rechazamos la hipotesis nula, por lo que la proporción de monocigóticos convictos es mayor que la de dicigóticos.
rm(list = ls())


# Ejercicio 6

data("infert")
head(infert)
dim(infert)

# aislamos las columnas con las que trabajar.
datos <- table(infert[c("education", "induced")])
datos

# proporciones para cada una de las variables de educación, mostrando los abortos inducidos.
plot(datos)

# H0 <- Las variables Education e Induced son independientes.
# H1 <- Ambas variables son dependientes.

chisq.test(datos)
# Como p-valor es menor a 0.05 rechazamos la hipotesis nula, por lo que la sí existen dependencia entre las variables de Educación y Abortos inducidos.

# Miramos que grado de asociación existe mediante Medidas de asociación global.
library(vcd)
assocstats(datos)
# Muestra una asocación débil (Se aprecia en la V de Cramer
rm(list = ls())