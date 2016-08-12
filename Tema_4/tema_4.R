library(quantreg)
library(car)
data("engel")
head(engel)
attach(engel)

boxplot(engel)
# Se aprecia que no siguen normalidad, pero existen outliers en ambas variables.

shapiro.test(engel$income)
# p < 0.05. Se rechaza la hipótesis nula de existencia de normalidad.
shapiro.test(engel$foodexp)
# p < 0.05. Se rechaza la hipótesis nula de existencia de normalidad.



# Hallar la correlación entre variables. Se usará correlaciones robustas al existir outliers.

plot(income, foodexp, pch=20)
abline(coef(lm(foodexp~income)), col='green')
# En principio sí parece haber linealidad.

library(WRS)

pbcor(income, foodexp)
wincor(income, foodexp)
# La correlación entre ambas variables es muy fuerte y positiva.
# p-valor menor a 0.05 por lo que rechazamos la hipótesis nula (correlación es 0). Existe una relación significativa entre ambas variables.



# Ajuste del modelo regresión lineal simple clásico y robusta.
# Alimentos (VD) en función de la renta (VI). Es decir, foodexp en función de income.

summary(fitLS <- lm(foodexp~income))
# Como el intercepto (beta 0) e income(beta1), poseen un p-valor < 0.05 rechazamos la hipótesis nula por lo que ambas son significativos, deben aparecer en la ecuación.
# La ecuación quedaría como: foodexp = 0.48518 * income + 147.47539
# F-statistic (ANOVA), el p-value es < 0.05 nos aporta que el modelo es significativo.
# Los parámetros del modelo intercepto y pendiente son ambos significativos (p-valor < 0.001)

# Ajuste global del Modelo y Bondad de Ajuste...
# Se nos muestra en F-statistic del paso anterior. Seria: F(1,233)=1141, p-value < 0.001
# La bondad de ajuste viene dada por R2 que en este caso es de 0.8304, el modelo es muy eficaz (entre 0 y 1).

# graficamos...
plot(foodexp~income, pch=20)
abline(fitLS)

# Evalúe gráficamente los supuestos del modelo (diagnóstico) y la presencia de datos atípicos (outliers)

# Gráfico genérico del modelo:
par(mfrow=c(2,2))
plot(fitLS)

# La linealidad se puede apreciar en el gráfico de los datos anterior. Existe linealidad.

# Independencia de residuos
library(lmtest)
dwtest(fitLS)
# Se usó prueba de Durbin-Test. Se obtuvo: Autocorrelación positiva dwtest < 2.

# Vemos los gráficos anteriores uno a uno...

# Residuals vs. Fitted.
par(mfrow=c(1,1))
plot(fitLS, which=1)
# No aparece ningún tipo de patrón. Parece media 0. No hay heterocesdasticidad ni problemas de linealidad. No hay correlacionalidad en los residuos.
# Se aprecian 3 outliers 59, 105 y 138.

# Normal Q-Q
plot(fitLS, which=2)
# A excepción de los outliers comentados anteriormente, parece que se ajustan los datos a la recta.
# La hipótesis de que los errores siguen una distribución normal parece cumplirse.

# Scale - Location
plot(fitLS, which=3)
# Como en el caso del gráfico primero, no aparece ningún tipo de patrón.

# Residuals vs. Leverage
plot(fitLS, which=4)
# Muestra los 3 outliers que se deben estudiar.

# graficamos distancias de cook
cooks <- cooks.distance(fitLS)
hat <- lm.influence(fitLS)$hat

par(mfrow = c(2,1))
plot(income,foodexp, col=ifelse(cooks > quantile(cooks, .90), 'red', 'black'), pch=20)
plot(income,foodexp, col=ifelse(hat > quantile(hat, .90), 'green', 'black'), pch=20)

par(mfrow = c(2,1))
plot(income,predict(fitLS), col=ifelse(cooks > quantile(cooks, .90), 'red', 'black'), pch=20)
plot(income,predict(fitLS), col=ifelse(hat > quantile(hat, .90), 'green', 'black'), pch=20)
# Parece que hay outliers que influyen en la recta. Se debería trabajar mediante el uso de técnicas robustas (Ya comentado al inicio).

# Uso de validación cruzada para ver el ajuste del modelo.
library(caret)
datos <- engel
train(foodexp~income, data = datos, method = "lm")
# R2 es 0.84. Posee un buen ajuste.



# -------------------------------------------------------------------------------
# AJUSTE USANDO ltsReg PARA ROBUSTOS (Regresión por mínimos cuadrados recortados)
#--------------------------------------------------------------------------------

# Graficamos los outliers
library(car)
plot(income, foodexp, pch=20, main = "Regresión Paramétrica")
abline(coef(lm(foodexp~income)), col='green')
out <- influence.measures(fitLS)
showLabels(engel$foodexp, engel$income, row.names(engel), 
           id.method=row.names(engel)[row.names(engel) %in% names(which(apply(out$is.inf, 1, any)))])

# Existencia de outliers. Lo conveniente es aplicar regresión robusta.

# library(rrcov)
fitLSR <- ltsReg(income, foodexp)
summary(fitLSR)
# Como el income(beta1), poseen un p-valor < 0.05 rechazamos la hipótesis nula por lo que solo él debe aparecer en la ecuación. El intercepto no pues aceptamos la hipótesis nula.
# La ecuación quedaría como: foodexp = 0.68244 * income
# F-statistic (ANOVA), el p-value es < 0.001 nos aporta que el modelo es significativo.

# Ajuste global del Modelo y Bondad de Ajuste...
# Se nos muestra en F-statistic del paso anterior. Seria: F(1,209)=2930, p-value < 0.001
# La bondad de ajuste viene dada por R2 que en este caso es de 0.9334, el modelo es muy eficaz (entre 0 y 1).

# graficamos...
par(mfrow=c(1,2))
plot(income, foodexp, pch=20, main = "Regresión Paramétrica")
abline(coef(lm(foodexp~income)), col='green')
plot(foodexp~income, pch=20, main = "Regresión Robusta")
abline(fitLSR, col='red')

# ¿Qué curva elegir?
summary(fitLS)$r.squared
summary(fitLS)$coefficients[,4]

summary(fitLSR)$r.squared
summary(fitLSR)$coefficients[,4]

# La regresión robusta obtuvo un mejor valor RSaquare, pasando del 0.83 a 0.93 con lo que ofrece un mayor ajuste.
detach(engel)
rm(list = ls())






# EJERCICIO NUMERO 2







library(faraway)
data("chicago")
head(chicago)
attach(chicago)

# preparación de los datos. Aplicamos logaritmos sobre la columna income, eliminamos la columna no usada en la regresión.

datos <- chicago
datos$income <- log(datos$income)
datos <- datos[,-5]
head(datos)

# Halle la matriz de correlaciones para todas las variables y el gráfico de dispersión. 

# graficamos a nivel general...
pairs(datos, pch=20)

# graficamos con correlación paramétrica
library(corrplot)
par(mfrow = c(1,1))
corrplot(cor(datos), type = "lower")
corrplot(cor(datos), method = "number", type = "lower")

# Gráfico resumen con distribucion (diagonal), linea de ajuste y scatter plot bivariable (inferior) y valor de correlacion y significancia (superior).
library(PerformanceAnalytics)
chart.Correlation(datos, histogram=TRUE, pch=19)

# En el gráfico se aprecia una correlación fuerte y positiva entre involact y fire (0.70), involact y race (0.71) y poco significativa pero positiva entre involact y age (0.48)
# Una correlacion fuerte y negativa entre involact e income (-0.70).

# Como se aprecian outliers, realizamos la correlación robusta (en esta caso winsorizada).
library(WRS)
winall(datos)
# Se aprecian datos ligeramente superiores aunque siguen el mismo patrón. También en el caso de involact con income es fuerte y neamtiva (-0.73)
# Se aprecia también mirando p-value, que en todos los casos, p-value < 0.05, por lo que rechazamos la hipótesis nula (correlación es 0). Existe una relación significativa entre las variables.


# Ajuste el modelo de regresión lineal múltiple clásico y robusto.

summary(fitLM0 <- lm(involact ~ ., data = datos))

# El intercepto posee un p-value > 0.05 por lo que lo excluimos al no ser significativo. Volvemos a ajustar...

summary(fitLM <- lm(involact ~ .-1, data = datos))

# Ahora todas las variables son significativas, a excepción de income con p-value > 0.05. Pero le añadimos a la ecuación.
# La ecuación quedaría como: involact = 0.0095*race + 0.0398*fire - 0.010*theft + 0.008 * age - 1.64 * income

# F-statistic (ANOVA), el p-value es < 0.001 nos aporta que el modelo es significativo.

# Se nos muestra en F-statistic del paso anterior. Seria: F(2,42)=101.2, p-value < 0.001
# La bondad de ajuste viene dada por R2 que en este caso es de 0.87, el modelo es muy eficaz (entre 0 y 1). El modelo es significativo y bueno.
anova(fitLM0, fitLM)
# p-valor > 0.05 nos dice que eliminar el intercepto no tiene efectos significativos en el modelo.
avPlots(fitLM)

# Colinealidad:
vif(fitLM)
# Como los datos que arroja la prueba no son mayores de 10, nos dice que los predictores NO están correlacionados. No hay colinealidad.

# Seleccion automática de predictores
summary(step(fitLM, action = F, direction = "forward"))

# En esta caso las variables se incluyen en el modelo final.

# Importancia relativa de los predictores.
library(relaimpo)
(crlm <- calc.relimp(fitLM0, type = c("lmg", "last", "first", "pratt"), rela=T))
plot(crlm)
# En la mayoría de los métodos se aprecia como "race" y "fire" los que más importancia tienen y de ellos "fire" el más imporante. (que además son las de mayor significancia con ***)

# Diagnóstico
par(mfrow = c(2,2))
plot(fitLM)

# Resiuals vs. Fitted.
par(mfrow=c(1,1))
plot(fitLM, which=1)
# Están bastante cercanos a 0, no hay casi dispersion, aparecen 3 outliers.

# Normal Q-Q
plot(fitLM, which=2)
# Se aprecia que se ajustan lso residuos a la curva por lo que existe normalidad. Se muestran los 3 outliers.

# Scale-Location
plot(fitLM, which=3)
# No se aprecian enfecto estilo embudo.

# Resiuals vs. Leverage
plot(fitLM, which=4)
# Se muestra la existencia de los 3 outliers.

# Test de datos atípicos
outlierTest(fitLM)
# Solo está detectando uno, el más alejado el 60610.

# Evaluación de observaciones influyentes.
avPlots(fitLM)

influencePlot(fitLM, id.method = "noteworthy", main="Grafico Influencia")
# Los más influyentes son los outliers 60610 y 60607 aunque con el outlierTest solo detectó el 60610.

# Normalidad de los residuos
qqPlot(fitLM)
# Ya comentado anteriormente que se ajusta bastante a la recta.

# Distribución de los residuos estunderizados
library(MASS)
sresid <- studres(fitLM)
hist(sresid, freq = F, main = "Distribucion de los residuos")
xfitLM <- seq(min(sresid), max(sresid), length=40)
yfitLM <- dnorm(xfitLM)
lines(xfitLM, yfitLM)
# Distribución bastante normal, aunque ligeramente sesgada a izda.

# Evaluación de homocedasticidad
ncvTest(fitLM)
# Se rechaza Ho por lo que no hay varianza constante en los errores.
spreadLevelPlot(fitLM)

# Evaluación de la linealidad
crPlots(fitLM)
# Race, Fire y age son lo que más se mantienen en la linealidad.

# Evaluación de la independencia de errores
durbinWatsonTest(fitLM)
# p-valor al límite 0.05.

# Por existencia de outliers, se debería trabajar  aplciando técnicas robustas. Que es el siguiente paso.




# -------------------------------------------------------------------------------
# AJUSTE USANDO ltsReg PARA ROBUSTOS (Regresión por mínimos cuadrados recortados)
#--------------------------------------------------------------------------------

# Graficamos los outliers
pairs(datos)
out <- influence.measures(fitLM)
plot(fitLM, which=4)
# Existencia de outliers. Lo conveniente es aplicar regresión robusta.

# library(rrcov)
fitLMR0 <- ltsReg(involact ~ race + fire + theft + age + income)
summary(fitLMR0)

# Todos menos income, poseen un p-valor < 0.05 rechazamos la hipótesis nula por lo que debeb aparecer en la ecuación. El intercepto no pues aceptamos la hipótesis nula.

fitLMR <- ltsReg(involact ~ race + fire + theft + age + income - 1)
summary(fitLMR)

# Ahora todas las variables son significativas, a excepción de income con p-value > 0.05. Pero le añadimos a la ecuación.
# Respecto al modelo paramétrico comentar que da menor significancia a la variable "race" e income no aparece con su estimado en netativo.
# La ecuación quedaría como: involact = 0.00458*race + 0.0453*fire - 0.00751*theft + 0.0019 * age - 0.0000059 * income

# F-statistic (ANOVA), el p-value es < 0.001 nos aporta que el modelo es significativo.

# Ajuste global del Modelo y Bondad de Ajuste...
# Se nos muestra en F-statistic del paso anterior. Seria: F(5,37)=72.38, p-value < 0.001
# La bondad de ajuste viene dada por R2 que en este caso es de 0.90, el modelo es muy eficaz (entre 0 y 1). El modelo es significativo y bueno.
# Comentar que R2 el ligeramente inferior al obtenido por el anterior (paramétrico)

# Estimacion de la precisión de los modeloss
# La regresión robusta obtuvo un mejor valor RSquare, pasando del 0.87 a 0.90 con lo que ofrece un mayor ajuste.

detach(chicago)
rm(list = ls())

