#------------------------------------------------------------------------#
# Máster en Estadística Aplicada con R software
# Máxima Formación
# TEMA 6. MODELOS AVANZADOS
# marzo 2015
#------------------------------------------------------------------------# 

    # - Modelo aditivo (AM): modelan la relación no-lineal entre X e Y
    # - Modelo lineal generalizado (GLM): modelan los errores (no necesariamente normal)
    # - Modelo mixto (MM): modelan la varianza e independencia, para 
    # datos anidados y estructuras de correlación.
    

#------------------------------------------------------------------------#
## MODELO ADITIVO (AM)
#------------------------------------------------------------------------#
    # Utilizaremos los datos de bioluminiscencia de organismos pelágicos en un 
    # gradiente trazado en el noreste del Océano Atlántico. Queremos modelar la 
    # relación bioluminiscencia-profundidad. Solo trabajaremos con la estación 16.
    setwd("~/Google Drive/Master Estadistica Aplicada/Tema 6/Codigo")
    ISIT<-read.table("ISIT.txt",header=T)
    head(ISIT)
    str(ISIT)
    
    # Graficamos los datos
    Sources16<-ISIT$Sources[ISIT$Station==16]
    Depth16<-ISIT$SampleDepth[ISIT$Station==16]
    plot(Depth16,Sources16,type="p")
    # Bioluminiscencia en función de la profundidad.
    # A menor profundidad, mayor numero de organismos con luminiscencia.
    # Se aprecia que NO es lineal. Se aplican modelos aditivos. Para ver que tipo de curva se ajusta mejor a los datos.
    
    # Se introducen los datos en formato fórmula ("y~x") 
    # "s(x)" función de suavizado para "x". La variable a suavizar. 
    # "fx = FALSE, k=-1" no fijamos de antemano ningún valor para el suavizado 
    # sino que éste se estima mediante validación cruzada. 
    # "bs = "cr"" ajustamos una regresión spline cúbica. hay más para probar, pero no mejoran mucho, dice...
    # detach("package:gam")
    library(mgcv)
    M3 <- gam(Sources16 ~ s(Depth16, fx = FALSE, k=-1,bs = "cr"))
    summary(M3)
    
    # Solo temenos el intercepto (pues hay una sola variable). Es significativo.
    # La devianza explicada (R2) es bastante alta: 97.5%, la varianza de los residuales es 7.0236, 
    # el término de suavizado para la profundidad es significativa al nivel del 5%, 
    # los grados de libertad (u orden del polinonio efectivo) estimados para el término de suavizado son 8.817 
    # (son los que  determinan la curvatura de la función; 1 para una línea recta, 
    #  4 para una curva similar a un polinomio cúbico). En esta caso tenemos un polinonio de mayor grado.
    
    # NOTA: La prueba para el término de suavizado (p-valor) no nos dice si el término 
    # es significativo, sino si el término se debe modelar con una función de suavizado. En este caso sí.
    
    anova(M3)
    # Esta salida es útil cuando tenemos variables nominales con más de 2 niveles 
    # ya que nos da un p-valor para todos los niveles utilizando el estadístico F.
    
    # Graficamos el ajuste del suavizado.
    op=par(mfrow=c(1,2))
    plot(M3, se = TRUE)
    
    # Hacemos otro gráfico con los valores ajustados o las prediciones.
    M3pred <- predict(M3, se = TRUE, type = "response")
    plot(Depth16, Sources16, type = "p")
    I1 <- order(Depth16)
    lines(Depth16[I1], M3pred$fit[I1], lty=1)
    lines(Depth16[I1], M3pred$fit[I1]+2*M3pred$se[I1],lty=2)
    lines(Depth16[I1], M3pred$fit[I1]-2*M3pred$se[I1],lty=2)
    par(op)
    
    # Diagnóstico gráfico del modelo. Evaluamos:  
    # - normalidad (QQ-plot)
    # - homogeneidad (residuales vs valores ajustados)
    # - el ajuste del modelo (valores ajustados vs valores observados)
    # El diagnóstico del modelo se lleva a cabo en esta opción:
    par(mfrow = c(2,2))
    gam.check(M3)
    
    # Los errores se desvían de la normalidad.
    # Parece haber un poco de heterogeneidad en el predictor. (arriba dcha.)
    # Puede ser que si añadimos más variables explicativas mejore, o interacciones o modelos más complejos.
    # No parece haber un buen ajuste, tal vez necesitamos incluir más variables 
    # explicativas, interacciones o modelos más complejos como el modelo aditivo mixto.
    
    # Vamos a comparar nuestro modelo con el modelo de regresión lineal simple y 
    # el modelo de regresión de un polinomio de tercer grado:
    M1 <- lm(Sources16 ~ Depth16)
    M2 <- lm(Sources16 ~ poly(Depth16, 3))
    AIC(M1, M2, M3)
    # Un valor de 2 de diferencia ya nos haría decantarnos por él.
    # Vemos que el menor AIC corresponde al modelo AM.
    
    # NOTA: si queremos ajustar un GAM basta con especificar la familia de la 
    # distribución de los errores (e.g. para una distribución de Poisson 
    # gam(y ~ s(x), family = poisson).
                       
    
#------------------------------------------------------------------------#    
# MODELO LINEAL GENERALIZADO (GLM)    
#------------------------------------------------------------------------#  
    
    
    #------------------------------------------------------------------------# 
    ### Modelo de Poisson
    #------------------------------------------------------------------------# 
    # - Para datos de conteo (enteros positivos). Errores distintos de la distr. normal.
    # Ejemplo: número de individuos, número de empresas, etc..
    
    # Utilizaremos los datos de atropellos de anfibios (variable respuesta) en 
    # 52 sitios a lo largo de una carretera en Portugal (Zuur et al. 2013). 
    # La base de datos contiene 16 variables explicativas:  
    #   
    # - OPEN.L, terrenos abiertos (ha)
    # - OLIVE, olivares (ha)
    # - MONT.S, Dehesa con arbustos (ha) 
    # - MONT, Dehesa sin arbustos (ha) 
    # - POLIC, policultivo (ha)
    # - SHRUB, Arbustos (ha)
    # - URBAN, Urbano (ha)
    # - WAT.RES, depósitos de agua (ha)
    # - L.WAT.C, longitud de los cursos de agua (km)
    # - L.D.ROAD, longitud de la carretera sucia (m)
    # - L.P.ROAD, longitud de la carretera asfaltada (km)
    # - D.WAT.RES, Distancia a los depósitos de agua 
    # - D.WAT.COUR, Distancia a cursos de agua
    # - D.PARK, Distancia al parque natural (m) 
    # - N.PATCH, Número de parches de hábitat
    # - P.EDGE, bordes del perímetro
    # - L.SDI, índice de diversidad de Shannon del Paisaje 
    
    
    # Comenzaremos considerando tan solo el efecto de la variable D.PARK. 
    # Queremos saber cómo es la relación del número de atropellos y la 
    # distancia al parque natural (variable explicativa D.PARK). 
    
    # Observamos los datos
    RK <- read.table("RoadKills.txt", header=T)
    head(RK)
    str(RK)
    
    # Graficamos los datos
    plot(RK$D.PARK, RK$TOT.N, 
     xlab = "Distancia al parque", ylab = "Atropellos")
    # Vemos que decae, además es muy dispersa para distancias cortas.
    # La relación parece no ser lineal y además hay variación. Se ve bien entre 0 y 10.000 que está muy dispersa.
    # la variación es mayor al aumentar el número de atropellos (aumenta la varianza al aumentar la media). Esto parece indicar que el GLM Poisson es la mejor opción.
    
    # Modelo 
    M1 <- glm(TOT.N ~ D.PARK, family = poisson, data = RK) 
    summary(M1)  
    # - el modelo ajustado.
    # - informacion sobre los residuales.
    # - estimaciones para el intercepto (4.31) y la pendiente (-0.000106).
    # - el estadístico z y su p-valor para la inferecia sobre los parámetros (H0=betai=0).
    # - la devianza del modelo nulo (D0, que solo incluye el intercepto) y de los residuales del modelo ajustado (Dr), que nos dan información sobre la bondad de ajuste del modelo. Cuanto menor sea la devianza de los residuales mejor es el modelo.
    # - el AIC, útil para la selección de modelos. Cuanto menor sea el AIC mejor es el modelo.
    
    # Para estimar un pseudo-R2 del GLM podemos utilizar la devianza explicada: 100*(D0-Dr)/D0
    
    (100*(1071.4-390.9)/1071.4)
    
    # Sale cercana al 64% (63.51) que no está nada mal.
    # La distancia al parque explica el 63.51% de la variabilidad en los atropellos. 
    
    # Graficamos el modelo ajustado.
    MyData=data.frame(D.PARK=seq(from=0,to=25000,by=1000))
    G<-predict(M1,newdata=MyData,type="link",se=T)
    F<-exp(G$fit)
    # intervalos de confianza
    FSEUP<-exp(G$fit+1.96*G$se.fit)
    FSELOW<-exp(G$fit-1.96*G$se.fit)
    plot(RK$D.PARK, RK$TOT.N, 
     xlab = "Distancia al parque", ylab = "Atropellos")
    lines(MyData$D.PARK,F,lty=1)
    lines(MyData$D.PARK,FSEUP,lty=2)
    lines(MyData$D.PARK,FSELOW,lty=2)
    # Modelo se comporta bastante bien.
    # Intervalos de confianza más amplios a menores distancias.
    # Vemos que los intervalos de confianza al 95% se ensanchan para distancias cortas, donde algunos puntos se escapan. La curva es exponencial. 
    
    # Como quedan bastantes valores por fuera, indica que tal vez con otra variable explicativa a mayores. El modelo mejoraría.
    
    # Selección de variables
    
    # De las 16 variables de la base de datos solo utilizaremos algunas
    # Vamos a transformar las variables que presentan valores muy altos lo que nos llevaría a problemas de variabilidad.
    # En este caso la transformación seria con raíz cuadrada.
    # Luego, ajustamos el modelo GLM Poisson sin interacciones (suponemos que no hay...)
    RK$SQ.POLIC<-sqrt(RK$POLIC)
    RK$SQ.WATRES<-sqrt(RK$WAT.RES)
    RK$SQ.LPROAD<-sqrt(RK$L.P.ROAD)
    RK$SQ.SHRUB<-sqrt(RK$SHRUB)
    RK$SQ.DWATCOUR<-sqrt(RK$D.WAT.COUR)
    
    M2<-glm(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
      SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
      SQ.DWATCOUR+D.PARK,family=poisson,data=RK)
    summary(M2)
    
    # Algunas variables explicativas no son significativas, 
    # por lo que tenemos que seleccionar el mejor modelo. Tenemos varias opciones:
    
    # - selección automática mediante la función "step(M2)".
    # - selección manual; eliminar de uno en uno aquellos términos menos significativos 
    # según el estadístico z. Luego, utilizar la función "drop1(M2, test="Chi")" para 
    # realizar la prueba de devianza uno en uno, o la función "anova(M2)" para la prueba 
    # de devianza secuencial.
    
    # Opción "step" -------------------------------
    summary(step(M2, trace=0)) # trace=0 para que no nos diga todos los pasos intermedios.
    # o
    # library(MASS); summary(stepAIC(M2, trace=0))
    # Ambos modelos propuestos son iguales al anterior M2.
    
    # Opción "drop" -------------------------------
    drop1(M2,test="Chi")
    # Los p-valores no son los mismos que los del comando "summary(M2)" porque 
    # ambas pruebas son aproximadas y utilizan estadísticos distintos. 
    
    # Siguiendo los resultados de "drop" eliminaríamos el menos significativo, la variable SQ.DWATCOUR.
    M2b<-glm(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
        SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
        D.PARK,family=poisson,data=RK)
    summary(M2b)
    # Ahora todos los términos son significativos. Les incluimos en el modelo.
    
    # Opción "anova" -------------------------------
    # Si se tienen dudas de qué opción seguir, podemos utilizar la función "anova(M2)", 
    # que tiene la ventaja de que también nos da p-valores cuando analizamos variables 
    # nominales.
    # Nuevamente comenzamos evaluando el término menos significativo, el SQ.DWATCOUR. 
    # Como ya hemos ajustado un modelo quitando este término, vamos a aplicar el anova directamente.
    anova(M2, M2b, test = "Chi")
    # No existen diferencias significativas entre ambos modelos H0 es "no existen diferencias" (p=0.132, Dif deviance=2.269) 
    # por lo que nos quedamos con el modelo más parsimonioso (el que tenga menos variables), el que no incluye la variable SQ.DWATCOUR.
    
    # A continuación, para evaluar el término OPEN.L (que no era del todo significativo en el
    # modelo inicial) debemos quitarlo del modelo y volver a ajustarlo.
    M3 <- glm(TOT.N ~ MONT.S + SQ.POLIC + D.PARK +
    SQ.SHRUB + SQ.WATRES + L.WAT.C + SQ.LPROAD +
    SQ.DWATCOUR, family = poisson, data = RK)
    # se compara en el resultado de la anterior comparación M2 y M2b)
    anova(M2, M3, test = "Chi")
    # Aquí sí detectamos una diferencia significativa (marginal) por lo que no quitamos este término del modelo. 
    # Como mencionamos anteriormente, debemos chequear posibles problemas de sobredispersión.
    
    
    # Sobredispersión ("overdispersion")
    #   La sobredispersión ocurre cuando la varianza es mayor que la media: 
    #   phi=D/(n-p) con n el número de casos y p el número de parámetros en la regresión 
    #  >1 indica sobredispersión.  Un problema. Menos o cercana a 1, bueno, no hay sobredispersión.
    #   Para el modelo M2 sería 
    (270.23/42) # que da 6.43
    
    # Las causas de la sobredispersión pueden ser:
    # - aparentes, debido a la mala especificación del modelo (e.g. faltan covariables o interacciones, o debido a la presencia de *outliers* en la respuesta, o a efectos no lineales de las covariables, o a elegir erróneamente la función link).
    # - real, debido a que la variación de los datos realmente es mayor que la media (también puede deberse a que hay muchos ceros, o que hay observaciones agrupadas o a la correlación entre observaciones).
    # 
    # El procedimiento a seguir es agregar covariables e interacciones (Lo más facil). 
    # Si esto no funciona podemos utilizar un GLM cuasi-Poisson (se cambia en la familia) o un GLM binomial negativo.
    
    
    # Validación del modelo (aunque hay sobredispersión y no se seguiría el proceso...)
    
    # Vamos a graficar los residuos que seleccionemos en función de: 
    # - cada variable explicativa incluida en el modelo, si encontramos un patrón deberíamos pensar en incluir términos cuadráticos o utilizar GAM o concluir que estamos incumpliendo el supuesto de independencia.
    # - cada variable explicativa no incluida en el modelo, si encontramos un patrón deberíamos incluir esta variable explicativa al modelo.
    # - en función del tiempo y/o coordenadas espaciales, si encontramos un patrón deberíamos concluir que estamos incumpliendo el supuesto de independencia.
    # - en función de los valores ajustados, si encontramos un patrón de dispersión estaríamos ante un caso donde hay sobredispersión o hemos elegido mal la distribución.  
    
    # gráficos por defecto
    op <- par(mfrow=c(2,2))
    plot(M2b)
    par(op)
    
    # Aumento de la variabilidad (embudo) arriba izda.
    # Diferencias en la normalidad.
    
    # para obtener los distintos residuales
    EP=resid(M2b,type="pearson")
    ED=resid(M2b,type="deviance")
    
    # como la función resid no considera la sobredispersión
    # debemos dividir los residuales de Pearson por la raíz cuadrada 
    # del valor del valor de sobredispersión (272.5/43=6.337209)
    mu=predict(M2b,type="response")
    E=RK$TOT.N-mu
    EP2=E/sqrt(6.337209*mu)
    
    op <- par(mfrow=c(2,2))
    plot(x=mu,y=E,main="Response residuals")
    plot(x=mu,y=EP,main="Pearson residuals")
    plot(x=mu,y=EP2,main="Pearson residuals scaled")
    plot(x=mu,y=ED,main="Deviance residuals")
    par(op)
    
    # No deberíamos encontrar ningún patrón en estos datos pero los hay (agrupamiento de datos a la izda...), 
    # seguramente debido a los problemas de sobredispersión que hemos detectado previamente. 
    
    # Posible solución: modelo binomial negativo.
    
  #------------------------------------------------------------------------# 
  ### Modelo Binomial negativo
  #------------------------------------------------------------------------# 
    # - Para datos de conteo (enteros positivos) con mayor variabilidad que la 
    # permitida por errores con distribución Poisson.
    
    # Continuamos con el ejemplo anterior y volvemos a partir con las 11 variables explicativas.
    # Utilizaremos la función "glm.nb" del paquete "MASS".
    library(MASS)
    M6<-glm.nb(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
     SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
     SQ.DWATCOUR+D.PARK,link="log",data=RK) # con la función link="log"
    summary(M6, cor=FALSE)
    # Vemos que ahora la mayoría de los parámetros no son significativos al nivel de significación del 5%, por lo que necesitamos utilizar técnicas de selección de variables. 
    # También observamos la información del parámetro "theta" de la distribución binomial. 
    
    
    # Opción step ------------------
    summary(step(M6, trace=0))
    summary(stepAIC(M6, trace=0))
    # Ambos seleccionan las variables OPEN.L, L.WAT.C, SQ.LPROAD y D.PARK , 
    # pero SQ.LPROAD no es significativa y L.WAT.C es marginalmente significativa. 
    
    # Opción drop -------------------
    M7<-glm.nb(TOT.N~OPEN.L+L.WAT.C+SQ.LPROAD+D.PARK,link="log",data=RK)
    drop1(M7,test="Chi")
    # Quitaríamos ambas variables (recuerda que los p-valores son aproximados y 
    # 0.044 está cerca de nuestro 0.05 de referencia).
    
    # Opción anova ----------------
    M7b<-glm.nb(TOT.N~OPEN.L+L.WAT.C+D.PARK,link="log",data=RK)
    M7c<-glm.nb(TOT.N~OPEN.L+D.PARK,link="log",data=RK)
    anova(M7b, M7c)
    # La significación es 0.02, no demasiado significativo, por lo que decidimos 
    # eliminar a ambas del modelo (pero podríamos dejar a L.WAT.C).
    
    # Sobredispersión
    summary(M7c)
    (51.839/49) # 1.05
    # Es cercana a 1, por lo que hay poca sobredispersión.
    
    
    # Validación del modelo
    M8<-glm.nb(TOT.N~OPEN.L+D.PARK,link="log",data=RK)
    summary(M8)
    drop1(M8,test="Chi")
    # Ambas variables son significativas.
    
    par(mfrow=c(2,2))
    plot(M8)
    # No observamos problemas. Salvo algún valor atípico (outliers)
    
    # Interpretación
    # OPEN.L y D.PARK poseen negativos, por lo que....
    # Menor número de atropellos de anfibios en zonas abiertas y a menores distancias del parque natural.
    
    # Cuanto más lejos del parque (D.PARK), menor número de atropellos a anfibios 
    # (TOT.N). Además, los terrenos abiertos (OPEN.L) disminuyen el número de atropellos.
    
    
  #------------------------------------------------------------------------#   
  ### Modelo logístico (binaria)
  #------------------------------------------------------------------------# 
    # - Para datos de binarios (0 y 1) como los datos de presencia y ausencia o 
    # pruebas con dos resultados posibles.
    
    #### Odds ratio
    # odds (O, u "oportunidad"): el número de casos en los que el evento ocurre 
    # dividido por el número de casos en que no ocurre).
    
    # odds ratio (OR, también llamado "oportunidad relativa" o "tasa de oportunidad relativa"): 
    #   OR=O2/O1. Si el OR=1 entonces no hay factor de riesgo. No hay diferencia entre las dos condiciones.
    
    
    # Utilizaremos los datos "Boar" de Zuur et al. 2013 para evaluar la probabilidad 
    # que tienen los osos en el sur-centro de España de contraer tuberculosis (Tb) 
    # en función su tamaño (LengthCT, distancia cabeza-tronco) y edad (age). 
    Boar<-read.table("Boar.txt",header=T)
    head(Boar)
    str(Boar)
    summary(Boar)
    # Vemos que hay valores perdidos y variables que pasar a factor porque están como int (sexo).
    
    # Primero exploramos los datos gráficamente.
    op<-par(mfrow=c(1,2))
    boxplot(LengthCT~Tb, data=Boar, ylab="LengthCT", xlab="Tb")
    plot(Tb~LengthCT, data=Boar, ylab="LengthCT", xlab="Tb")
    abline(lm(Tb~LengthCT, data=Boar))
    par(op)
    # Observamos que los osos con tuberculosis presentan tamaños mayores (suelen ser más altos). 
    # También graficamos un modelo lineal simple para observar los inconvenientes de este tipo de ajuste.
    # Se ve que no es lineal, es tipo Sí, No.
    
    # Modelo
    
    #por defecto binomial(link="logit")
    B1=glm(Tb~LengthCT, family = binomial, data = Boar)
    summary(B1)
    
    # En el GLM Bernoulli no tenderemos problemas de sobredispersión. No calculamos el coeficiente.
    # Tenemos los valores de la devianza y el AIC.
    # El intercepto y la pendiente valen -3.89 y 0.032, respectivamente, y ambos 
    # son significativos al nivel del 5%. 
    
    # El modelo ajustado sería: logit(phi) = -3.89 + 0.032 * LengthCT
    
    # Esto significa que para una unidad de aumento en el tamaño de los osos 
    # aumenta 0.032 el log odds de la tuberculosis.
    
    # Si en lugar de utilizar los log odds queremos interpretar los odds ratio debemos aplicar la exponencial:
    exp(coef(B1)) #odds ratio
    exp(confint(B1)) #IC
    # Para una unidad de aumento en el tamaño de los osos aumentaría un uno la ventaja relativa de tener la enfermedad.
    
    # Graficamos el modelo. 
    range(Boar$LengthCT, na.rm=T)
    MyData <- data.frame(LengthCT = seq(from = 46.5, to = 165, by = 1))
    Pred <- predict(B1,newdata = MyData, type = "response")
    plot(x=Boar$LengthCT, y = Boar$Tb,xlab="Length", ylab="Tb")
    lines(MyData$LengthCT,Pred)
    # Como se trata de datos que solo contienen ceros y unos, vemos dos bandas de puntos.
    # Esta no sería la forma más fácil de validarlo. Pero mostramos los gráficos de validación.
    par(mfrow=c(2,2))
    plot(B1)
    
    summary(B1)
    drop1(B1, test="Chi")
    # Con la función  drop1 vemos que el término LengthCT es significativo. Como hay un solo término no haría falta. Mejor cuando hay más.
    
    # VALIDACIÓN DE ESTOS TIPOS DE MODELOS.
    
    #   Una forma de obtener la significación de todo el modelo es probar la hipótesis nula de que el modelo no es mejor (en términos del likelihood) que el modelo nulo (un modelo ajustado con solo el término intercepto, donde todos los beta valen cero).
    # 1-pchisq(D0-Dr, df0-dfr)
    # D0 es la devianza del modelo nulo y Dr la del modelo ajustado o residual, Df los grados de libertad (nulo y ajustado).
    
    1 - pchisq(700.76-663.56, 507-506)
    # El modelo es significativo.
    
    # Bondad de ajuste pseudo-R2 de Nagelkerke.
    library(rms)
    lrm(Tb~LengthCT, data=Boar)
    
    

#------------------------------------------------------------------------#
## Modelo lineal mixto (LMM)
#------------------------------------------------------------------------#
    # - modelar la varianza en lugar de asumir que es constante y 
    # - modelar las observaciones correlacionadas. 
    # - Combinan efectos fijos (factores fijas determinados por el investigador) y aleatorios
    # - Conclusiones para fijos se van a remitir solo a los niveles que estamos evaluando...
    # y para aleatorios: son respecto a toda la población o posibilidades para toda la variable tratamiento.
    
    
    #### Protocolo 
    #     1)  Estructura aleatoria óptima. 
    #       - construir un modelo saturado (fijo). Con todas las variables.
    #       - comparar modelos con distinta estructura aleatoria (REML y anova, AIC o BIC).
    #     2)  Comparar modelos con Estructura fija óptima (ML).
    #     3)  Ajuste del modelo final con REML. Interpretar y validar. 
    
    #     Utilizaremos datos sobre la fauna bentónica de sistemas intermareales de 
    #     la costa holandesa en el Instituto *RIKZ* en 2002 (Zuur et al. 2007). 
    #     Se desea estudiar el cambio en el nivel del mar sobre la fauna de los 
    #     fondos costeros. Se recogieron datos del macrobentos (riqueza o número 
    #     de especies, S) y de variables abióticas (e.g. altura del sitio de muestreo 
    #     respecto a la marea, NAP) en 9 playas donde en cada una se seleccionaron 
    #     5 puntos de muestreo. 
    
    RIKZ<-read.table("RIKZ.txt", header=T)
    head(RIKZ)
    str(RIKZ)
    
    # Pasamos los factores a factores.
    RIKZ$fExp<-factor(RIKZ$Exposure)
    table(RIKZ[,"fExp"])
    #     Toma 3 valores 8, 10 y 11 con frecuencias 5, 20 y 20, respectivamente. 
    #     Como la categoría 8 tiene pocos casos vamos a unirla a la categoría 10.
    RIKZ$fExp[RIKZ$fExp==8]<-10
    # Eliminamos el factor 8 para que no aparezca.
    RIKZ$fExp<-droplevels(RIKZ$fExp)
    RIKZ$fBeach<-factor(RIKZ$Beach)
    
    #     1)  Estructura aleatoria óptima.
    # Creamos el modelo saturado (todas las variables, el intercepto (de ahí el 1) y todas las interacciones)
    # Comparamos un modelo sin componente aleatoria.
    m1<-gls(Richness~1+NAP*fExp, method="REML", data=RIKZ)
    # Otro que considera las playas.
    m2<-lme(Richness~1+NAP*fExp, method="REML", data=RIKZ, random=~1|fBeach)
    # Otro que considera dentro de las playas, la altura de la marea.
    m3<-lme(Richness~1+NAP*fExp, method="REML", data=RIKZ, random=~NAP|fBeach)
    
    AIC(m1,m2,m3)
    anova(m1,m2,m3)
    #     El mejor modelo es m2. (Ya tenemos seleccionada la estructura aleatoria)
    
    #     2)  Estructura fija óptima. (En random ponemos lo que pusimos en el m2)
    # Solo con NAP (la parte fija)
    m2a<-lme(Richness~1+NAP, data=RIKZ, random=~1|fBeach, method="ML")
    # NAP y Exposición
    m2b<-lme(Richness~1+NAP+fExp, data=RIKZ, random=~1|fBeach, method="ML")
    # NAP, Exposición y también con la interacción
    m2c<-lme(Richness~1+NAP+fExp+NAP:fExp, data=RIKZ, random=~1|fBeach, method="ML")
    # Sólo con la Exposición
    m2d<-lme(Richness~1+fExp, data=RIKZ, random=~1|fBeach, method="ML")
    
    # Hay que tener cuidado, pues el ANOVA es secuencial como vimos en otro capítulo.
    anova(m2a,m2b,m2c)
    anova(m2d,m2c)
    #     El mejor modelo es el m2c. El que posee todas las variables y todos los efectos de interacción.
    
    
    #     3)  Ajuste del modelo final con REML 
    m2c_final<-lme(Richness~1+NAP+fExp+NAP:fExp, data=RIKZ, random=~1|fBeach, method="REML")
    summary(m2c_final)
    #     Todos los términos son significativos.
    
    #      Validación 
    #     Debemos detectar homogeneidad de varianza (en residuos vs. predicciones), 
    #     ausencia de outliers (se considera outlier un residuo estandarizado cuyo 
    #     valor es superior a 2 en valor absoluto) y ausencia de patrones al graficar 
    #     los residuos frente a las variables explicativas.
    plot(m2c)
    Res<-residuals(m2c_final, type="normalized")
    Fit<-fitted(m2c_final) #level=1
    
    # Residuos vs. predicciones
    op<-par(mfrow=c(2,2))
    plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted")  
    abline(h=0)
    # Vemos cierta dispersión en los datos (embudo), falta de variabilidad, aunque debido a los valores superiores que se marchan.
    # Eliminando dichos valores, se ajusta mejor. Queda entre 2 y -2 por lo que estaría tan mal.
    
    # Residuos vs. predictores (variables explicativas)
    plot(Res~RIKZ$fBeach, xlab="Beach", ylab="Residuals")
    abline(h=0)
    plot(Res~RIKZ$fExp, xlab="Exposure", ylab="Residuals")
    abline(h=0)
    plot(Res~RIKZ$NAP, xlab="NAP", ylab="Residuals")
    abline(h=0)
    # En el boxplot de arriba sí parece haber un efecto de Playa.
    # Boxplot de abajo se ve efecto de exposición.
    # NAP a excepción de los dos valores que se marchan tb. parece estar bien.
    
    # Normalidad de los residuos
    hist(Res)
    qqnorm(Res)
    qqline(Res)
    # Se ve un problema de normalidad en el Q-Q Plot
    
    par(op)
    
    #     Detectamos algo de heterocedasticidad, podríamos utilizar un modelo 
    #     mixto generalizado (GLMM) con distribución de errores de Poisson. 
    #     Además, algunas observaciones pueden ser posibles outliers (atípicos) 
        
    #### Bondad de ajuste R2
    #     R2 de Nakagawa & Schielzeth (2013). 
    #     R2 marginal que solo tiene en cuenta el componente fijo  
    #     R2 condicional que tiene en cuenta ambos componentes 
    #     la diferencia entre ambos indica la variabilidad en los efectos aleatorios.
    library(lme4)
    m0<-lmer(Richness~1+(1|fBeach), data=RIKZ)
    mF<-lmer(Richness~1+NAP+fExp+NAP:fExp+(1|fBeach), data=RIKZ)
    library(piecewiseSEM)
    sem.model.fits(list(m0,mF))
    # Lo que nos dice que nuestro modelo explica un 66% la variabilidad.
    
#------------------------------------------------------------------------#
## Modelo lineal mixto generalizado (GLMM)
#------------------------------------------------------------------------#
    #     Continuamos con nuestro ejemplo de playas.
    # Usamos poisson debido que la variable riqueza es una variable de conteo.
    
    GLMM1 <- glmer(Richness ~ NAP + (NAP|Beach), data = RIKZ, family=poisson)
    summary(GLMM1)
    # Vemos que el intercepto y el NAP son significativos.
    
    # Graficamos los residuos.
    plot(GLMM1)
    # Hay algún dato que se escapa, pero no vemos un patrón claro.
    
    # Realizamos más gráficos...
    Res <- residuals(GLMM1, type="pearson")
    Fit <- fitted(GLMM1)
    
    op<-par(mfrow=c(2,2))
    plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals",
         main="Residuals vs. fitted")
    abline(h=0)
    plot(Res ~ RIKZ$NAP, xlab="NAP", ylab="Residuals", main = "NAP")
    abline(h=0)
    hist(Res, main="Histogram of residuals", xlab="Residuals")
    qqnorm(Res)
    qqline(Res)
    par(op)
    # Se ve una desviación de la normalidad (abajo), y algunos residuos que se salen
    # aunque no hay un patrón claro.
    # Como es poisson, los de abajo no nos interesan, solo los de arriba!!.
    
    AIC(GLMM1)
    AIC(m2c)
    # Vemos que el AIC de este modelo es bastante menor que el del anterior (el no generalizado).
    
    sem.model.fits(list(m2c,GLMM1))
    
    #     Observamos una pequeña mejoría en los residuos y los valores de AIC 
    #     también son mejores.(El anterior 66 y este un 81%)
    
#------------------------------------------------------------------------#    
## Modelo aditivo mixto generalizado (GAMM)
#------------------------------------------------------------------------#
    #     Continuamos con nuestro ejemplo de playas.
    
    # Aditivo porque vamos a considerar que NAP se puede modelar con un función de suavizado.
    # Solo vamos considerar efectos principales no interacciones.
    library(mgcv)
    GAM1<-gamm(Richness~s(NAP)+fExp, random=list(fBeach=~1), data=RIKZ)
    summary(GAM1$lme)
    summary(GAM1$gam)
    # Correspondería a una recta, da 1.
    anova(GAM1$gam)
    
    #     Graficamos. 
    plot(GAM1$gam)
    plot(GAM1$lme)
    
    Res<-resid(GAM1$lme, type="normalized")
    # Vemos el ajuste de la recta y los residuos que alguna se escapa, pero el resto están bien.
    
    plot(Res~RIKZ$Richness)
    # En este sí se ve un patrón.
    
    boxplot(Res~RIKZ$fExp)
    # En el boxplot no se ve tanta diferencia respecto a la playa.
    
    AIC(GAM1$lme)
    AIC(m2c)
    # El AIC es ligeramente mayor al del glm mixto.
    #     El suavizado indica que la relación es lineal (Ref.df=1) pero que 
    #     no mejora el AIC del modelo. No mejora mucho haciendo el suavizado.
    
    #     Intentemos ahora un modelo aditivo mixto generalizado es decir que sea poisson y además de tener suavizado.
    library(mgcv)
    GAMM1<-gamm(Richness~s(NAP)+fExp, random=list(fBeach=~1), data=RIKZ, family=poisson)
    summary(GAMM1$lme)
    summary(GAMM1$gam)
    # Todos los terminos son significativos.
    # El suavizado es significativo y con grado entre 2 y 3. Entre un polinomio de segundo y tercer grado.
    
    anova(GAMM1$gam)
    #     El término paramétrico es significativo: hay un efecto de exposición de las playas.
    #     El término negativo de fExp11 (-0.8752) indica que las playas de 
    #     exposición 11 tienen menor riqueza que las de exposición 10 (el intercepto). 
    #     El término de suavizado con NAP es significativo. 
    #     La varianza del intercepto aleatorio vale 0.27^2.
    
    #     Graficamos el nuevo modelo. 
    plot(GAMM1$gam)
    plot(GAMM1$lme)
    
    Res<-resid(GAMM1$lme, type="normalized")
    # Solo hay un punto en la esquina superior que se escapa. Los demás no enseñan un patrón demasiadamente claro.
    
    plot(Res~RIKZ$Richness)
    boxplot(Res~RIKZ$fExp)
    # En el 10 recordamos que está el 8 incluido. Por lo que hay algo más de variabilidad.
    
    AIC(GAMM1$lme)
    AIC(m2c)
    #     El modelo ha mejorado considerablemente, el AIC ha disminuido y  
    #     solo hay un residuo >|2|.
    # Un aditivo mixto generalizado mejora bastante respecto al mixto.
    