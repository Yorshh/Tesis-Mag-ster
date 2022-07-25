library(quantmod)
library(tseries)
library(rmgarch)
library(rugarch)
library(parallel)
library(fImport)
library(fPortfolio)
library(gtools)

#Esto es porque cuando cargué los datos puse que me agregara la columna Vx [x: 1-> 1131]
#names(LogReturns2) <- LogReturns2[1,]    
#LogReturns2 <- LogReturns2[-1,]
#LogReturns2$BIO.P <- NULL       #No se pudo calcular el DCC por problemas de convergencia


LogReturns2$"BIO-B" <- NULL
LogReturns2$"TAP-A" <- NULL
LogReturns2$"WSO-B" <- NULL




omega <- c()
alfa <- c()
beta <- c()
a <- c()
b <- c()
box_test <- c()

fechas <- LogReturns2[2:nrow(LogReturns2),1]    #Elimina la primera fecha para las correlaciones
fechas_2 <- LogReturns2[1:nrow(LogReturns2), 1] #Todas las fechas, para el de los residuales

nombres <- c()    #Para guardar los nombres de los activos que analizo (Lo hare así más que nada para cuando haga pruebas y no utilice todos los activos)
df_correlaciones <- as.data.frame(fechas)

data_residuales <- as.data.frame(fechas_2)  #Dataframe con los residuales estandarizados
ljung_test <- c()                           #Lista para guardar el test Ljung-Box de cada serie

bine_size <- c(20,30,40,50)                 #Samplesizes para el test de Godness of fit (por default, estos son los tamaños)
df_gof <- data.frame()                      #DataFrame con los test de Goodness of fit (adjusted Pearson) para los distintos samplesizes


#Aplicar los modelos para todas las series... Extrae las correlaciones condicionales, los parámetros, los tests y los residuales estandarizados
for (i in 3:5){
  print(colnames(LogReturns2)[i])
  #Guardando el nombre del activo a analizar
  nombres <- append(nombres, colnames(LogReturns2)[i])
  
  
  #seteando x e y
  x <- LogReturns2[,2]           #x es los retornos del sp500
  y <- LogReturns2[,i]           #y es los retornos de cada uno de los otros activos
  
  
  #Uniendolo en un DataFrame
  datos <- as.data.frame(x) #Nota, creo que e
  datos <- cbind(datos, y)
  
  
  
  
  #Especificaciones del modelo GARCH
  spec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                     variance.model =list(garchOrder = c(1,1), model = "sGARCH"),
                     distribution = "sstd")
  
  
  #Primero voy a calcular el GARCH(1,1) Para cada serie, esto lo haré para obtener los LjungBox test y GOF
  if (i == 3){
    garch_fit <- ugarchfit(spec, data = x, fit.control=list(scale=TRUE)) #Primera Iteración, GARCH S&P 500
    residuales_garch <- residuals(garch_fit, standardize = T)          #Residuales Estandarizados
    residuales_garch_2 <- residuales_garch*residuales_garch              #Residuales Estandarizados al Cuadrado
    #Guardando los residuales estandarizados (por si acaso necesito ver algo en python)
    data_residuales <- cbind(data_residuales, residuales_garch)          #Guarda los que NO estan al cuadrado
    test_ljung <- Box.test(residuales_garch_2, lag = 1, type = c("Ljung-Box"), fitdf = 0)$p.value  #Test Ljung-Box
    ljung_test <- append(ljung_test, test_ljung)
    df_gof <- rbind(df_gof, gof(garch_fit, bine_size)[,3])                    #Test Goodness of fit
    
    
    #Los mismo para el otro activo de la primera iteración
    garch_fit <- ugarchfit(spec, data = y, fit.control=list(scale=TRUE))
    residuales_garch <- residuals(garch_fit, standardize = T)          #Residuales Estandarizados
    residuales_garch_2 <- residuales_garch*residuales_garch              #Residuales Estandarizados al Cuadrado
    data_residuales <- cbind(data_residuales, residuales_garch)
    test_ljung <- Box.test(residuales_garch_2, lag = 1, type = c("Ljung-Box"), fitdf = 0)$p.value  #Test Ljung-Box
    ljung_test <- append(ljung_test, test_ljung)
    df_gof <- rbind(df_gof, gof(garch_fit, bine_size)[,3])                    #Test Goodness of fit
    
  }else{
    #Ahora calcula para las otras series 
    garch_fit <- ugarchfit(spec, data = y, fit.control=list(scale=TRUE))
    residuales_garch <- residuals(garch_fit, standardize = T)          #Residuales Estandarizados
    residuales_garch_2 <- residuales_garch*residuales_garch              #Residuales Estandarizados al Cuadrado
    data_residuales <- cbind(data_residuales, residuales_garch)
    test_ljung <- Box.test(residuales_garch_2, lag = 1, type = c("Ljung-Box"), fitdf = 0)$p.value  #Test Ljung-Box
    ljung_test <- append(ljung_test, test_ljung)
    df_gof <- rbind(df_gof, gof(garch_fit, bine_size)[,3])                    #Test Goodness of fit
  }
  
  
  #Especificaciones del modelo DCC
  dcc_spec1 = dccspec(uspec = multispec(replicate(2, spec)), dccOrder = c(1,1),
                      distribution = "mvt")
  
  
  #Ajustando el modelo con los datos en base a las especificaciones establecidas
  
  garchdccfit = dccfit(dcc_spec1, data = datos, fit.control=list(scale=TRUE))
  
  
  #Nota! la siguiente línea entrega la matriz de correlaciones
  #Debo extraer el valor de la correlación entre activos únicamente
  matriz_correlaciones <- rcor(garchdccfit)
  
  correlaciones <- c()
  for(j in c(1:nrow(LogReturns2)-1)){
    correlaciones <- append(correlaciones, matriz_correlaciones[2,1,j])
  }
  
  df_correlaciones <- cbind(df_correlaciones, correlaciones)
  
  omega <- append(omega, as.double(coef(garchdccfit)[8]))
  alfa <- append(alfa, as.double(coef(garchdccfit)[9]))
  beta <- append(beta, as.double(coef(garchdccfit)[10]))
  
  a <- append(a, as.double(coef(garchdccfit)[13]))
  b <- append(b, as.double(coef(garchdccfit)[14]))
  
}


#Nota: A partir de abajo hay varias variables llamadas names, nombres, y así.
#Fue producto del desorden en su momento y la hora jaja (son casi las 6am)
#nombres: todos los activos, sin el SP500
#nombres_prueba: Nombre de todo... fecha + sp500 + el resto
#nombres_prueba_2: lo mismo PERO SIN FECHA
#names: Fecha + el resto... SIN EL SP500



nombres_prueba <- c("Date", "SP500")
nombres_prueba <- append(nombres_prueba, nombres)
colnames(data_residuales) <- nombres_prueba
data_residuales



df_ljung_test <- data.frame()
df_ljung_test <- cbind(ljung_test)
nombres_prueba_2 <- c("sp500")
nombres_prueba_2 <- append(nombres_prueba_2, nombres)
df_ljung_test <- cbind(nombres_prueba_2, df_ljung_test)
df_ljung_test



str_bine_sizes <- c("20", "30", "40", "50")
colnames(df_gof) <- str_bine_sizes
df_gof <- cbind(nombres_prueba_2, df_gof)
df_gof




resultados <- data.frame(nombres)
resultados <- cbind(resultados, omega)
resultados <- cbind(resultados, alfa)
resultados <- cbind(resultados, beta)
resultados <- cbind(resultados, a)
resultados <- cbind(resultados, b)



#Me falto agregar la columna "Date"
names <- c("Date")
names <- append(names, nombres)
colnames(df_correlaciones) <- names




write.csv(resultados, "coeficientes_dccgarch_final.csv")
write.csv(df_correlaciones,"sp_correlation_final.csv")
write.csv(data_residuales,"data_residuales.csv")
write.csv(df_ljung_test,"ljung_test.csv")
write.csv(df_gof,"goodness_test.csv")


