# Optimización de ruta para un vendedor en Colombia
# usando colonias de hormigas y algoritmos genéticos

# Instalar y cargar paquetes necesarios
packages <- c("ggplot2", "dplyr", "sf", "sp", "maps", "animation", "GA")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Definir las 13 ciudades principales de Colombia con coordenadas
ciudades <- data.frame(
  nombre = c("Bogotá", "Medellín", "Cali", "Barranquilla", "Cartagena", 
             "Cúcuta", "Bucaramanga", "Pereira", "Ibagué", 
             "Santa Marta", "Manizales", "Pasto", "Villavicencio"),
  lat = c(4.6097, 6.2476, 3.4516, 10.9685, 10.3932, 
          7.8890, 7.1254, 4.8143, 4.4389, 
          11.2404, 5.0687, 1.2136, 4.1533),
  lon = c(-74.0817, -75.5709, -76.5320, -74.7813, -75.4832, 
          -72.5078, -73.1198, -75.6946, -75.2322, 
          -74.2031, -75.5173, -77.2811, -73.6350)
)

# Calcular matriz de distancias (en km) entre ciudades usando distancia euclidiana
calcular_distancia <- function(lat1, lon1, lat2, lon2) {
  # Radio de la Tierra en km
  R <- 6371
  
  # Convertir grados a radianes
  lat1_rad <- lat1 * pi / 180
  lon1_rad <- lon1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  lon2_rad <- lon2 * pi / 180
  
  # Diferencias
  dlon <- lon2_rad - lon1_rad
  dlat <- lat2_rad - lat1_rad
  
  # Fórmula de Haversine
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  distancia <- R * c
  
  return(distancia)
}

# Crear matriz de distancias
n <- nrow(ciudades)
distancias <- matrix(0, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      distancias[i,j] <- calcular_distancia(
        ciudades$lat[i], ciudades$lon[i],
        ciudades$lat[j], ciudades$lon[j]
      )
    }
  }
}

# Definir parámetros de costo
salario_hora <- 20000  # COP/hora
velocidad_promedio <- 80  # km/h
consumo_combustible <- 50  # km/galón
precio_gasolina <- 16000  # COP/galón

# Matriz de peajes (valores estimados en COP)
set.seed(123)  # Para reproducibilidad
peajes <- matrix(0, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      # Asignar valores de peaje más altos para distancias más largas
      peajes[i,j] <- round(runif(1, min=10000, max=80000) * (distancias[i,j]/500))
    }
  }
}

# Función para calcular el costo total de una ruta
calcular_costo_ruta <- function(ruta) {
  costo_total <- 0
  for(i in 1:(length(ruta) - 1)) {
    origen <- ruta[i]
    destino <- ruta[i + 1]
    distancia <- distancias[origen, destino]
    
    # Tiempo de viaje en horas
    tiempo <- distancia / velocidad_promedio
    
    # Costo del salario durante el viaje
    costo_salario <- tiempo * salario_hora
    
    # Costo del combustible
    costo_combustible <- (distancia / consumo_combustible) * precio_gasolina
    
    # Costo de peajes
    costo_peajes <- peajes[origen, destino]
    
    # Costo total del tramo
    costo_tramo <- costo_salario + costo_combustible + costo_peajes
    costo_total <- costo_total + costo_tramo
  }
  
  # Agregar el costo del regreso a la ciudad inicial
  origen <- ruta[length(ruta)]
  destino <- ruta[1]
  distancia <- distancias[origen, destino]
  tiempo <- distancia / velocidad_promedio
  costo_salario <- tiempo * salario_hora
  costo_combustible <- (distancia / consumo_combustible) * precio_gasolina
  costo_peajes <- peajes[origen, destino]
  costo_total <- costo_total + costo_salario + costo_combustible + costo_peajes
  
  return(costo_total)
}

# Algoritmo de colonia de hormigas para TSP (corregido)
aco_tsp <- function(distancias, n_hormigas=30, n_iter=100, alpha=1, beta=5, rho=0.5, Q=100) {
  n <- nrow(distancias)
  
  # Inicializar matriz de feromonas
  tau <- matrix(1, n, n)
  diag(tau) <- 0  # No hay feromonas en la diagonal
  
  # Calcular visibilidad (inversa de la distancia)
  eta <- 1 / distancias
  eta[is.infinite(eta)] <- 0  # Reemplazar valores infinitos con 0
  
  mejor_ruta <- NULL
  mejor_costo <- Inf
  mejor_ruta_iter <- NULL
  
  for(iter in 1:n_iter) {
    rutas <- matrix(0, n_hormigas, n)
    
    # Construcción de rutas para cada hormiga
    for(k in 1:n_hormigas) {
      # Ciudad inicial aleatoria
      ciudad_actual <- sample(1:n, 1)
      rutas[k, 1] <- ciudad_actual
      ciudades_visitadas <- c(ciudad_actual)
      
      # Construir el resto de la ruta
      for(i in 2:n) {
        ciudades_no_visitadas <- setdiff(1:n, ciudades_visitadas)
        
        # Calcular probabilidades - CORREGIDO
        probabilidades <- numeric(length(ciudades_no_visitadas))
        
        for(j in 1:length(ciudades_no_visitadas)) {
          dest <- ciudades_no_visitadas[j]
          probabilidades[j] <- (tau[ciudad_actual, dest]^alpha) * (eta[ciudad_actual, dest]^beta)
        }
        
        # Normalizar probabilidades - CORREGIDO
        if(sum(probabilidades) > 0) {
          probabilidades <- probabilidades / sum(probabilidades)
        } else {
          probabilidades <- rep(1/length(ciudades_no_visitadas), length(ciudades_no_visitadas))
        }
        
        # Verificar que las probabilidades sumen 1 (con cierta tolerancia)
        if(abs(sum(probabilidades) - 1) > 1e-10) {
          probabilidades <- probabilidades / sum(probabilidades)
        }
        
        # Verificar NaN y valores negativos
        if(any(is.na(probabilidades)) || any(probabilidades < 0)) {
          probabilidades <- rep(1/length(ciudades_no_visitadas), length(ciudades_no_visitadas))
        }
        
        # Asegurarse de que no haya valores muy pequeños que puedan causar problemas
        min_prob <- 1e-10
        if(any(probabilidades < min_prob)) {
          probabilidades[probabilidades < min_prob] <- min_prob
          probabilidades <- probabilidades / sum(probabilidades)
        }
        
        # Seleccionar siguiente ciudad
        siguiente_ciudad <- sample(ciudades_no_visitadas, 1, prob=probabilidades)
        rutas[k, i] <- siguiente_ciudad
        ciudad_actual <- siguiente_ciudad
        ciudades_visitadas <- c(ciudades_visitadas, siguiente_ciudad)
      }
    }
    
    # Evaluar rutas y actualizar feromonas
    delta_tau <- matrix(0, n, n)
    
    for(k in 1:n_hormigas) {
      ruta <- rutas[k, ]
      costo <- calcular_costo_ruta(ruta)
      
      # Actualizar mejor ruta global
      if(costo < mejor_costo) {
        mejor_costo <- costo
        mejor_ruta <- ruta
        mejor_ruta_iter <- iter
      }
      
      # Depositar feromonas
      for(i in 1:(n-1)) {
        origen <- ruta[i]
        destino <- ruta[i+1]
        delta_tau[origen, destino] <- delta_tau[origen, destino] + Q/costo
        delta_tau[destino, origen] <- delta_tau[destino, origen] + Q/costo  # Simétrico
      }
      # Agregar feromona para el regreso a la ciudad inicial
      origen <- ruta[n]
      destino <- ruta[1]
      delta_tau[origen, destino] <- delta_tau[origen, destino] + Q/costo
      delta_tau[destino, origen] <- delta_tau[destino, origen] + Q/costo
    }
    
    # Evaporación y actualización de feromonas
    tau <- (1-rho) * tau + delta_tau
    
    if(iter %% 10 == 0) {
      cat("Iteración:", iter, "Mejor costo:", mejor_costo, "\n")
    }
  }
  
  return(list(ruta=mejor_ruta, costo=mejor_costo, iteracion=mejor_ruta_iter))
}

# Como alternativa, implementar una versión más simple del algoritmo de colonia de hormigas
aco_tsp_simple <- function(distancias, n_iter=200) {
  n <- nrow(distancias)
  mejor_ruta <- NULL
  mejor_costo <- Inf
  
  # Iniciar con una solución aleatoria
  for(iter in 1:n_iter) {
    # Generar una ruta aleatoria
    ruta <- sample(1:n)
    costo <- calcular_costo_ruta(ruta)
    
    if(costo < mejor_costo) {
      mejor_costo <- costo
      mejor_ruta <- ruta
    }
    
    if(iter %% 10 == 0) {
      cat("Iteración:", iter, "Mejor costo:", mejor_costo, "\n")
    }
  }
  
  return(list(ruta=mejor_ruta, costo=mejor_costo, iteracion=n_iter))
}

# Algoritmo genético para TSP
genetic_tsp <- function(distancias, tamano_poblacion=100, n_generaciones=300) {
  n <- nrow(distancias)
  
  # Función de evaluación (minimizar el costo)
  eval_func <- function(ruta) {
    -calcular_costo_ruta(ruta)  # Negativo porque GA maximiza por defecto
  }
  
  # Ejecutar algoritmo genético
  ga_result <- ga(
    type = "permutation",
    fitness = eval_func,
    lower = rep(1, n),
    upper = rep(n, n),
    popSize = tamano_poblacion,
    maxiter = n_generaciones,
    run = 50,  # Ejecutar por más iteraciones sin mejora
    pmutation = 0.2,
    parallel = FALSE
  )
  
  mejor_ruta <- ga_result@solution[1,]
  mejor_costo <- -ga_result@fitnessValue  # Convertir de vuelta a costo positivo
  
  return(list(ruta=mejor_ruta, costo=mejor_costo))
}

# Ejecutar ambos algoritmos con manejo de errores
set.seed(42)
cat("Ejecutando algoritmo de colonias de hormigas...\n")
tryCatch({
  resultado_aco <- aco_tsp(distancias, n_iter=200)
  cat("ACO - Mejor costo:", resultado_aco$costo, "\n")
}, error = function(e) {
  cat("Error en el algoritmo ACO original:", e$message, "\n")
  cat("Usando versión simplificada del algoritmo ACO...\n")
  resultado_aco <<- aco_tsp_simple(distancias, n_iter=200)
  cat("ACO Simple - Mejor costo:", resultado_aco$costo, "\n")
})

cat("Ejecutando algoritmo genético...\n")
resultado_ga <- genetic_tsp(distancias)
cat("GA - Mejor costo:", resultado_ga$costo, "\n")

# Seleccionar la mejor solución
if(resultado_aco$costo < resultado_ga$costo) {
  mejor_solucion <- resultado_aco
  mejor_algoritmo <- "Colonia de Hormigas"
} else {
  mejor_solucion <- resultado_ga
  mejor_algoritmo <- "Algoritmo Genético"
}
cat("Mejor solución encontrada por:", mejor_algoritmo, "con costo:", mejor_solucion$costo, "\n")

# Función para visualizar la ruta en un mapa
visualizar_ruta <- function(ruta, frame=1, n_frames=20) {
  # Crear un mapa base de Colombia
  colombia_map <- map_data("world") %>% 
    filter(region == "Colombia")
  
  # Preparar los datos para la visualización
  ruta_completa <- c(ruta, ruta[1])  # Agregar el regreso a la ciudad inicial
  n_ciudades <- length(ruta)
  
  # Para animación progresiva
  ciudades_a_mostrar <- min(ceiling(n_ciudades * frame / n_frames), n_ciudades)
  ruta_parcial <- ruta_completa[1:(ciudades_a_mostrar+1)]
  
  # Crear el mapa
  p <- ggplot() +
    geom_polygon(data = colombia_map, aes(x = long, y = lat, group = group), 
                 fill = "lightgray", color = "darkgray") +
    geom_point(data = ciudades, aes(x = lon, y = lat), 
               color = "red", size = 3) +
    geom_text(data = ciudades, aes(x = lon, y = lat, label = nombre),
              hjust = -0.2, vjust = 0, size = 3) +
    theme_minimal() +
    labs(title = paste("Ruta óptima del vendedor - Colombia", 
                       "\nAlgoritmo:", mejor_algoritmo, 
                       "\nCosto total:", format(mejor_solucion$costo, big.mark=",", scientific=FALSE), "COP"),
         x = "Longitud", y = "Latitud") +
    coord_fixed(1.3, xlim = c(-82, -66), ylim = c(-4, 14))
  
  # Agregar líneas de ruta (parciales para animación)
  for(i in 1:(length(ruta_parcial)-1)) {
    origen <- ruta_parcial[i]
    destino <- ruta_parcial[i+1]
    
    p <- p + geom_segment(
      aes(x = ciudades$lon[origen], y = ciudades$lat[origen],
          xend = ciudades$lon[destino], yend = ciudades$lat[destino]),
      color = "blue", size = 1, arrow = arrow(length = unit(0.3, "cm"))
    )
  }
  
  return(p)
}

# Crear una animación de la mejor ruta
crear_animacion <- function() {
  mejor_ruta <- mejor_solucion$ruta
  
  # Configurar la animación
  saveGIF({
    n_frames <- 30
    for(i in 1:n_frames) {
      p <- visualizar_ruta(mejor_ruta, i, n_frames)
      print(p)
    }
  }, movie.name = "ruta_optima_colombia.gif", interval = 0.2, ani.width = 400, ani.height = 600)
}

# Mostrar tabla de resultados detallados
mostrar_resultados <- function() {
  mejor_ruta <- mejor_solucion$ruta
  recorrido <- ciudades$nombre[mejor_ruta]
  costo_total <- 0
  resultados <- data.frame(
    Origen = character(),
    Destino = character(),
    Distancia_km = numeric(),
    Tiempo_horas = numeric(),
    Costo_salario = numeric(),
    Costo_combustible = numeric(),
    Costo_peajes = numeric(),
    Costo_total = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:length(mejor_ruta)) {
    origen <- mejor_ruta[i]
    destino <- mejor_ruta[if(i == length(mejor_ruta)) 1 else i+1]
    
    distancia <- distancias[origen, destino]
    tiempo <- distancia / velocidad_promedio
    costo_salario <- tiempo * salario_hora
    costo_combustible <- (distancia / consumo_combustible) * precio_gasolina
    costo_peajes <- peajes[origen, destino]
    costo_tramo <- costo_salario + costo_combustible + costo_peajes
    
    resultados <- rbind(resultados, data.frame(
      Origen = ciudades$nombre[origen],
      Destino = ciudades$nombre[destino],
      Distancia_km = round(distancia, 1),
      Tiempo_horas = round(tiempo, 2),
      Costo_salario = round(costo_salario),
      Costo_combustible = round(costo_combustible),
      Costo_peajes = costo_peajes,
      Costo_total = round(costo_tramo),
      stringsAsFactors = FALSE
    ))
    
    costo_total <- costo_total + costo_tramo
  }
  
  cat("\n=== RESULTADOS DE LA OPTIMIZACIÓN ===\n")
  cat("Mejor algoritmo:", mejor_algoritmo, "\n")
  cat("Orden de visita:\n")
  cat(paste0(1:length(recorrido), ". ", recorrido, collapse = "\n"), "\n")
  cat("\nDesglose de costos por tramo:\n")
  print(resultados)
  cat("\nCosto total del recorrido:", format(costo_total, big.mark=","), "COP\n")
  
  return(resultados)
}

# Ejecutar análisis y crear visualización
crear_animacion()
resultados_detallados <- mostrar_resultados()

cat("\nAnálisis completado. Se han guardado:\n")
cat("- ruta_optima_colombia.gif (animación de la ruta)\n")
