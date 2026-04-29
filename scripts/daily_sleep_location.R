##############################################
# Deetction of the daily rest location
#load packages
library(data.table)
library(RPostgreSQL)
library(sp) 
library(lubridate)
library(maptools)
library(geosphere)
library(rgeos)
library(rgdal)
library(raster)
library(mgcv)
library(ggplot2)

# upload positions from database
positions <- positions[ up_validpos==TRUE & up_gamres < 75, ]
positions[, date := as.Date(up_timestamp_utc)]

# calculation of swimming distance
comp.dist <- function (easting, northing, depth, timestamp){
  diff.time <- c( NA, diff(as.numeric(timestamp), lag = 1))   # difference in time
  diff.east <- c(NA, diff(easting, lag = 1))                  # difference of easting
  diff.north <- c(NA, diff(northing, lag = 1))                # difference of northing
  diff.depth <- c(NA,diff(depth, lag=1))                      # difference in depth
  swim.dist2D <- sqrt(diff.east^2+diff.north^2)               # swimmed distance in 2D not including depth 
  swim.dist3D <- sqrt(swim.dist2D^2+diff.depth^2)             # swimmed distance in 3D not including depth 
  return(list(diff.depth = diff.depth, diff.time = diff.time, swim.dist2D = swim.dist2D, swim.dist3D = swim.dist3D))
}

setkey(positions, fi_fishid, up_timestamp_utc)
positions[,c("diff.depth", "diff.time", "swim.dist2D_e", "swim.dist3D_e") := comp.dist(up_egam, up_ngam, up_depth, up_timestamp_utc), by = fi_fishid]
positions[, ss_2d := swim.dist2D_e/diff.time]

# kalkulace distance to shore
positions_df <- as.data.frame(positions)

# Převedeme na sf body
positions_sf <- st_as_sf(positions_df, coords = c("up_easting", "up_northing"), crs = st_crs(shape.rimov))
# Výpočet vzdáleností (v metrech, pokud CRS je metrický)
#shoreline <- st_boundary(shape.rimov)
#vzdalenosti <- st_distance(positions_sf, shoreline)
#vzdalenosti_m <- as.numeric(vzdalenosti)

# description of individuals
# I would start working with these three 
# T449313_1 – week 24 -39 – map 53200 - stacionary behaviour 
# T449215_1 – week 24 – 38 – map 43400 - from week 24 - 27 migration through tha whole reservoir, after week 27 stacionary behaviour 
# T449317_1 – low horizontal movement just vertical, since 2017−07−19 – map 53600

# other individulas: short decription 
# T449202_1 – fish fully available from 2017−10−09 – similar behavior – map 42100
# T449213_1 - fish fully available from 2017−08−21 – higher movements – map 43200 
# T449310_1 - fish fully available from 2017−09−18 – similar behavior – map 52900
# T449208_1 – od 3.cervence


fish_data_raw <- positions[date > "2017-06-12" & date < "2017-11-30"]

# Přejmenování sloupců pro usnadnění práce
setnames(fish_data_raw, c("up_easting", "up_northing", "up_timestamp_utc"), c("x", "y", "time"))

#ggplot(fish_data_raw, aes(x = x, y = y))+geom_point()
# Výpočet rychlosti (pokud `ss_2d` není dostatečně spolehlivý nebo chcete přepočítat)
# Použijeme data.table syntaxi pro efektivitu
fish_data_prep <- fish_data_raw  %>%
  group_by(fi_fishid) %>% # Ujistěte se, že rychlost počítáte pro každou rybu zvlášť
  mutate(
    # Použijte diff.time přímo, pokud je spolehlivý. Jinak počítáme rozdíl časů.
    # dt = as.numeric(time - lag(time), units = "secs"),
    dx = x - lag(x),
    dy = y- lag(y),
    # speed = sqrt(dx^2 + dy^2) / ifelse(dt > 0, dt, NA) # Pokud chcete přepočítat rychlost
    speed = ss_2d # Použijeme existující `ss_2d` pro rychlost
  ) %>%
  ungroup()


library(RcppRoll)
# Vyhlazení rychlosti pomocí pohyblivého průměru (rollmean z balíčku zoo)
# To pomáhá odstranit šum a lépe identifikovat období klidu.
# Vyhlazení rychlosti pomocí pohyblivého průměru (rollmean z balíčku zoo)
# To pomáhá odstranit šum a lépe identifikovat období klidu.
window_size <- 5 # počet bodů pro průměr (můžete experimentovat)
fish_data_prep <- fish_data_prep %>%
  group_by(fi_fishid) %>%
  mutate(
    smoothed_speed = RcppRoll::roll_meanr(speed, n = window_size, fill = NA)
  ) %>%
  ungroup()

str(fish_data_prep )

# Vizualizace vyhlazené rychlosti pro první rybu
find_sleep_location_day <- function(data_subset, minPts_dbscan = 50, eps_val_dbscan = 20, sleep_speed_threshold = 0.01) {
  if (nrow(data_subset) < minPts_dbscan) {
    return(data.frame(x_sleep = NA, y_sleep = NA, cluster_size = 0, dbscan_eps = eps_val_dbscan, status = "Not enough points in subset"))
  }
  
  sleeping_candidates <- data_subset %>%
    filter(!is.na(smoothed_speed), smoothed_speed < sleep_speed_threshold, !is.na(x), !is.na(y))
  
  if (nrow(sleeping_candidates) < minPts_dbscan) {
    return(data.frame(x_sleep = NA, y_sleep = NA, cluster_size = 0, dbscan_eps = eps_val_dbscan, status = "No low-speed points for clustering"))
  }
  
  clustering_result <- tryCatch({
    dbscan::dbscan(sleeping_candidates[, c("x", "y")], eps = eps_val_dbscan, minPts = minPts_dbscan)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(clustering_result) || is.null(clustering_result$cluster)) {
    return(data.frame(x_sleep = NA, y_sleep = NA, cluster_size = 0, dbscan_eps = eps_val_dbscan, status = "DBSCAN failed"))
  }
  
  sleeping_candidates$cluster <- clustering_result$cluster
  
  valid_clusters <- sleeping_candidates %>%
    filter(cluster != 0) %>%
    group_by(cluster) %>%
    tally() %>%
    arrange(desc(n))
  
  if (nrow(valid_clusters) == 0) {
    return(data.frame(x_sleep = NA, y_sleep = NA, cluster_size = 0, dbscan_eps = eps_val_dbscan, status = "No valid sleep cluster found"))
  }
  
  largest_cluster_id <- valid_clusters$cluster[1]
  largest_cluster_size <- valid_clusters$n[1]
  
  sleep_coords <- sleeping_candidates %>%
    filter(cluster == largest_cluster_id) %>%
    summarise(x_sleep = mean(x), y_sleep = mean(y))
  
  return(data.frame(
    x_sleep = sleep_coords$x_sleep,
    y_sleep = sleep_coords$y_sleep,
    cluster_size = largest_cluster_size,
    dbscan_eps = eps_val_dbscan,
    status = "Detected"
  ))
}

# --- Aplikace funkce pro každý den a každou rybu pomocí nest() a map() ---
# 1. Vnoříme data podle fi_fishid a date
nested_data <- fish_data_prep %>%
  group_by(fi_fishid, date) %>%
  nest() # Toto vytvoří sloupec 'data', který obsahuje data.frame pro každou skupinu

# 2. Použijeme purrr::map_dfr k aplikaci funkce na každý vnořený data.frame
# map_dfr automaticky spojí vrácené data.framy dohromady
daily_sleep_locations <- nested_data %>%
  mutate(
    sleep_info = map(data, ~ find_sleep_location_day(., minPts_dbscan = 50, eps_val_dbscan = 20, sleep_speed_threshold = 0.015))
  ) %>%
  unnest(sleep_info) %>% # Rozbalí sloupec 'sleep_info' do jednotlivých sloupců
  ungroup() %>%
  filter(!is.na(x_sleep) & cluster_size > 0) # Odfiltrujte dny, kde nebyl spánek detekován

# vypocet zda ryba spala na danem miste i daalsi den 
# perimetr 50 m 
# Výpočet vzdálenosti k místu předešlé noci
daily_sleep_locations[, `:=`(
  x_prev = data.table::shift(x_sleep),
  y_prev = data.table::shift(y_sleep),
  date_prev = data.table::shift(date)
), by = fi_fishid]

# Eukleidovská vzdálenost k předešlému místu (jen pokud stejná ryba a noc předchozí)
daily_sleep_locations[, dist_prev := ifelse(date == date_prev + 1,
  sqrt((x_sleep - x_prev)^2 + (y_sleep - y_prev)^2),
  NA_real_
)]

# Nový sloupec, zda spali v perimetru 50 m
daily_sleep_locations[, same_place_50m := dist_prev <= 30]

# Výsledek
daily_sleep_locations <- daily_sleep_locations[, .(fi_fishid, date, x_sleep, y_sleep, dist_prev, same_place_50m)]

