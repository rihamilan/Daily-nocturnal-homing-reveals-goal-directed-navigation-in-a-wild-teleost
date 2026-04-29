##############################################
# VISUALIZATION of gam model results to weekly plots
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


con <-  dbConnect(drv = PostgreSQL(), dbname ="teridb", host="172.21.3.20", user= "teriuser", password = "t3r1us3r!")
# getting positions for two selected pikeperch
#select species
species <- "pikeperch"

# make selection
select.fish.qu <- 
  paste("SELECT fi_fishid
FROM teri.fish
WHERE
fi_species = '", species,"' ;", sep = "")

# getting tag numbers
tagid.list <- data.table(dbGetQuery(con,select.fish.qu ))
tagid.list <- tagid.list[,fi_fishid]

# query for positions
select.position.qu <- paste("WITH tmp as (SELECT  * FROM teri.positions 
                            where 
                            up_fishvalid AND
                            up_timestamp_utc BETWEEN '2017-04-27 00:00:00' AND '2017-11-20 23:59:59' AND 
                            tu_tagmode = 'fish' AND
                            fi_fishid IN ('", paste(tagid.list, collapse = "','"),"')
                            )

                            SELECT up_timestamp_utc, b.lt_lake, up_easting, up_northing,up_egam, up_ngam, ST_X(ST_Transform(up_geom, 4326)) as lon,ST_Y(ST_Transform(up_geom, 4326)) as lat, up_depth, up_validpos, fi_fishid, up_gamres, up_fishvalid,  up_bottom_depth,
                                   ST_LineLocatePoint(lcl_centerline, up_geom_mid) * ST_length(lcl_centerline) distfromdam 
                                FROM tmp b INNER JOIN teri.lakecenterline a ON b.lt_lake = a.lt_lake ", sep = "")


# upload positions from database
positions <- data.table(dbGetQuery(con, select.position.qu))
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

# getting shape of the Rimov
library(raster)
slope_rast <- raster("~/Teri/shp_files/rimov_slope_raster/tmp_slope_raster.tif")
plot(slope_rast)
slope_rast_lowres <- aggregate(slope_rast, fact=10, fun=mean)

shape.rimov  <- st_read("~/Teri/shp_files/rimov_469/Rimov_pol_469m_UTM33.shp") 
shape_points_rim <- fortify(shape.rimov)[,1:2]
names(shape_points_rim) <- c("x", "y")
shape_points_rim <- shape_points_rim[1:(length(shape_points_rim$x)),]

raster.sub <- crop(slope_rast, extent(shape_points_rim))
raster.sub <- mask(raster.sub, shape.rimov)
raster.sub_lowres <- aggregate(raster.sub, fact=2, fun=mean)

# Lake Most
#Lat 48.836447
#Lon 14.482155

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
#ggplot(fish_data_prep , aes(x = time, y = smoothed_speed )) +
# geom_line(alpha = 0.7) +
# labs(title = "Vyhlazená rychlost ryby v čase", y = "Rychlost (m/s)") +
# theme_minimal() +
# scale_y_log10() # Logaritmická škála pro lepší zobrazení nízkých rychlostí


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

# Zobrazení výsledků
print(daily_sleep_locations)

# Vizualizace denních spánkových lokací
ggplot(daily_sleep_locations, aes(x = date, y = x_sleep, color = fi_fishid)) +
  geom_line() +
  geom_point() +
  labs(title = "Denní X souřadnice místa spánku", x = "Datum", y = "X souřadnice (Easting)") +
  theme_minimal()+
  
  ggplot(daily_sleep_locations, aes(x = date, y = y_sleep, color = fi_fishid)) +
  geom_line() +
  geom_point() +
  labs(title = "Denní Y souřadnice místa spánku", x = "Datum", y = "Y souřadnice (Northing)") +
  theme_minimal()

ggplot()+geom_sf(data= shape.rimov)+
  geom_point(data = daily_sleep_locations, aes(x = x_sleep, y = y_sleep, color = date), size = 3) +
  geom_path(data = daily_sleep_locations, aes(x = x_sleep, y = y_sleep, color = date, group = fi_fishid), arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  labs(title = "Trajektorie místa spánku v prostoru (po dnech)", x = "X souřadnice (Easting)", y = "Y souřadnice (Northing)") +
  theme_minimal() +
  scale_color_gradientn(colors = viridisLite::viridis(10))#+
#coord_sf(xlim = c(462000, 462700), ylim =c(5410300, 5411000))

library(plotly)
library(viridisLite)

p <- ggplot() +
  geom_sf(data = shape.rimov, fill = "grey95", color = "black") +
  geom_point(data = daily_sleep_locations, 
             aes(x = x_sleep, y = y_sleep, color = date, text = paste("Fish ID:", fi_fishid)), 
             size = 3) +
  geom_path(data = daily_sleep_locations, 
            aes(x = x_sleep, y = y_sleep, color = date, group = fi_fishid, 
                text = paste("Fish ID:", fi_fishid, "<br>Date:", date)), 
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  labs(
    title = "Trajektorie místa spánku v prostoru (po dnech)", 
    x = "X souřadnice (Easting)", 
    y = "Y souřadnice (Northing)"
  ) +
  theme_minimal() +
  scale_color_gradientn(colors = viridis(10))

# převod na interaktivní plotly graf
ggplotly(p, tooltip = "text")

#write_csv(daily_sleep_locations, "~/Teri/pikeperch_navigation/data/daily_sleep_locations.csv")

# dalsi rozsireni 
#Chceš, abych ti připravil skript, který:
# shromáždí výsledky těchto denních clusterů,
# agreguje je podle jedince,
# vyhodnotí stabilitu a opakovatelnost míst odpočinku?

daily_sleep_locations <- data.table(read_csv("~/Teri/pikeperch_navigation/data/daily_sleep_locations.csv"))
daily_sleep_locations <- daily_sleep_locations[date < "2017-09-24"]
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

daily_sleep_locations[, count := .N, by = .(fi_fishid, x_sleep, y_sleep)]

# Odstraníme NA
daily_sleep_na <- na.omit(daily_sleep_locations[, .(fi_fishid, x_sleep, y_sleep, date)])

# Výstupní tabulka pro všechny ryby
daily_sleep_clusters <- rbindlist(lapply(split(daily_sleep_na, by = "fi_fishid"), function(dt) {
  # Výběr souřadnic
  coords <- dt[, .(x_sleep, y_sleep)]
  
  # DBSCAN (např. eps = 20 m, minPts = 3)
  clust <- dbscan(coords, eps = 20, minPts = 3)
  
  # Výstup s přiřazeným clusterem
  dt[, cluster_id := clust$cluster]
  return(dt)
}))

sleep_cluster_unique<- unique(daily_sleep_clusters[,.(fi_fishid, x_sleep, y_sleep, cluster_id)])

ggplot(sleep_cluster_unique[fi_fishid == "T449215_1"], aes(x = x_sleep, y = y_sleep, col =as.factor(cluster_id)))+geom_point()


daily_sleep_clusters[cluster_id != 0, mean_x_sleep := mean(x_sleep), by = .(fi_fishid, cluster_id)]
daily_sleep_clusters[cluster_id == 0, mean_x_sleep := x_sleep]
daily_sleep_clusters[cluster_id != 0, mean_y_sleep := mean(y_sleep), by = .(fi_fishid, cluster_id)]
daily_sleep_clusters[cluster_id == 0, mean_y_sleep := y_sleep]

daily_sleep_clusters[, count := length(unique(date)), by = .(fi_fishid,  cluster_id)]
daily_sleep_clusters[cluster_id == 0, count_true := 1]
daily_sleep_clusters[cluster_id != 0, count_true := count]

daily_sleep_clusters <- merge(daily_sleep_clusters, daily_sleep_locations[,.(fi_fishid, date, same_place_50m)], by = c("fi_fishid", "date"))
daily_sleep_clusters[is.na(same_place_50m ), same_place_50m  := F]


daily_sleep_clusters[fi_fishid == "T449215_1" & same_place_50m == F, ]

# cluster colour, 0 
# Aplikace podmínky
daily_sleep_clusters[
  cluster_id != 0,  # mimo cluster 0
  same_place_50m := if (any(same_place_50m)) TRUE else same_place_50m,
  by = .(fi_fishid, cluster_id)
]

# 
daily_sleep_clusters[, same_place_50m := factor(same_place_50m, levels = c(TRUE, FALSE))]

library(ggspatial)
ggplot(daily_sleep_clusters[fi_fishid == "T449313_1"]) +
  geom_path(aes(x = mean_x_sleep, y = mean_y_sleep), alpha = 0.1) +
  geom_point(aes(x = mean_x_sleep, y = mean_y_sleep, size = count_true, color = same_place_50m), alpha = 0.7) +
  #scale_color_viridis_c(option = "D", name = "Den (číselně)") +
  #scale_size_continuous(name = "Počet nocí") +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(~fi_fishid) +
  labs(title = "Spánková místa: četnost a časová distribuce", x = "Easting [m]", y = "Northing [m]")
library(ggrepel)
daily_sleep_clusters_5 <- daily_sleep_clusters[fi_fishid == "T449215_1"]
daily_sleep_clusters_5[cluster_id == 1 , cluster_1 := T]
daily_sleep_clusters_5[cluster_id != 1 , cluster_1 := F]
daily_sleep_clusters_5[, cluster_1 := factor(cluster_1)]

daily_sleep_clusters_3 <- daily_sleep_clusters[fi_fishid == "T449313_1"]
daily_sleep_clusters_3[cluster_id == 1 , cluster_1 := T]
daily_sleep_clusters_3[cluster_id != 1 , cluster_1 := F]
daily_sleep_clusters_3[, cluster_1 := factor(cluster_1)]


plot_m_5 <- ggplot() +
  geom_sf(data = shape.rimov, fill = NA, color = alpha("dodgerblue3", 0.5), linewidth = 0.5) +
  geom_point(
    data = daily_sleep_clusters_5,
    aes(x = mean_x_sleep, y = mean_y_sleep, size = count_true,
        fill = same_place_50m, shape = cluster_1),
    color = "black", stroke = 0.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # <- tvary s výplní
  scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "#E41A1C")) +
  coord_sf()  +
  annotate("text", x = min(shape.rimov$geometry[[1]][[1]][,1]) + 200, 
           y = max(shape.rimov$geometry[[1]][[1]][,2]) - 200,
           label = "#2", size = 4, fontface = "bold") +
  geom_label(data = daily_sleep_clusters_5[same_place_50m == TRUE],
             aes(x = mean_x_sleep, y = mean_y_sleep, label = count_true),
             size = 3.5, fill = "white", color = "red", label.size = 0, fontface = "bold",
             label.padding = unit(0, "lines"), nudge_y = 0, nudge_x = 500) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  labs(title = "b)", x = "Easting [m]", y = "Northing [m]") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
    #plot.title = element_blank()
  )

plot_m_5
plot_m_3 <- ggplot() +
  geom_sf(data = shape.rimov, fill = NA, color = alpha("dodgerblue3", 0.5), linewidth = 0.5) +
  geom_point(
    data = daily_sleep_clusters_3,
    aes(x = mean_x_sleep, y = mean_y_sleep, size = count_true,
        fill = same_place_50m, shape = cluster_1),
    color = "black", stroke = 0.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # <- tvary s výplní
  scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "#E41A1C")) +
  coord_sf() +
  annotate("text", x = min(shape.rimov$geometry[[1]][[1]][,1]) + 200, 
           y = max(shape.rimov$geometry[[1]][[1]][,2]) - 200,
           label = "#1", size = 4, fontface = "bold") +
  geom_label(data = daily_sleep_clusters_3[same_place_50m == TRUE],
             aes(x = mean_x_sleep, y = mean_y_sleep, label = count_true),
             size = 3.5, fill = "white", color = "red", label.size = 0, fontface = "bold",
             label.padding = unit(0, "lines"), nudge_y = 0, nudge_x = 500) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  labs(title = "a)", x = "Easting [m]", y = "Northing [m]") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
    #plot.title = element_blank()
  )
plot_m_3

plot_t_5 <- ggplot() +
  geom_path(data = daily_sleep_clusters_5, aes(x = date, y = mean_y_sleep), alpha = 0.1) +
  geom_point(data = daily_sleep_clusters_5, 
             aes(x = date, y = mean_y_sleep, fill = same_place_50m, shape = cluster_1),
             color = "black", stroke = 0.5, size = 2) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # <- tvary s výplní
  scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "#E41A1C")) +
  #scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) substr(format(x, "%B"), 1, 3)) +
  ylab("Northing")+
  xlab("Time")+
  ggtitle("d)")+
  annotate("text",
           x = min(daily_sleep_clusters_5$date) + 5,
           y = min(daily_sleep_clusters_5$mean_y_sleep) + 50,
           label = "#2", size = 4, fontface = "bold") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", colour = "grey70", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "none",
    #plot.title = element_blank(),
    #axis.title.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12)
  )
plot_t_5
plot_t_3 <- ggplot() +
  geom_path(data = daily_sleep_clusters_3, aes(x = date, y = mean_y_sleep), alpha = 0.1) +
  geom_point(data = daily_sleep_clusters_3, 
             aes(x = date, y = mean_y_sleep, fill = same_place_50m, shape = cluster_1),
             color = "black", stroke = 0.5, size = 2) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # <- tvary s výplní
  scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "#E41A1C")) +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) substr(format(x, "%B"), 1, 3)) +
  ylab("Northing")+
  xlab("Time")+
  ggtitle("c)")+
  annotate("text",
           x = min(daily_sleep_clusters_3$date) + 5,
           y = min(daily_sleep_clusters_3$mean_y_sleep) + 50,
           label = "#1", size = 4, fontface = "bold") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", colour = "grey70", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "none",
   # plot.title = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12)
  )
plot_t_3

  
library(patchwork)
# Tři ukázkové grafy

# Definice rozložení přes area:
design <- "
AABBCC
AADDCC
"

# Složení:
final_plot <- (plot_m_3 + plot_t_3 + plot_m_5 + plot_t_5+ plot_layout(design = design))
final_plot 
ggsave("~/Teri/pikeperch_navigation/output/resting_place_summary.png", plot = final_plot, width = 8, height = 6, dpi = 300, units = "in")
