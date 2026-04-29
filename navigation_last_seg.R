library(dbscan)
library(gdistance)
library(patchwork)
library(scales)
library(sf)
library(suncalc)
library(zoo)
library(geosphere)
library(gdistance)

# getting shape of the Rimov
shape.rimov  <- st_read("~/Teri/shp_files/rimov_469/Rimov_pol_469m_UTM33.shp") 
fish_length <- data.table(fishid = c("T449215_1", "T449313_1"), fishid_new =  c("#2", "#1"), tl = c(0.49, 0.495))
# load data
nav_sum_sub <- data.table(read_csv("~/Teri/pikeperch_navigation/data/nav_sum_sub.csv"))
nav_last_seg  <- data.table(read_csv("~/Teri/pikeperch_navigation/data/navigation_segments.csv"))
# filtrovani
nav_seg_summary <- nav_sum_sub[fishid %in% c("T449215_1","T449313_1") ] 
last_seg_sub <- nav_last_seg[navigation_event_id %in% nav_seg_summary$navigation_event_id &
                               transition_phase %in% c("before_first", "after_last") ]
# coverage 
traj_cov <- unique(nav_seg_summary[,.(navigation_event_id, observation_coverage)])
last_seg_sub <- merge(last_seg_sub, traj_cov, by = c("navigation_event_id"))
# Definuj datum/daty
dates <- seq(as.Date("2017-06-01"), as.Date("2017-12-07"), by = "day")

# Definuj polohu Rimova
lat <- 48.85
lon <- 14.49

# Spočítej časy
sun_data <- data.table(getSunlightTimes(date = dates, lat = lat, lon = lon, keep = c("sunrise", "sunset")))

# filter out corrupted segments
#Before_first – 68, 82, 455, 8, 16, 69, 280
#After_last  -4, 26, 29
bad_bf <- c(68, 82, 455, 8,  16, 23, 69, 460, 280, 374, 233, 362, 368, 298,196, 142, 138, 131, 130, 141, 137)
bad_al <- c(4, 20,26, 374, 213, 29, 230, 203, 267)
whole_traj <- c(10, 17, 18, 249, 374)

# Odstranění problematických podle typu a ID
last_seg_sub <- last_seg_sub[!(transition_phase == "before_first" & navigation_event_id %in% bad_bf)]
last_seg_sub <- last_seg_sub[!(transition_phase == "after_last" & navigation_event_id %in% bad_al)]
last_seg_sub <- last_seg_sub[!(navigation_event_id %in% whole_traj)]

# bottom depth below zero assign to 0 
last_seg_sub[bottom_depth <0, bottom_depth := 0]
# rest place cluster
# Předpoklad: last_seg_sub obsahuje sloupce fishid, x_sleep, y_sleep
# Odstraníme NA
sleep_coords <- na.omit(last_seg_sub[, .(fishid, x_sleep, y_sleep)])

# Výstupní tabulka pro všechny ryby
sleep_clusters <- rbindlist(lapply(split(sleep_coords, by = "fishid"), function(dt) {
  # Výběr souřadnic
  coords <- dt[, .(x_sleep, y_sleep)]
  
  # DBSCAN (např. eps = 20 m, minPts = 3)
  clust <- dbscan(coords, eps = 20, minPts = 3)
  
  # Výstup s přiřazeným clusterem
  dt[, cluster_id := clust$cluster]
  return(dt)
}))

sleep_cluster_unique<- unique(sleep_clusters[,.(fishid, x_sleep, y_sleep, cluster_id)])

# now calulate distance to shore ####
last_seg_sub[, n_points := length(timestamp), by = .(fishid, navigation_event_id)]
range(last_seg_sub$n_points)
last_seg_sub_df <- data.frame(last_seg_sub)
# Převedeme na sf body
last_seg_sub_sf <- st_as_sf(last_seg_sub_df, coords = c("easting", "northing"), crs = st_crs(shape.rimov))
# Výpočet vzdáleností (v metrech, pokud CRS je metrický)
shoreline <- st_boundary(shape.rimov)
dist2shore_df <- st_distance(last_seg_sub_sf, shoreline)
dist2shore_m <- as.numeric(dist2shore_df)
# assing distance to shore 
last_seg_sub[, dist_to_shore := dist2shore_m]

n_segments <- last_seg_sub[,.(n_seg = length(unique(navigation_event_id))), by =.(fishid, transition_phase)]
n_segments

# rozdeleni na segmenty  > 50 m od brehu
last_seg_sub[dist_to_shore > 50 ,shore_50 := F]
last_seg_sub[dist_to_shore <= 50 ,shore_50 := T]

# interpolace hloubky
last_seg_sub[, depth_interp := na.approx(depth, x = timestamp, na.rm = FALSE), by = .(fishid, navigation_event_id)]

# vypocet distance to bottom 
last_seg_sub[,d2b := bottom_depth-depth_interp]

# Determination of shore and offhore segments 
# Výchozí pořadí
setorder(last_seg_sub, navigation_event_id, timestamp)

# # 1. Převod shore logiky na číselnou podobu
# last_seg_sub[, shore_50_fix := fifelse(shore_50 == TRUE, 1L, 0L)]
# 
# # 2. Identifikace přechodů mezi shore/offshore
# last_seg_sub[, change := as.logical(shore_50_fix) != data.table::shift(as.logical(shore_50_fix), fill = as.logical(shore_50_fix)[1]), 
#              by = .(fishid, navigation_event_id)]
# 
# # 3. Identifikace bloků (souvislé úseky stejných typů)
# last_seg_sub[, shore_block := cumsum(change), by = .(fishid, navigation_event_id)]
# 
# # 4. Spočítání délky bloků
# block_sizes <- last_seg_sub[, .N, by = .(fishid, navigation_event_id, shore_block)]
# 
# # 5. Najdi krátké bloky (obě kategorie: shore i offshore)
# short_blocks <- block_sizes[N < 6, .(fishid, navigation_event_id, shore_block)]
# 
# # 6. Nastav tyto krátké bloky jako NA (k vyplnění)
# last_seg_sub[short_blocks, on = .(fishid, navigation_event_id, shore_block), shore_50_fix := NA]
# 
# # 7. Vyplň chybějící hodnoty dopředu i dozadu
# last_seg_sub[, shore_50_merged := nafill(shore_50_fix, type = "locf")]
# last_seg_sub[, shore_50_merged := nafill(shore_50_merged, type = "nocb")]
# 
# # 8. Vytvoř nové segmentové ID podle výsledného merged stavu
# last_seg_sub[, shore_group_final := rleid(shore_50_merged), by = .(fishid, transition_phase, navigation_event_id)]
# last_seg_sub[, shore_merged_group := rleid(shore_50_fix), by = .(fishid, transition_phase, navigation_event_id)]

# 
# Ujisti se o pořadí uvnitř fáze
setorder(last_seg_sub, fishid, navigation_event_id, transition_phase, timestamp)

# 1) Převod na 0/1
last_seg_sub[, shore01 := as.integer(shore_50)]  # 1 = nearshore, 0 = offshore

# 2) RLE po FÁZÍCH (důležité!) a délky běhů
last_seg_sub[, run_id := data.table::rleid(shore01),
             by = .(fishid, navigation_event_id, transition_phase)]

runs <- last_seg_sub[, .(n = .N, val = first(shore01)),
                     by = .(fishid, navigation_event_id, transition_phase, run_id)]

# 3) Info o sousedech pro každou fázi
runs[, `:=`(
  prev_val = data.table::shift(val),
  next_val = data.table::shift(val, type = "lead"),
  prev_n   = data.table::shift(n),
  next_n   = data.table::shift(n, type = "lead")
), by = .(fishid, navigation_event_id, transition_phase)]

# 4) Urči cílovou hodnotu pro krátké běhy (<6): delší ze sousedů (fallback na dostupného)
runs[, replace_with :=
       fifelse(n >= 6, val,
               fifelse(!is.na(prev_n) & (is.na(next_n) | prev_n >= next_n), prev_val, next_val)
       )]

# 5) Mapuj zpět do bodů
last_seg_sub <- runs[last_seg_sub,
                     on = .(fishid, navigation_event_id, transition_phase, run_id),
                     nomatch = 0L][
                       , shore01_clean := replace_with][]

# 6) Pokud byl běh už dost dlouhý, ponech původní
last_seg_sub[n >= 6, shore01_clean := val]

# 7) Finální sloupec nearshore/offshore po vyhlazení + nové skupiny
last_seg_sub[, shore_50_merged := as.integer(shore01_clean)]
last_seg_sub[, shore_group_final := data.table::rleid(shore_50_merged),
             by = .(fishid, navigation_event_id, transition_phase)]

# 9. Vzdálenost od předchozího bodu
last_seg_sub[, distance_from_previous := sqrt((easting - data.table::shift(easting))^2 + 
                                                (northing - data.table::shift(northing))^2), 
             by = .(fishid, transition_phase, navigation_event_id)]

# 10. Merge s odpočinkovým clusterem
last_seg_sub <- merge(last_seg_sub, sleep_cluster_unique, by = c("fishid", "x_sleep", "y_sleep"), all.x = TRUE)

# rychlost 
# Časový rozdíl (v sekundách)
last_seg_sub[, time_diff_sec := as.numeric(difftime(timestamp, data.table::shift(timestamp), units = "secs")), by = .(fishid, transition_phase, navigation_event_id)]

# Rychlost (v m/s)
last_seg_sub[, speed_m_per_s := distance_from_previous / time_diff_sec]

# add sleeping cluster 
sleep_clust_nav_id <- unique(last_seg_sub[,.(navigation_event_id,cluster_id)])
#fwrite(sleep_clust_nav_id[cluster_id == 1,], "~/Teri/pikeperch_navigation/data/sleep_clust_nav_id5_3.csv")
#fwrite(last_seg_sub, "~/Teri/pikeperch_navigation/data/firstlast_seg_final.csv")
# calculation of segments summary for each shore/ offshore segmants
segments_fragm <- last_seg_sub[, .(
  start_time = min(timestamp),
  end_time = max(timestamp),
  n_points = .N,
  from_x = first(easting),
  from_y = first(northing),
  to_x = last(easting),
  to_y = last(northing),
  dist_total = sum(distance_from_previous, na.rm = TRUE),
  d_direct = sqrt((last(easting) - first(easting))^2 + (last(northing) - first(northing))^2),
  shore_50 = first(shore_50),
  sd_d2b = sd(d2b, na.rm = TRUE),
  mean_d2b = mean(d2b, na.rm = TRUE),
  x_sleep = mean(x_sleep,na.rm = T),
  y_sleep = mean(y_sleep,na.rm = T),
  n_points = length(timestamp) 
), by = .(fishid,transition_phase, navigation_event_id, shore_group_final)][
  , `:=`(
    tortuosity = fifelse(d_direct > 0, dist_total / d_direct, NA_real_),
    azimuth_deg = (atan2(to_x - from_x, to_y - from_y) * 180 / pi) %% 360
  )
]

segments_fragm[, duration_seg := as.numeric(difftime(end_time,start_time, units = "sec"))]
segments_fragm[, cov_pos := n_points/(duration_seg/15)]
segments_fragm[, length_ev := length(shore_group_final), by = .(fishid, navigation_event_id, transition_phase)]

segments_fragm[navigation_event_id == 91]

# distance to bottom 
ggplot(last_seg_sub[  transition_phase == "after_last"], aes(x = dist_home , y = depth_interp, group = navigation_event_id, col = shore_50))+
  geom_point()+
  facet_wrap(~fishid)

# speed
ggplot(last_seg_sub[  transition_phase == "after_last"], aes(x = dist_home , y = speed_m_per_s,  col = shore_50))+
  geom_smooth()+
  facet_wrap(~fishid)+
  xlim(0,500)+
  ylim(0, 2)
  

# detour nova verze s pocatecnim bearinggem ####

# --- Funkce: výpočet průměrného bearingu po zadané vzdálenosti ---
mean_bearing <- function(coords_utm, dist_limit = 100, crs_utm = 32633) {
  if (nrow(coords_utm) < 2) return(NA_real_)
  
  sf_line <- st_linestring(coords_utm)
  sf_geom <- st_sfc(sf_line, crs = crs_utm)
  
  # ⬇️ Oprava: výběr jen prvních dvou sloupců (lon, lat)
  coords_lonlat <- st_coordinates(st_transform(sf_geom, crs = 4326))[, 1:2]
  
  cumdist <- c(0, cumsum(sqrt(rowSums(diff(coords_utm)^2))))
  idx <- which(cumdist <= dist_limit)
  if (length(idx) < 2) return(NA_real_)
  
  bearings <- numeric(length(idx) - 1)
  for (i in seq_along(bearings)) {
    bearings[i] <- geosphere::bearing(coords_lonlat[idx[i], ], coords_lonlat[idx[i+1], ])
  }
  
  bearings_circ <- circular::circular(bearings, units = "degrees", template = "geographics")
  as.numeric(circular::mean.circular(bearings_circ))
}

# --- Převod bodů na sf ---
last_seg_sf <- st_as_sf(last_seg_sub[!is.na(x_sleep)], coords = c("easting", "northing"), crs = 32633, remove = FALSE)

last_seg_sf_ordered <- last_seg_sf %>%
  arrange(fishid, navigation_event_id, transition_phase, timestamp)

last_seg_dt <- as.data.table(last_seg_sub)
valid_navs <- last_seg_sub[, .N, by = .(fishid, navigation_event_id, transition_phase)][N >= 2]
last_seg_sub_filtered <- last_seg_sf_ordered %>% 
  semi_join(valid_navs, by = c("fishid", "navigation_event_id", "transition_phase"))

traj_list <- last_seg_sub_filtered %>%
  group_by(fishid, navigation_event_id, transition_phase) %>%
  group_split() %>%
  lapply(function(df) {
    if (nrow(df) < 2) return(NULL)
    coords <- st_coordinates(df)
    st_linestring(coords)
  })

valid_idx <- which(!sapply(traj_list, is.null))
traj_list <- traj_list[valid_idx]

group_ids <- last_seg_sub_filtered %>%
  distinct(fishid, navigation_event_id, transition_phase) %>%
  slice(valid_idx)

traj_sf <- st_sf(group_ids, geometry = st_sfc(traj_list, crs = 32633))
traj_sf <- traj_sf %>%
  filter(lengths(geometry) > 0)

# --- Raster vody a cost surface ---
#buffer_dist <- 30
#shape.rimov_buffered <- st_buffer(shape.rimov, dist = buffer_dist)
#bbox_rimov <- st_bbox(shape.rimov_buffered)
# bufferovaný polygon

# water_rast <- raster(
#   xmn = bbox_rimov["xmin"],
#   ymn = bbox_rimov["ymin"],
#   xmx = bbox_rimov["xmax"],
#   ymx = bbox_rimov["ymax"],
#   resolution = 5,
#   crs = st_crs(shape.rimov)$wkt
# )
# values(water_rast) <- 1
# water_rast <- mask(water_rast, as(shape.rimov_buffered, "Spatial"))
# 
# tr <- transition(water_rast, transitionFunction = function(x) 1, directions = 16)
# tr <- geoCorrection(tr, type = "c")

shape.rimov_buffered <- readRDS("~/Teri/pikeperch_navigation/data/shape.rimov_buffered.rds")
water_rast <- readRDS("~/Teri/pikeperch_navigation/data/rimov_water_rast.rds")
tr <- readRDS("~/Teri/pikeperch_navigation/data/cost_surface_tr.rds")

# --- Výpočet trajektorií a bearingů ---
nav_dt_for_join <- as.data.table(last_seg_dt)

start_end_pts <- nav_dt_for_join[, 
                                 .(
                                   start_easting = easting[which.min(timestamp)],
                                   start_northing = northing[which.min(timestamp)],
                                   end_easting = easting[which.max(timestamp)],
                                   end_northing = northing[which.max(timestamp)]
                                 ),
                                 by = .(fishid, navigation_event_id, transition_phase)
]

results_detour <- traj_sf %>%
  left_join(start_end_pts, by = c("fishid", "navigation_event_id", "transition_phase")) %>%
  rowwise() %>%
  mutate(
    traj_length = st_length(geometry),
    start_pt_sf = st_sfc(st_point(c(start_easting, start_northing)), crs = 32633),
    end_pt_sf = st_sfc(st_point(c(end_easting, end_northing)), crs = 32633)
  ) %>%
  mutate(
    # Empirický průměrný bearing (počáteční úsek trajektorie)
    empirical_bearing = {
      traj_coords <- st_coordinates(geometry)[, 1:2]
      mean_bearing(traj_coords, dist_limit = 100)
    },
    
    # Cost distance
    cost_distance = {
      if (!is.na(start_pt_sf) && !is.na(end_pt_sf) && !is.null(tr)) {
        sp_start <- as(start_pt_sf, "Spatial")
        sp_end <- as(end_pt_sf, "Spatial")
        dist <- tryCatch(
          costDistance(tr, sp_start, sp_end),
          error = function(e) {
            message(paste("Error calculating costDistance for fishid:", fishid, "nav_id:", navigation_event_id, "-", e$message))
            NA_real_
          }
        )
        as.numeric(dist)
      } else {
        NA_real_
      }
    },
    
    # Least-cost path bearing (průměrný počáteční bearing LCP)
    lcp_bearing = {
      if (!is.na(start_pt_sf) && !is.na(end_pt_sf) && !is.null(tr)) {
        sp_start <- as(start_pt_sf, "Spatial")
        sp_end <- as(end_pt_sf, "Spatial")
        path <- tryCatch(
          shortestPath(tr, sp_start, sp_end, output = "SpatialLines"),
          error = function(e) return(NULL)
        )
        if (!is.null(path)) {
          coords <- coordinates(path)[[1]][[1]]
          mean_bearing(coords, dist_limit = 100)
        } else {
          NA_real_
        }
      } else {
        NA_real_
      }
    }
  ) %>%
  mutate(
    bearing_diff = ifelse(!is.na(empirical_bearing) & !is.na(lcp_bearing),
                          abs(empirical_bearing - lcp_bearing),
                          NA_real_),
    detour_ratio = ifelse(cost_distance > 0 & !is.infinite(cost_distance),
                          as.numeric(traj_length) / cost_distance,
                          NA_real_)
  ) %>%
  ungroup() %>%  # ← přesunuto před select
  dplyr::select(fishid, navigation_event_id, transition_phase, traj_length, cost_distance,
                detour_ratio, empirical_bearing, lcp_bearing, bearing_diff)

results_detour <- results_detour %>%
  mutate(
    bearing_diff_adj = {
      diff_raw <- abs(empirical_bearing - lcp_bearing) %% 360
      ifelse(diff_raw > 180, 360 - diff_raw, diff_raw)
    }
  )
# merge with sleep cluster
results_detour <- merge(results_detour, sleep_clust_nav_id, by = c("navigation_event_id"))
# as numeric 
results_detour$traj_length <- as.numeric(results_detour$traj_length)
results_detour$cost_distance <- as.numeric(results_detour$cost_distance)

results_detourcl1 <- results_detour %>% filter(cluster_id == 1 & transition_phase == "after_last")
sleep_spot_cluster1 <- last_seg_sub[cluster_id == 1,.(x = mean(x_sleep), y= mean(y_sleep))]
#2. Kolik metrů strávila ve volné vodě vs. u břehu
#Např. pomocí shore_50 masky, vypočítáš délku trajektorie ve volné vodě před návratem:

last_seg_sub[, segment_length := shift(distance_from_previous, type = "lead"), by = .(fishid, navigation_event_id, transition_phase)]

shore_segment_summary <- last_seg_sub[!is.na(segment_length), .(
  total_length = sum(as.numeric(segment_length), na.rm = TRUE),
  shore_length = sum(as.numeric(ifelse(shore_50_merged == 1, segment_length, 0)), na.rm = TRUE),
  openwater_length = sum(as.numeric(ifelse(shore_50_merged == 0, segment_length, 0)), na.rm = TRUE)
), by = .(fishid, navigation_event_id, transition_phase)][, shore_ratio := shore_length / total_length]

# Přidej podíl podél břehu
shore_segment_summary[, shore_ratio := shore_length / total_length]

# merge s predchozim datasetem
results_shore <- merge(results_detour, shore_segment_summary, by = c("fishid", "navigation_event_id", "transition_phase"), all.x = TRUE)

# summary for whole last segments
segments_sum <- last_seg_sub[, .(
  start_time = min(timestamp),
  end_time = max(timestamp),
  n_points = .N,
  from_x = first(easting),
  from_y = first(northing),
  to_x = last(easting),
  to_y = last(northing),
  from_light = first(light_lux),
  to_light = last(light_lux),
  dist_total = sum(distance_from_previous, na.rm = TRUE),
  d_direct = sqrt((last(easting) - first(easting))^2 + (last(northing) - first(northing))^2),
  sd_d2b = sd(d2b, na.rm = TRUE),
  mean_d2b = mean(d2b, na.rm = TRUE),
  x_sleep = mean(x_sleep,na.rm = T),
  y_sleep = mean(y_sleep,na.rm = T)
), by = .(fishid,transition_phase, navigation_event_id)][
  , `:=`(
    tortuosity = fifelse(d_direct > 0, dist_total / d_direct, NA_real_),
    azimuth_deg = (atan2(to_x - from_x, to_y - from_y) * 180 / pi) %% 360,
    duration_min = as.numeric(difftime(end_time,start_time, units = "mins" ))
  )
]
names(segments_sum)

segments_sum <- as_tibble(segments_sum)

results_detour_dt <- data.table(results_detour)
results_detour_dt[,.(n_seg = length(unique(navigation_event_id))), by =.(fishid, transition_phase)]
    
###################
# TABLE #####
library(dplyr)
# Přidání duration do results_shore
segments_sum <- as_tibble(segments_sum)

results_joined <- results_shore %>%
  dplyr::left_join(
    segments_sum %>% 
      dplyr::select(fishid, transition_phase, navigation_event_id, duration_min),
    by = c("fishid", "transition_phase", "navigation_event_id")
  ) %>%
  dplyr::mutate(SI = cost_distance / traj_length)

# Přidání duration do results_shore
results_joined <- results_shore %>%
  left_join(
    as.data.frame(segments_sum)[, c("fishid", "transition_phase", "navigation_event_id", "duration_min")],
    by = c("fishid", "transition_phase", "navigation_event_id")
  ) %>%
  mutate(SI = cost_distance / traj_length)

# Shrnutí po rybě a transition fázi
summary_table <- results_joined %>%
  filter(cluster_id == 1)%>%
  st_drop_geometry() %>%
  group_by(fishid, transition_phase) %>%
  summarise(
    N_events = n_distinct(navigation_event_id),
    N_segments = n(),
    
    traj_length_mean = mean(traj_length),
    traj_length_sd = sd(traj_length),
    traj_length_range = paste0(round(min(traj_length), 1), " - ", round(max(traj_length), 1)),
    
    cost_distance_mean = mean(cost_distance),
    cost_distance_sd = sd(cost_distance),
    cost_distance_range = paste0(round(min(cost_distance), 1), " - ", round(max(cost_distance), 1)),
    
    duration_mean = mean(duration_min, na.rm = TRUE),
    duration_sd = sd(duration_min, na.rm = TRUE),
    duration_range = paste0(round(min(duration_min, na.rm = TRUE), 1), " - ", round(max(duration_min, na.rm = TRUE), 1)),
    
    detour_ratio_mean = mean(detour_ratio),
    detour_ratio_sd = sd(detour_ratio),
    detour_ratio_range = paste0(round(min(detour_ratio), 2), " - ", round(max(detour_ratio), 2)),
    
    shore_ratio_mean = mean(shore_ratio),
    shore_ratio_sd = sd(shore_ratio),
    shore_ratio_range = paste0(round(min(shore_ratio), 2), " - ", round(max(shore_ratio), 2)),
    
    SI_mean = mean(SI),
    SI_sd = sd(SI),
    SI_range = paste0(round(min(SI), 2), " - ", round(max(SI), 2))
  ) %>%
  ungroup()

# Výstup
print(summary_table)
library(flextable)
library(officer)
# Zaokrouhlení čísel pro lepší čitelnost
summary_fmt <- summary_table %>%
  mutate(across(where(is.numeric), ~round(., 2)))

# Flextable
ft <- flextable(summary_fmt)
ft <- autofit(ft)
ft <- theme_booktabs(ft)
ft <- set_caption(ft, "Summary of navigation metrics per fish and transition phase")

# Export do Wordu
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
#print(doc, target = "~/Teri/pikeperch_navigation/output/navigation_summary_table_segments.docx")

metrics <- c("traj_length", "cost_distance", "duration_min", "detour_ratio", "shore_ratio", "SI")
fish_ids <- unique(results_joined$fishid)

# Funkce na t-test mezi transition_phase u jedné ryby
run_tests_per_fish <- function(fish_id) {
  df <- results_joined %>% filter(fishid == fish_id)
  
  res <- lapply(metrics, function(metric) {
    x <- df %>% filter(transition_phase == "before_first") %>% pull(!!sym(metric))
    y <- df %>% filter(transition_phase == "after_last") %>% pull(!!sym(metric))
    
    test <- t.test(x, y)
    data.frame(
      fishid = fish_id,
      metric = metric,
      p_value = round(test$p.value, 4),
      mean_before = round(mean(x, na.rm = TRUE), 2),
      mean_after = round(mean(y, na.rm = TRUE), 2),
      direction = ifelse(mean(y, na.rm = TRUE) > mean(x, na.rm = TRUE), "↑ after", "↓ after")
    )
  })
  do.call(rbind, res)
}

# Aplikace pro všechny ryby
test_results <- lapply(fish_ids, run_tests_per_fish) %>% bind_rows()

# Výstup tabulky
print(test_results)

# vizualization prezentation ####
results_detourcl1_all <- results_joined %>% filter(cluster_id == 1)
results_detourcl1_all <- data.table(results_detourcl1_all )
results_detourcl1_all[transition_phase == "before_first",transition_phase_id := "odchod"]
results_detourcl1_all[transition_phase == "after_last",transition_phase_id := "navrat"]
results_detourcl1_all[, transition_phase_id := factor(transition_phase_id, levels = c("odchod","navrat"))]
results_detourcl1_all[fishid == "T449313_1", fishid_new := "#1"]
results_detourcl1_all[fishid == "T449215_1", fishid_new := "#2"]

pr_gg1 <- ggplot(results_detourcl1_all, aes(
  x = transition_phase_id,
  y = traj_length,
  col = fishid_new
)) +
  geom_violin(
    scale = "width",
    trim = FALSE, size =1.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.3,
    alpha = 0.5
  ) +
  ## >>> PRÁVĚ TADY PŘIDÁME group = fishid_new
  stat_summary(
    aes(group = fishid_new),
    fun = "mean",
    geom = "point",
    shape = 23,
    size = 3,
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  labs(x = "", y = "delka trajektorie (m)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    text = element_text(size = 16)
  )

pr_gg2 <-  ggplot(results_detourcl1_all, aes(
  x = transition_phase_id,
  y = cost_distance,
  col = fishid_new
)) +
  geom_violin(
    scale = "width",
    trim = FALSE, size =1.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.3,
    alpha = 0.5
  ) +
  ## >>> PRÁVĚ TADY PŘIDÁME group = fishid_new
  stat_summary(
    aes(group = fishid_new),
    fun = "mean",
    geom = "point",
    shape = 23,
    size = 3,
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  labs(x = "", y = "prima vzdalenost (m)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    text = element_text(size = 16)
  )

pr_gg3 <-  ggplot(results_detourcl1_all, aes(
  x = transition_phase_id,
  y = duration_min,
  col = fishid_new
)) +
  geom_hline(yintercept = c(60, 120, 180, 240))+
  geom_violin(
    scale = "width",
    trim = FALSE, size =1.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.3,
    alpha = 0.5
  ) +
  ## 👉 jen popisky pro hodnoty z geom_hline
  scale_y_continuous(
    breaks = c(60, 120, 180, 240),
    labels = c("60", "120", "180", "240")
  ) +
  ## >>> PRÁVĚ TADY PŘIDÁME group = fishid_new
  stat_summary(
    aes(group = fishid_new),
    fun = "mean",
    geom = "point",
    shape = 23,
    size = 3,
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  labs(x = "", y = "trvani odchodu/navratu (min)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    text = element_text(size = 16)
  )

pr_gg4 <-  ggplot(results_detourcl1_all, aes(
  x = transition_phase_id,
  y = shore_ratio,
  col = fishid_new
)) +
  #geom_hline(yintercept = c(60, 120, 180, 240))+
  geom_violin(
    scale = "width",
    trim = FALSE, size =1.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.3,
    alpha = 0.5
  ) +
  ## >>> PRÁVĚ TADY PŘIDÁME group = fishid_new
  stat_summary(
    aes(group = fishid_new),
    fun = "mean",
    geom = "point",
    shape = 23,
    size = 3,
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  labs(x = "", y = "pomer breh/volna voda") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    text = element_text(size = 16)
  )

(pr_gg1+pr_gg2)/(pr_gg3+pr_gg4)
###############
#vizualizace ######

# Fig 2 - descriptive parameters ####
# vizualizace casoveho rozlozeni segmentu ####
# Převod data ve sun_data
# Převod data ve sun_data
sun_data[, date := as.Date(date)]

# Přidání data do segments_sum
segments_sum <- data.table(segments_sum)
segments_sum[, date := as.Date(start_time)]

# Merge sunrise/sunset časy
segments_sum_light <- merge(segments_sum, sun_data[, .(date, sunrise, sunset)], by = "date", all.x = TRUE)

# Funkce na převod času na decimální hodiny
to_hour <- function(x) hour(x) + minute(x)/60 + second(x)/3600

# Výpočet hodin
segments_sum_light[, `:=`(
  start_hour = to_hour(start_time),
  end_hour = to_hour(end_time),
  sunrise_hour = to_hour(sunrise),
  sunset_hour = to_hour(sunset)
)]

# Výběr before_first a after_last segmentů
before <- segments_sum_light[transition_phase == "before_first"]
after <- segments_sum_light[transition_phase == "after_last"]

# Spojení do jedné tabulky
cycles <- merge(
  before[, .(fishid, navigation_event_id, T0 = start_hour, T1 = end_hour, sunrise_hour, sunset_hour)],
  after[, .(fishid, navigation_event_id, T2 = start_hour, T3 = end_hour)],
  by = c("fishid", "navigation_event_id")
)
# Přidej logický sloupec: návrat před východem slunce
cycles[, t3_before_sunrise := T3 < sunrise_hour]

# Spočítej pravděpodobnost celkově
# Přidej indikátor návratu před východem slunce
cycles[, t3_before_sunrise := T3 < sunrise_hour]

# Spočítej pravděpodobnost pro každou rybu
prob_by_fish <- cycles[
  , .(prob_t3_before_sunrise = mean(t3_before_sunrise, na.rm = TRUE)),
  by = fishid
]

# Ujistěte se, že navigation_event_id je faktor
# Předpokládáme, že cycles je data.table
cycles[, navigation_event_id := as.factor(navigation_event_id)]

# Připrav data
plot_data_5 <- cycles %>%
  filter(fishid == "T449215_1") %>%
  mutate(navigation_event_id = factor(navigation_event_id),
         nav_id_num = as.numeric(navigation_event_id))  # Převod jen pro vykreslení obdélníků

plot_data_3 <- cycles %>%
  filter(fishid == "T449313_1") %>%
  mutate(navigation_event_id = factor(navigation_event_id),
         nav_id_num = as.numeric(navigation_event_id))  # Převod jen pro vykreslení obdélníků


# Data pro první rybu
plot_data_5_long <- plot_data_5 %>%
  pivot_longer(cols = c(T0, T1, T2, T3), names_to = "phase", values_to = "time") %>%
  mutate(phase = recode(phase,
                        "T0" = "leave RP",
                        "T1" = "start foraging",
                        "T2" = "stop forage_start to RP",
                        "T3" = "arrival to RP"
  ))

# Data pro druhou rybu
plot_data_3_long <- plot_data_3 %>%
  pivot_longer(cols = c(T0, T1, T2, T3), names_to = "phase", values_to = "time") %>%
  mutate(phase = recode(phase,
                        "T0" = "leave RP",
                        "T1" = "start foraging",
                        "T2" = "stop forage_start to RP",
                        "T3" = "arrival to RP"
  ))

# Barvy fází
phase_colors <- c(
  "leave RP" = "blue",
  "start foraging" = "green",
  "stop forage_start to RP" = "orange",
  "arrival to RP" = "red"
)

# První graf
plot2 <- ggplot(plot_data_5_long) +
  geom_rect(aes(xmin = 0, xmax = sunrise_hour, ymin = nav_id_num - 0.5, ymax = nav_id_num + 0.5),
            fill = "grey80", alpha = 0.4) +
  geom_rect(aes(xmin = sunset_hour, xmax = 24, ymin = nav_id_num - 0.5, ymax = nav_id_num + 0.5),
            fill = "grey80", alpha = 0.4) +
  geom_point(aes(x = time, y = as.factor(navigation_event_id), color = phase), size = 2) +
  scale_color_manual(values = phase_colors) +
  scale_x_continuous(breaks = seq(0, 24, 2), limits = c(0, 24)) +
  annotate("text", x = 22, y = 10, label = "#2", size = 12)+
  labs(y = "Forage event", title = "b)") +
  theme_minimal() +
  scale_color_manual(
    breaks = c("leave RP", "start foraging", "stop forage_start to RP", "arrival to RP"),
    values = c(
      "leave RP" = "blue",
      "start foraging" = "green",
      "stop forage_start to RP" = "orange",
      "arrival to RP" = "red"
    ),
    labels = c(
      "leave RP" = "leave RP",
      "start foraging" = "start foraging",
      "stop forage_start to RP" = "stop foraging",
      "arrival to RP" = "arrival to RP"
    )
  )+
  theme_minimal() +
  theme(text = element_text(size = 15), 
        #plot.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

plot2
# Druhý graf
plot1 <- ggplot(plot_data_3_long) +
  geom_rect(aes(xmin = 0, xmax = sunrise_hour, ymin = nav_id_num - 0.5, ymax = nav_id_num + 0.5),
            fill = "grey80", alpha = 0.4) +
  geom_rect(aes(xmin = sunset_hour, xmax = 24, ymin = nav_id_num - 0.5, ymax = nav_id_num + 0.5),
            fill = "grey80", alpha = 0.4) +
  geom_point(aes(x = time, y = navigation_event_id, color = phase), size = 2) +
  scale_color_manual(values = phase_colors) +
  scale_x_continuous(breaks = seq(0, 24, 2), limits = c(0, 24)) +
  annotate("text", x = 22, y = 20, label = "#1", size = 12)+
  scale_color_manual(
    breaks = c("leave RP", "start foraging", "stop forage_start to RP", "arrival to RP"),
    values = c(
      "leave RP" = "blue",
      "start foraging" = "green",
      "stop forage_start to RP" = "orange",
      "arrival to RP" = "red"
    ),
    labels = c(
      "leave RP" = "leave RP",
      "start foraging" = "start foraging",
      "stop forage_start to RP" = "stop foraging",
      "arrival to RP" = "arrival to RP"
    )
  )+
  labs(x = "Day time (hours)", y = "Forage event", title = "a)") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        #plot.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),  # Skrýt čísla na ose x
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) # Skrýt popisek osy x
#axis.title.y = element_text(vjust = -2, hjust = 0)) # Skrýt popisek osy x

  # theme_minimal() +
  # theme(text = element_text(size = 15), 
  #   #plot.title = element_blank(),
  #       panel.grid.major.x = element_blank(),
  #       panel.grid.major.y = element_line(linewidth = 0.1),
  #       panel.grid.minor.x = element_blank(),
  #       legend.position = "bottom",
  #       legend.title = element_blank(),
  #       axis.title.y = element_blank(),
  #       axis.text.y = element_blank())+

# Kombinace a společná legenda
combined_plot_timing_phases <- plot1 / plot2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Výsledek
combined_plot_timing_phases
# Uložení grafu do PNG ve vysokém rozlišení
# ggsave("odchod_navrat.png", plot = combined_plot_timing_phases, width = 16, height = 12, units = "cm", dpi = 600)

# duration of each phase 
# Vyfiltrujeme potřebné sloupce pro každý typ fáze
before <- segments_sum[transition_phase == "before_first", .(navigation_event_id, fishid, start_time)]
after  <- segments_sum[transition_phase == "after_last",  .(navigation_event_id, end_time)]

# Sloučíme dohromady podle navigation_event_id
event_durations <- merge(before, after, by = "navigation_event_id")

# Spočítáme trvání v minutách (nebo libovolné jednotce)
event_durations[, duration_mins := as.numeric(difftime(end_time, start_time, units = "mins"))]

# Výstup
event_durations[]
event_durations[fishid == "T449313_1" , fishid_new := "#1"]
event_durations[fishid == "T449215_1" , fishid_new := "#2"]
event_durations[, fishid_new := factor(fishid_new, levels =c("#1", "#2"))]
#event_durations[ transition_phase == "after_last" , tra_phase_new := "Return"]
#event_durations[ transition_phase == "before_first" , tra_phase_new := "Before forage"]

duration_gg <- ggplot(event_durations, aes(x = fishid_new, y = duration_mins/60, fill = fishid_new)) +
  geom_violin(trim = FALSE, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, aes(color = fishid)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels = function(x) gsub("_(T\\d+)", "\n\\1", x)) +  # hezčí osy
  labs(
    x = NULL,
    y = "Duration of forage events (hours)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 14)
  )
duration_gg 

# load max distance from the sleeping spot
results_length_traj <- data.table(read_csv("~/Teri/pikeperch_navigation/data/results_length_traj.csv"))
length_traj_sub <- results_length_traj[fishid %in% c("T449215_1", "T449313_1")] 
length_traj_sub[fishid == "T449215_1", new_fishid := "#2"]
length_traj_sub[fishid == "T449313_1", new_fishid := "#1"]
length_traj_sub[, new_fishid := factor(new_fishid, levels = c("#1", "#2"))]

cost_distance_gg <- ggplot(length_traj_sub, aes(x = new_fishid, y = max_dist_home, fill = fishid)) +
  # Violin plot bez ořezu a s černým outline
  geom_violin(trim = FALSE, alpha = 0.3, color = "black") +
  
  # Body s jitterem, zabarvené podle fishid
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, aes(color = fishid)) +
  
  # Medián jako bílý bod (shape 23 = čtverec s výplní)
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  
  # Přímé přiřazení barev (můžeš upravit podle potřeby)
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  
  labs(
    x = NULL,
    y = "Maximal distance from RP [m]"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )

# sunset and sunrise and navigation phases  ####


# Příprava dat
# Předpokládáme, že after_last obsahuje sloupce:
# - navigation_event_id
# - transition_phase (before_first, foraging, after_last, arrival to RP)
# - light_lux
# - fishid


# Parametry polohy
lat <- 48.85
lon <- 14.49
last_seg_sub[transition_phase == "after_last"]
# 1. Výběr klíčových časových bodů pro každou navigační událost
light_summary <- last_seg_sub[cluster_id == 1, .(
  departure_time = if (.N > 0 && any(transition_phase == "before_first")) 
    min(timestamp[transition_phase == "before_first"], na.rm = TRUE) else NA_POSIXct_,
  
  foraging_start_time = if (.N > 0 && any(transition_phase == "before_first")) 
    max(timestamp[transition_phase == "before_first"], na.rm = TRUE) else NA_POSIXct_,
  
  foraging_end_time = if (.N > 0 && any(transition_phase == "after_last")) 
    min(timestamp[transition_phase == "after_last"], na.rm = TRUE) else NA_POSIXct_,
  
  arrival_time = if (.N > 0 && any(transition_phase == "after_last")) 
    max(timestamp[transition_phase == "after_last"], na.rm = TRUE) else NA_POSIXct_
), by = .(fishid, navigation_event_id)]


# 2. Výpočet časů sunrise/sunset pro každou událost (podle arrival_time)
# 1. Získání dat pro východ a západ zvlášť podle každé fáze
light_summary[, `:=`(
  date_departure = as.Date(departure_time),
  date_foraging_start = as.Date(foraging_start_time),
  date_foraging_end = as.Date(foraging_end_time),
  date_arrival = as.Date(arrival_time)
)]

# 2. Výpočet časů slunce
light_summary[, `:=`(
  sunset_departure = getSunlightTimes(date_departure, lat, lon, keep = "sunset")$sunset,
  sunset_foraging_start = getSunlightTimes(date_foraging_start, lat, lon, keep = "sunset")$sunset,
  sunrise_foraging_end = getSunlightTimes(date_foraging_end, lat, lon, keep = "sunrise")$sunrise,
  sunrise_arrival = getSunlightTimes(date_arrival, lat, lon, keep = "sunrise")$sunrise
)]

# 3. Určení, jestli fáze byla před/po sluncem
light_summary[, `:=`(
  departure_after_sunset = departure_time > sunset_departure,
  foraging_start_after_sunset = foraging_start_time > sunset_foraging_start,
  foraging_end_before_sunrise = foraging_end_time < sunrise_foraging_end,
  arrival_before_sunrise = arrival_time < sunrise_arrival
)]


# Zjisti, kolik záznamů máš pro T449313_1 v arrival_before_sunrise
light_summary[fishid == "T449313_1", .(
  n_total = .N,
  n_valid = sum(!is.na(arrival_before_sunrise)),
  n_true = sum(arrival_before_sunrise, na.rm = TRUE)
)]


# Výpočet trvání výpravy v minutách
light_summary[, departure_sunset := as.numeric(difftime(departure_time,sunset_foraging_start,  units = "mins"))]
light_summary[, arrival_sunrise := as.numeric(difftime(arrival_time,sunrise_arrival,  units = "mins"))]
light_summary[, duration_event := as.numeric(difftime(arrival_time,departure_time,  units = "mins"))]
light_summary[, duration_departure := as.numeric(difftime(foraging_start_time,departure_time,  units = "mins"))]
light_summary[, duration_arrival := as.numeric(difftime(arrival_time ,foraging_end_time,  units = "mins"))]
light_summary[, stop_forage := as.numeric(difftime(sunrise_arrival ,foraging_end_time,  units = "mins"))]

mean_light <- light_summary[,.(mean_time = mean(departure_sunset, na.rm = T), sd_time = sd(departure_sunset, na.rm = T)), by = .(departure_after_sunset, fishid)]
mean_light <- light_summary[,.(mean_time = mean(stop_forage, na.rm = T), sd_time = sd(stop_forage, na.rm = T)), by = .( fishid)]



light_phase_probabilities <- light_summary[
  !is.na(departure_after_sunset),
  .(
    p_departure_after_sunset = mean(departure_after_sunset, na.rm = TRUE),
    p_foraging_start_after_sunset = mean(foraging_start_after_sunset, na.rm = TRUE),
    p_foraging_end_before_sunrise = mean(foraging_end_before_sunrise, na.rm = TRUE),
    p_arrival_before_sunrise = mean(arrival_before_sunrise, na.rm = TRUE)
  ),
  by = fishid
]

# Připrav data do long formátu
light_phase_long <- melt(
  light_phase_probabilities,
  id.vars = "fishid",
  variable.name = "phase",
  value.name = "probability"
)

# Připrav data do long formátu a převedeme zpět na data.table
light_phase_long <- melt(
  light_phase_probabilities,
  id.vars = "fishid",
  variable.name = "phase",
  value.name = "probability"
)
setDT(light_phase_long)  # Převod na data.table

# Převeď názvy fází na srozumitelnější štítky
light_phase_long[, phase := factor(phase, levels = c(
  "p_departure_after_sunset",
  "p_foraging_start_after_sunset",
  "p_foraging_end_before_sunrise",
  "p_arrival_before_sunrise"
), labels = c(
  "Departure\n(after sunset)",
  "Start foraging\n(after sunset)",
  "End foraging\n(before sunrise)",
  "Arrival\n(before sunrise)"
))]
# Vykresli stacked bar plot
library(scales)

light_phase_long <- light_phase_long %>%
  mutate(
    fishid_new = case_when(
      fishid == "T449215_1" ~ "#2",
      fishid == "T449313_1" ~ "#1",
      TRUE ~ fishid  # fallback
    )
  )
light_phase_long<- light_phase_long %>%
  mutate(fishid_new = factor(fishid_new, levels = c("#1", "#2")))

timing_gg <- ggplot(light_phase_long, aes(x = phase, y = probability, fill = fishid_new)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.6, width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Trip phase and timing",
    y = "Proportion of events",
    title = "c)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text( hjust = 0.5, size = 10),
    #axis.text.y = element_text(size = 14),
   # axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
   # legend.text = element_text(size = 14),
    #plot.title = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

# fig2 composition ####

timing_gg_full <-  ((plot1 / plot2) | timing_gg) +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom")
timing_gg_full
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/timing.png",
  plot = timing_gg_full,
  width = 12,      # šířka v palcích
  height = 8,     # výška v palcích
  dpi = 300
)

# Převod do long formátu
# Převeď do long formátu – jednoduše:
# Zajisti, že výchozí tabulka je data.table
light_summary <- as.data.table(light_summary)

# Melt pomocí data.table verze (NE z reshape2 nebo base R)
light_long <- data.table(melt(
  light_summary,
  measure.vars = c("departure_sunset", "arrival_sunrise"),
  variable.name = "event_type",
  value.name = "minutes"
))

# Kontrola (volitelná)
print(class(light_long))  # mělo by vrátit: "data.table" "data.frame"

# Teď tohle už bezpečně projde:
light_long[, event_type := factor(
  event_type,
  levels = c("departure_sunset", "arrival_sunrise"),
  labels = c("Departure", "Arrival")
)]

# Filtrování ne-NA hodnot
light_long <- light_long[!is.na(minutes)]
light_long <-  light_long%>%
  mutate(
    fishid_new = case_when(
      fishid == "T449215_1" ~ "#2",
      fishid == "T449313_1" ~ "#1",
      TRUE ~ fishid  # fallback
    )
  )

light_long <-  light_long %>%
  mutate(fishid_new = factor(fishid_new, levels = c("#1", "#2")))

ggplot(light_long, aes(x = event_type, y = minutes, fill = fishid_new)) +
  geom_violin(trim = FALSE, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, aes(color = fishid)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  facet_wrap(~fishid_new)+
  coord_cartesian(y = c(-150, 250))+
  scale_x_discrete(labels = function(x) gsub("_(T\\d+)", "\n\\1", x)) +  # hezčí osy
  labs(
    x = NULL,
    y = "Minutes from sunset / to sunrise"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13)
  )

# duration of event vs length ####
# Převedení wide -> long
light_long_traj <- data.table(melt(
  data = light_summary,
  id.vars = c("fishid", "navigation_event_id"),
  measure.vars = c("duration_departure", "duration_arrival"),
  variable.name = "transition_phase",
  value.name = "duration"
))
# Přejmenování faktorových úrovní transition_phase
light_long_traj[transition_phase == "duration_departure", transition_phase := "before_first"]
light_long_traj[transition_phase == "duration_arrival", transition_phase := "after_last"]

# Spojení s results_shore
# Pokud results_shore je sf objekt, konvertujeme data.frame část na data.table
results_shore_dt <- as.data.table(results_shore)

# Merge
light_long_traj <- merge(light_long_traj, results_shore_dt, 
                         by = c("fishid", "navigation_event_id", "transition_phase"), 
                         all.x = TRUE)
light_long_traj  <-  light_long_traj %>%
  mutate(
    fishid_new = case_when(
      fishid == "T449215_1" ~ "#2",
      fishid == "T449313_1" ~ "#1",
      TRUE ~ fishid  # fallback
    )
  )


light_long_traj[fishid_new == "#1" & transition_phase == "after_last" & duration > 150]
# Předpoklad: máš tabulku `slopes` s hodnotami `fishid_new`, `transition_phase`, `speed_BL_per_sec`


# Zkopíruj si data
dt <- copy(light_long_traj)

# Odstraň NA a převeď duration na sekundy
dt <- dt[!is.na(duration) & !is.na(traj_length)]
dt[, duration_sec := duration * 60]

# Výpočet slope (tedy rychlosti v m/s, protože y = vzdálenost, x = čas)
slopes <- dt[, {
  fit <- lm(traj_length ~ duration_sec)
  list(
    slope = coef(fit)[2],              # rychlost (m/s)
    intercept = coef(fit)[1],
    slope_se = summary(fit)$coefficients[2, 2], # Standardní chyba sklonu
    r2 = summary(fit)$r.squared
  )
}, by = .(fishid_new, transition_phase)]

# Připojení délky ryb
slopes <- merge(slopes, fish_length, by = "fishid_new")

# Rychlost v m/s už je přímo slope
slopes[, avg_speed_ms := slope]
slopes[, avg_speed_kmh := avg_speed_ms * 3.6]
slopes[, speed_BL_per_sec := avg_speed_ms / tl]
slopes[, speed_BL_se := slope_se / tl]

# Výpis
slopes[]
# Hodinové značky po 60 minutách — např. do 5 hodin = 300 minut
vlines <- seq(60, 300, by = 60)

# === GRAF: #5 ===
p5_length_dur <- ggplot(light_long_traj[fishid_new == "#2" & !is.na(duration)],
                   aes(x = duration , y = traj_length, col = transition_phase)) +
  geom_point(alpha = 0.8, shape = 21, stroke = 0.3, fill = NA, size = 2)+
  geom_smooth(method = "lm") +
  geom_vline(xintercept = vlines, linetype = "dashed", color = "lightgrey", linewidth = 0.3)+
  scale_color_manual(
    values = c("before_first" = "#00BFC4", "after_last" = "#F8766D"),
    labels = c("before_first" = "Before foraging", "after_last" = "Return"),
    name = NULL  # nebo třeba name = "Phase"
  ) +
  scale_fill_manual(
    values = c("before_first" = "#00BFC4", "after_last" = "#F8766D"),
    labels = c("before_first" = "Before foraging", "after_last" = "Return"),
    name = NULL
  )+
  lims(y = c(0,max(light_long_traj$traj_length,na.rm = T)), x = c(0,max(light_long_traj$duration,na.rm = T)))+
  annotate("text", x = 50, y = 4000,
           label = "Speed: 0.77 BL/s\nR² = 0.97",
           color = "#F8766D", hjust = 0, size = 5) +
  annotate("text", x = 100, y = 2000,
           label = "Speed: 0.80 BL/s\nR² = 0.91",
           color = "#00BFC4", hjust = 0, size = 5) +
  labs(x = "Duration (minutes)", y = "Trajectory length (m)", title = "b) #2") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    text = element_text(size = 15),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# === GRAF: #3 ===
p3_length_dur <- ggplot(light_long_traj[fishid_new == "#1" & !is.na(duration)],
                   aes(x = duration , y = traj_length, col = transition_phase)) +
  geom_point(alpha = 0.8, shape = 21, stroke = 0.3, fill = NA, size = 2)+
  geom_vline(xintercept = vlines, linetype = "dashed", color = "lightgrey", linewidth = 0.3)+
  geom_smooth(method = "lm") +
  scale_color_manual(
    values = c("before_first" = "#00BFC4", "after_last" = "#F8766D"),
    labels = c("before_first" = "Outbound", "after_last" = "Return"),
    name = NULL  # nebo třeba name = "Phase"
  ) +
  scale_fill_manual(
    values = c("before_first" = "#00BFC4", "after_last" = "#F8766D"),
    labels = c("before_first" = "Outbound", "after_last" = "Return"),
    name = NULL
  )+
  lims(y = c(0,max(light_long_traj$traj_length,na.rm = T)), x = c(0,max(light_long_traj$duration,na.rm = T)))+
  annotate("text", x = 50, y = 4000,
           label = "Speed: 0.92 BL/s\nR² = 0.92",
           color = "#F8766D", hjust = 0, size = 5) +
  annotate("text", x = 100, y = 2000,
           label = "Speed: 0.75 BL/s\nR² = 0.90",
           color = "#00BFC4", hjust = 0, size = 5) +
  labs(x = "Duration (minutes)", y = "Trajectory length (m)", title = "a) #1") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid = element_blank(),
    legend.position = c(0.7, 0.1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 15)
    
  )
p3_length_dur
# Spojení obou plotů
length_dur_gg <-  p3_length_dur + p5_length_dur 
length_dur_gg 
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/length_dur_gg.png",
  plot = length_dur_gg,
  width = 8,      # šířka v palcích
  height = 6,     # výška v palcích
  dpi = 300
)


# last two segments ####

# Předpoklad: shore_segments je data.table
# 1. Vyber fáze
before_first <- last_seg_sub[transition_phase == "before_first"]
after_last <- last_seg_sub[transition_phase == "after_last"]

# 2. První dva segmenty ve before_first
first_segs <- before_first[
  , {
    segs <- sort(unique(shore_group_final))
    list(shore_group_final = head(segs, 2))
  },
  by = .(fishid, navigation_event_id)
]

# 3. Poslední dva segmenty v after_last
last_segs <- after_last[
  , {
    segs <- sort(unique(shore_group_final))
    list(shore_group_final = tail(segs, 2))
  },
  by = .(fishid, navigation_event_id)
]

######
# Přepočítej rle po fázi
last_seg_sub[, seg_id_phase := data.table::rleid(shore_50_merged),
             by = .(fishid, navigation_event_id, transition_phase)]

# Spočti délky běhů po fázi
seg_len_phase <- last_seg_sub[, .N,
                              by = .(fishid, navigation_event_id, transition_phase, seg_id_phase)]

# Zachovej jen dostatečně dlouhé běhy (>= 6 bodů)
long_segs <- seg_len_phase[N >= 6]

# OUTBOUND (before_first): vezmi 1. a 2. dlouhý běh
first2 <- long_segs[transition_phase == "before_first"
][order(seg_id_phase)
][, .SD[1:2], by = .(fishid, navigation_event_id)]

first_segment_points <- merge(
  last_seg_sub[transition_phase == "before_first"],
  first2[, .(fishid, navigation_event_id, seg_id_phase)],
  by = c("fishid", "navigation_event_id", "seg_id_phase")
)

# RETURN (after_last): vezmi poslední 2 dlouhé běhy
last2 <- long_segs[transition_phase == "after_last"
][order(seg_id_phase)
][, .SD[(.N-1L):.N], by = .(fishid, navigation_event_id)]

last_segment_points <- merge(
  last_seg_sub[transition_phase == "after_last"],
  last2[, .(fishid, navigation_event_id, seg_id_phase)],
  by = c("fishid", "navigation_event_id", "seg_id_phase")
)

test <- first_segment_points[navigation_event_id == 64]

unique(first_segment_points$navigation_event_id)
ggplot()+
geom_sf(data = shape.rimov, fill = "grey90", color = "grey60") +
  coord_sf(xlim = c(461900, 462800), ylim = c(5410550, 5411100), expand = FALSE) +
  geom_point(data =test, aes(x = easting, y = northing, col = shore_50))
  

# 5. Spojit obě části dohromady
selected_segments <- rbind(first_segment_points, last_segment_points)

# Předpoklad: selected_segments obsahuje first/last 2 segmenty každé události
movement_summary <- selected_segments[
  order(timestamp),
  .(shore_seq = {
    vals <- na.omit(unique(shore_50_merged))
    if (length(vals) == 0) NA_character_ else paste0(vals, collapse = "-")
  }),
  by = .(fishid, navigation_event_id, transition_phase, cluster_id)
][
  , move_type := fcase(
    transition_phase == "before_first" & shore_seq == "1-0", "shore to offshore",
    transition_phase == "before_first" & shore_seq == "0-1", "offshore to shore",
    transition_phase == "before_first" & shore_seq == "1-1", "only shore",
    transition_phase == "before_first" & shore_seq == "1",   "only shore",
    transition_phase == "before_first" & shore_seq == "0",   "only offshore",
    
    transition_phase == "after_last" & shore_seq == "0-1", "offshore to shore",
    transition_phase == "after_last" & shore_seq == "1-0", "shore to offshore",
    transition_phase == "after_last" & shore_seq == "1-1", "only shore",
    transition_phase == "after_last" & shore_seq == "1",   "only shore",
    transition_phase == "after_last" & shore_seq == "0",   "only offshore",
    
    default = "unclear"
  )
]
movement_summary[fishid == "T449313_1" , fishid_new := "#1"]
movement_summary[fishid == "T449215_1" , fishid_new := "#2"]
movement_summary[, fishid_new := factor(fishid_new, levels =c("#1", "#2"))]
movement_summary[ transition_phase == "after_last" , tra_phase_new := "Return"]
movement_summary[ transition_phase == "before_first" , tra_phase_new := "Outbound"]
movement_summary[, tra_phase_new := factor(tra_phase_new, levels =c("Outbound", "Return"))]
# Přehledná tabulka
mov_table <- movement_summary[cluster_id == 1, .N, by = .(fishid,transition_phase, move_type) ]
setkey(mov_table,fishid,transition_phase,move_type)
mov_table[,.(sumN = sum(N)), by = transition_phase]

# set2_colors <- c("shore_to_offshore" = "#66c2a5",
#                  "only_shore"        = "#fc8d62",
#                  "only_offshore"     = "#8da0cb",
#                  "offshore_to_shore" = "#e78ac3")

# colour-blind
set2_colors <- c(
  "shore to offshore" = "#009E73", # green
  "only shore"        = "#E69F00", # orange
  "only offshore"     = "#0072B2", # blue
  "offshore to shore" = "#CC79A7"  # magenta
)
movement_summary

types_gg <- ggplot(movement_summary, aes(x = fishid_new, fill = move_type)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  facet_wrap(~tra_phase_new)+
  scale_fill_manual(values = set2_colors) +
  labs( y = "Percentage", fill = "Move type") +
  theme_minimal()+
  theme(
    text = element_text(size = 14),
    axis.title.x = element_blank()
  )
types_gg
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/types_gg.png",
  plot = types_gg ,
  width = 8,      # šířka v palcích
  height = 5,     # výška v palcích
  dpi = 300
)

 # vypocitani miss point u after last ####
# 1. Vezmi movement_summary pro after_last návraty
after_movement <- movement_summary[transition_phase == "after_last"]

# 2. Základní tabulka všech návratů
return_accuracy_complete <- after_movement[, .(fishid, navigation_event_id, move_type)]
return_accuracy_complete[, dist_to_rp := NA_real_]

# 3. direct_home + only_offshore → 0
return_accuracy_complete[move_type %in% c("direct home", "only offshore"),
                         dist_to_rp := 0]

only_shore_points_coords <- selected_segments[
  transition_phase == "after_last" & shore_50_merged %in% c(0,1)
][
  return_accuracy_complete[move_type == "only shore"],
  on = .(fishid, navigation_event_id),
  nomatch = 0
][
  order(timestamp),
  .SD[1],  # první bod po návratu
  by = .(fishid, navigation_event_id)
][
  , .(fishid, navigation_event_id,
      dist_to_rp = sqrt((easting - x_sleep)^2 + (northing - y_sleep)^2),
      easting, northing)
]

# 5. offshore_to_shore → poslední offshore bod (shore_50_fix == 0)
offshore_points_coords <- selected_segments[
  transition_phase == "after_last" & shore_50_merged == 0
][
  return_accuracy_complete[move_type == "offshore to shore"],
  on = .(fishid, navigation_event_id),
  nomatch = 0
][
  order(timestamp),
  .SD[.N],  # poslední offshore bod
  by = .(fishid, navigation_event_id)
][
  , .(fishid, navigation_event_id, 
      dist_to_rp = sqrt((easting - x_sleep)^2 + (northing - y_sleep)^2),
      easting, northing)
]

# 6. Slouč obě nové výpočty zpět do hlavní tabulky
return_accuracy_complete <- merge(
  return_accuracy_complete,
  rbind(only_shore_points_coords, offshore_points_coords),
  by = c("fishid", "navigation_event_id"),
  all.x = TRUE,
  suffixes = c("", "_new")
)

# 7. Vyplň chybějící hodnoty (pokud dist_to_rp původně nebylo spočteno)
return_accuracy_complete[
  is.na(dist_to_rp) & !is.na(dist_to_rp_new),
  dist_to_rp := dist_to_rp_new
][
  , dist_to_rp_new := NULL  # uklidit
]

# 8. Ověření
return_accuracy_complete[, .N]  # mělo by být 128
summary(return_accuracy_complete$dist_to_rp)

return_accuracy_complete <- merge(return_accuracy_complete, movement_summary[transition_phase == "after_last",.(cluster_id, navigation_event_id)], by = c("navigation_event_id"))
return_accuracy_cl1 <- return_accuracy_complete[cluster_id == 1]



library(ggplot2)

ggplot(return_accuracy_cl1, aes(x = move_type, y = dist_to_rp)) +
  geom_violin(fill = "grey80", color = "black", scale = "width", trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  labs(
    x = "Return move type",
    y = "Distance from RP at end of return (m)",
    title = "Accuracy of return navigation by move type"
  ) +
  facet_wrap(~fishid)+
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Spojit oba typy návratových bodů
nav_return_points <- rbind(only_shore_points_coords, offshore_points_coords)
nav_return_points <- merge(nav_return_points, movement_summary[transition_phase == "after_last",.(cluster_id, navigation_event_id)], by = c("navigation_event_id"))

# Převést na sf
nav_return_sf <- st_as_sf(nav_return_points, 
                          coords = c("easting", "northing"),
                          crs = st_crs(shape.rimov))
nav_return_sf_cl1 <- nav_return_sf %>% filter(cluster_id == 1) 
# Zobrazit na mapě
ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  geom_sf(data = nav_return_sf_cl1, aes(color = dist_to_rp), size = 2) +
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y),
             shape = 21, fill = "white", color = "black", size = 3) +
  scale_color_viridis_c(option = "plasma", name = "Dist to RP (m)") +
  labs(title = "Points from which fish navigated around srore to RP") +
  facet_wrap(~fishid)+
  coord_sf(xlim = c(461200, 462800), ylim = c(5410200, 5411100), expand = FALSE) +
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 10),
    legend.position = "right"
  )


# delka poslednich inshore a offhore segmentu ####

selected_segments_sf <- st_as_sf(
  selected_segments,
  coords = c("easting", "northing"),
  crs = 32633
)

# SI index pro offhore
# Sloučí body do trajektorie podle segmentu
segment_lengths <- selected_segments_sf |>
  group_by(fishid, navigation_event_id, transition_phase, shore_group_final, shore_50_merged) |>
  summarise(
    geometry = st_cast(st_combine(geometry), "LINESTRING"),
    .groups = "drop"
  ) |>
  mutate(
    segment_length_m = as.numeric(st_length(geometry))
  ) |>
  st_drop_geometry()

segment_lengths <- as.data.table(segment_lengths)

ggplot(segment_lengths[!is.na(shore_50_merged),], aes(x = transition_phase, y = segment_length_m, col =  as.factor(shore_50_merged))) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~fishid+transition_phase, scales = "free_x") +
  ylim(0,2000)+
  labs(
    title = "Segment Lengths of Inshore and Offshore Trajectories",
    x = "Transition Phase",
    y = "Segment Length (m)"
  ) +
  theme_minimal(base_size = 14) 

# 1. Filtrovat offshore segmenty (shore_50_fix == 0)
# 1. Připrav offshore segmenty znovu jako data.table
offshore_segments <- selected_segments[shore_50_merged == 0]

# 2. Výpočet SI pro každý offshore segment
offshore_si <- offshore_segments[
  order(timestamp),  # zajistit správné pořadí bodů
  {
    # přímá vzdálenost z prvního bodu k domovu
    straight_dist <- sqrt((easting[1] - x_sleep[.N])^2 + (northing[1] - y_sleep[.N])^2)
    
    # délka celé trajektorie
    total_path <- sum(sqrt(diff(easting)^2 + diff(northing)^2), na.rm = TRUE)
    
    si <- if (total_path > 0) straight_dist / total_path else NA_real_
    
    list(
      straight_distance_m = straight_dist,
      total_path_length_m = total_path,
      straightness_index = si
    )
  },
  by = .(fishid, navigation_event_id, transition_phase, shore_group_final)
]

ggplot(offshore_si, aes(x = transition_phase, y = straightness_index, fill = transition_phase)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~fishid) +
  labs(
    title = "Straightness Index of Offshore Segments by Fish",
    x = "Transition Phase",
    y = "Straightness Index (SI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



# donavigace ####
library(data.table)
library(sf)
library(geosphere)
library(circular)
library(broom)

library(data.table)
library(sf)
library(geosphere)
library(circular)

 # last_offshore <- offshore_segments[!is.na(shore_50_fix) & shore_50_merged== 0 & transition_phase == "after_last"
 #   , .SD[shore_merged_group == max(shore_merged_group, na.rm = TRUE)], 
 #   by = .(fishid, navigation_event_id)
 # ]
# last_seg_sub[, change := c(TRUE, shore_50_fix[-1] != shore_50_fix[-.N]), 
#              by = .(fishid, navigation_event_id)]
# 
# last_offshore <- offshore_segments[shore_50_merged== 0 & transition_phase == "after_last"
#    , .SD[shore01_clean == max(shore01_clean, na.rm = TRUE)], 
#    by = .(fishid, navigation_event_id)
#  ]

last_offshore <- offshore_segments[
  shore_50_merged == 0 & transition_phase == "after_last",
  .SD[shore_group_final == max(shore_group_final, na.rm = TRUE)],
  by = .(fishid, navigation_event_id)
]
 ggplot(last_offshore[cluster_id ==1 & fishid == "T449313_1"], aes(x = easting, y = northing, group = navigation_event_id))+
   geom_path(size = 0.1)
 
# 1. Spočítej 3 posuny
last_offshore[, `:=`(
  lead5_easting = shift(easting, type = "lead", n = 5),
  lead5_northing = shift(northing, type = "lead", n = 5)
), by = .(fishid, navigation_event_id)]

# Bearing z bodu i → i+4
valid_rows <- !is.na(last_offshore$lead5_easting) & !is.na(last_offshore$lead5_northing)

# převod na sf a WGS
sf_from <- st_as_sf(last_offshore[valid_rows], coords = c("easting", "northing"), crs = 32633)
sf_to   <- st_as_sf(last_offshore[valid_rows], coords = c("lead5_easting", "lead5_northing"), crs = 32633)
sf_from <- st_transform(sf_from, 4326)
sf_to   <- st_transform(sf_to, 4326)

# výpočet bearingu
bearings <- geosphere::bearing(
  st_coordinates(sf_from),
  st_coordinates(sf_to)
)

# přiřazení zpět
last_offshore[valid_rows, bearing_movement := bearings]

# Funkce pro výpočet bearingu (ve stupních) mezi dvěma body
bearing_to_home <- function(x1, y1, x2, y2) {
  atan2(x2 - x1, y2 - y1) * 180 / pi
}

# Přidej směr ke "home_center"
last_offshore[, bearing_to_home := bearing_to_home(easting, northing, x_sleep, y_sleep)]


# 7. Odchylka od ideálního směru
last_offshore[, bearing_diff := (bearing_to_home - bearing_movement + 180) %% 360 - 180]

# 8. Čas od začátku návratu
last_offshore[, time_from_start := as.numeric(difftime(timestamp, min(timestamp), units = "secs")),
              by = .(fishid, navigation_event_id)]

# 9. Model: |bearing_diff| ~ čas
bearing_lm <- last_offshore[
  !is.na(bearing_diff),
  {
    model <- lm(abs(bearing_diff) ~ time_from_start)
    slope <- coef(model)[["time_from_start"]]
    list(slope = slope)
  },
  by = .(fishid, navigation_event_id)
]

# 10. Délka trajektorie
segment_lengths <- last_offshore[
  !is.na(bearing_diff),
  .(segment_length_m = sum(sqrt(diff(easting)^2 + diff(northing)^2), na.rm = TRUE)),
  by = .(fishid, navigation_event_id, cluster_id)
]

# 11. Spojení
bearing_lm <- merge(bearing_lm, segment_lengths, by = c("fishid", "navigation_event_id"))


ggplot(bearing_lm[cluster_id == 1] , aes(x = segment_length_m, y = slope)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~fishid)+
  labs(
    title = "Slope of |bearing_diff| vs. segment length",
    x = "Offshore segment length (m)",
    y = "Slope of |Δbearing| over time"
  ) +
  theme_minimal()
# 7. Vybrat slope
#bearing_slope <- bearing_lm[term == "time_rel_sec"]
#nav_return_points <- merge(nav_return_points, movement_summary[transition_phase == "after_last",.(cluster_id, navigation_event_id)], by = c("navigation_event_id"))

bearing_lm[cluster_id == 1 & fishid == "T449215_1"]

ggplot(bearing_lm[cluster_id == 1] , aes(x = fishid, y = slope, fill = fishid)) +
  geom_boxplot() +
  labs(
    title = "Slope of |Bearing Difference| over Time in Offshore Segments",
    x = "Fish ID",
    y = "Slope (change in |ΔBearing| per second)"
  ) +
  theme_minimal(base_size = 14)

ggplot(last_offshore[cluster_id == 1 & fishid == "T449215_1"], aes(x = time_from_start, y = abs(bearing_diff)))+
         geom_point()+
         geom_smooth(method= "lm")+
  facet_wrap(~navigation_event_id, scales = "free")
       
ggplot(last_offshore[navigation_event_id == 138], 
       aes(x = time_from_start, y = abs(bearing_diff))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, col = "blue") +
  labs(title = "Bearing difference over time (Event 138)", x = "Čas do RP (s)", y = "|bearing_diff|")
# vytvoř sf body pro pohyb
traj_sf <- st_as_sf(last_offshore[navigation_event_id == 72], coords = c("easting", "northing"), crs = 32633)

# spánkové místo
home_pt <- last_offshore[navigation_event_id == 72, .SD[1], by = fishid]
home_sf <- st_as_sf(home_pt, coords = c("x_sleep", "y_sleep"), crs = 32633)

# plot
ggplot() +
  geom_sf(data = traj_sf, aes(color = time_rel_sec)) +
  geom_sf(data = home_sf, shape = 21, fill = "black", size = 3) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  labs(title = "Trajektorie návratu a vzdálenost od cíle", color = "Čas do RP") +
  coord_sf()


library(data.table)
library(zoo)
library(pracma)

# ---------- Pomocné funkce ----------
ang_diff_signed <- function(a, b) {
  # podepsaný rozdíl úhlů (stupně) v intervalu (-180, 180]
  d <- a - b
  ((((d + 180) %% 360) - 180))
}

prep_bearing_error <- function(dt, window_size = 7) {
  dt <- copy(dt)
  
  # čas od začátku
  if (!"time_from_start" %in% names(dt)) {
    dt[, time_from_start := as.numeric(difftime(timestamp, min(timestamp), units = "secs")),
       by = .(fishid, navigation_event_id)]
  }
  
  # podepsaná chyba směru a její vyhlazení
  dt[, bearing_err_signed := ang_diff_signed(bearing_movement, bearing_to_home)]
  dt[, bearing_err_smooth := zoo::rollmean(bearing_err_signed, k = window_size,
                                           fill = NA, align = "center")]
  dt[, bearing_err_abs_s  := abs(bearing_err_smooth)]
  
  # vzdálenost mezi po sobě jdoucími body (pokud chybí)
  if (!"distance_from_previous" %in% names(dt)) {
    dt[, distance_from_previous :=
         sqrt((easting - shift(easting))^2 + (northing - shift(northing))^2),
       by = .(fishid, navigation_event_id)]
  }
  
  dt
}

detect_donavigation <- function(dt,
                                window_size = 10,
                                min_delta = 5,
                                min_duration = 30,
                                min_points = 3,
                                max_correction_duration = 600,
                                debug = TRUE) {
  dt <- data.table::copy(dt)
  
  # 0) vstupy – očekáváme distance_to_home a distance_from_previous
  if (!all(c("distance_to_home","distance_from_previous") %in% names(dt))) {
    stop("Chybí sloupce distance_to_home a/nebo distance_from_previous.")
  }
  
  # 1) bearing_diff_abs + smoothing
  if (!"bearing_diff_abs" %in% names(dt)) {
    dt[, bearing_diff_abs := abs(bearing_diff)]
  }
  dt[, bearing_diff_smooth := zoo::rollmean(bearing_diff_abs, k = window_size,
                                            fill = NA, align = "center")]
  
  # 2) čas od startu
  if (!"time_from_start" %in% names(dt)) {
    dt[, time_from_start := as.numeric(
      difftime(timestamp, min(timestamp), units = "secs")
    ), by = .(fishid, navigation_event_id)]
  }
  
  y <- dt$bearing_diff_smooth
  t <- dt$time_from_start
  
  # 3) extrémy
  peaks   <- pracma::findpeaks(y, minpeakdistance = 2)
  valleys <- pracma::findpeaks(-y, minpeakdistance = 2)
  if (is.null(peaks) | is.null(valleys)) {
    if (debug) message("Žádné peaky/valleys.")
    return(NULL)
  }
  peaks_idx   <- as.integer(peaks[,2])
  valleys_idx <- as.integer(valleys[,2])
  all_extrema <- sort(unique(c(peaks_idx, valleys_idx)))
  
  # 4) sběr úseků peak -> následující valley
  corrections <- list()
  phase_id <- 1
  i <- 1
  while (i < length(all_extrema)) {
    start_idx <- all_extrema[i]
    if (!(start_idx %in% peaks_idx)) { i <- i + 1; next }
    
    # kandidátní minima v časovém limitu
    possible_ends <- all_extrema[(i+1):length(all_extrema)]
    possible_ends <- possible_ends[possible_ends %in% valleys_idx]
    possible_ends <- possible_ends[ t[possible_ends] - t[start_idx] <= max_correction_duration ]
    
    if (!length(possible_ends)) { i <- i + 1; next }
    end_idx <- possible_ends[1]
    
    # hrubá filtrace
    delta <- y[start_idx] - y[end_idx]
    dur   <- t[end_idx] - t[start_idx]
    npt   <- end_idx - start_idx + 1
    if (is.na(delta) || is.na(dur)) { i <- which(all_extrema == end_idx) + 1; next }
    if (delta < min_delta || dur < min_duration || npt < min_points) {
      i <- which(all_extrema == end_idx) + 1
      next
    }
    
    # 5) METRIKY ÚČINNOSTI
    # vzdálenosti k RP
    start_dist <- dt$distance_to_home[start_idx]
    end_dist   <- dt$distance_to_home[end_idx]
    # celková uražená dráha v úseku
    total_dist <- sum(dt$distance_from_previous[start_idx:end_idx], na.rm = TRUE)
    # čistý zisk směrem k RP (m)
    home_gain_m <- start_dist - end_dist
    
    # podíl kroků, které zmenšují distance_to_home
    steps <- dt[(start_idx+1):end_idx,
                .(d = distance_to_home - shift(distance_to_home)), # změna vůči předchozímu
                by = .(idx = .I + start_idx)]
    steps <- steps[!is.na(d)]
    progress_ratio <- if (nrow(steps)) mean(steps$d < 0) else NA_real_
    
    # normalizace zisku (0–1+, záporné pokud se vzdálil)
    gain_ratio <- if (is.finite(start_dist) && start_dist > 0)
      home_gain_m / start_dist else NA_real_
    
    # efektivita trajektorie (kolik čistého přiblížení na 1 m dráhy)
    path_efficiency <- if (is.finite(total_dist) && total_dist > 0)
      home_gain_m / total_dist else NA_real_
    
    # jednočíselné skóre (váhy můžeš kdykoli upravit)
    score_efficiency <- 0.5 * gain_ratio +
      0.3 * progress_ratio +
      0.2 * path_efficiency
    
    # metriky korekce
    dist_m <- total_dist
    rate_deg_per_s <- delta / dur
    
    corrections <- append(corrections, list(data.table(
      fishid              = unique(dt$fishid),
      navigation_event_id = unique(dt$navigation_event_id),
      correction_phase_id = phase_id,
      start_idx, end_idx,
      start_time = t[start_idx],
      end_time   = t[end_idx],
      delta_bearing = delta,
      duration   = dur,
      distance_m = dist_m,
      rate_deg_per_s = rate_deg_per_s,
      # účinnost donavigace:
      start_dist = start_dist,
      end_dist   = end_dist,
      home_gain_m = home_gain_m,
      gain_ratio = gain_ratio,
      progress_ratio = progress_ratio,
      path_efficiency = path_efficiency,
      score_efficiency = score_efficiency
    )))
    
    phase_id <- phase_id + 1
    i <- which(all_extrema == end_idx) + 1
  }
  
  if (!length(corrections)) {
    if (debug) message("Nebyl detekován žádný korekční úsek.")
    return(NULL)
  }
  
  corrections <- data.table::rbindlist(corrections)
  
  # 6) flag do dat
  dt[, correction_phase_id := NA_integer_]
  for (j in seq_len(nrow(corrections))) {
    dt[corrections$start_idx[j]:corrections$end_idx[j],
       correction_phase_id := corrections$correction_phase_id[j]]
  }
  
  if (debug) {
    message("Detekováno ", nrow(corrections), " korekčních úseků.")
  }
  
  list(
    corrections = corrections[],
    dt = dt[]
  )
}

score_trajectory_global <- function(dt, res = NULL, hit_radius = 20, bearing_window = 7) {
  stopifnot(nrow(dt) > 0)
  dt_corr <- if (!is.null(res$dt)) data.table::copy(res$dt) else data.table::copy(dt)
  
  # pořadí a čas
  if (!"time_from_start" %in% names(dt_corr)) {
    dt_corr[, time_from_start := as.numeric(difftime(timestamp, min(timestamp), units = "secs")),
            by = .(fishid, navigation_event_id)]
  }
  data.table::setorder(dt_corr, time_from_start)
  
  # segmentové vzdálenosti
  if (!"distance_from_previous" %in% names(dt_corr)) {
    dt_corr[, distance_from_previous :=
              sqrt((easting - data.table::shift(easting))^2 +
                     (northing - data.table::shift(northing))^2)]
  }
  
  # bearing smooth (pro bearing_work)
  if (!"bearing_diff_abs" %in% names(dt_corr))
    dt_corr[, bearing_diff_abs := abs(bearing_diff)]
  dt_corr[, bearing_diff_smooth :=
            zoo::rollmean(bearing_diff_abs, k = bearing_window, fill = NA, align = "center")]
  
  # globální metriky
  start_dist <- dt_corr$distance_to_home[which.min(dt_corr$time_from_start)]
  end_dist   <- dt_corr$distance_to_home[which.max(dt_corr$time_from_start)]
  net_home_gain_m <- start_dist - end_dist
  
  total_path_m <- sum(dt_corr$distance_from_previous, na.rm = TRUE)
  path_efficiency_global <- if (is.finite(total_path_m) && total_path_m > 0)
    net_home_gain_m / total_path_m else NA_real_
  
  min_distance_to_home <- suppressWarnings(min(dt_corr$distance_to_home, na.rm = TRUE))
  hit_any <- any(dt_corr$distance_to_home <= hit_radius, na.rm = TRUE)
  t_first_hit <- if (hit_any)
    suppressWarnings(min(dt_corr$time_from_start[dt_corr$distance_to_home <= hit_radius], na.rm = TRUE))
  else NA_real_
  
  bearing_work <- sum(abs(diff(dt_corr$bearing_diff_smooth)), na.rm = TRUE)
  
  # korekční část (pokud je)
  corr_path_len_m <- NA_real_
  corr_home_gain_m <- NA_real_
  corr_contrib_gain <- NA_real_
  time_in_correction_ratio <- NA_real_
  
  if (!is.null(res) && !is.null(res$corrections) && nrow(res$corrections) > 0) {
    idx_mask <- rep(FALSE, nrow(dt_corr))
    for (j in seq_len(nrow(res$corrections))) {
      idx_mask[res$corrections$start_idx[j]:res$corrections$end_idx[j]] <- TRUE
    }
    dt_corr[, in_correction := idx_mask]
    
    corr_path_len_m <- sum(dt_corr[in_correction == TRUE]$distance_from_previous, na.rm = TRUE)
    
    # kladný příspěvek jen tehdy, když se vzdálenost zmenšuje
    if ("distance_to_home" %in% names(dt_corr)) {
      dth <- dt_corr$distance_to_home
      dth_corr <- dth[dt_corr$in_correction == TRUE]
      corr_home_gain_m <- sum(pmax(0, -diff(dth_corr)), na.rm = TRUE)
      if (is.finite(net_home_gain_m) && net_home_gain_m > 0)
        corr_contrib_gain <- corr_home_gain_m / net_home_gain_m
    }
    
    # podíl času v korekcích – bezpečně
    dt_corr[, time_step := c(NA_real_, diff(time_from_start))]
    time_in_correction <- sum(dt_corr[in_correction == TRUE]$time_step, na.rm = TRUE)
    total_time <- diff(range(dt_corr$time_from_start, na.rm = TRUE))
    time_in_correction_ratio <-
      if (is.finite(total_time) && total_time > 0) time_in_correction / total_time else NA_real_
  }
  
  data.table::data.table(
    fishid                  = unique(dt_corr$fishid),
    navigation_event_id     = unique(dt_corr$navigation_event_id),
    start_dist              = start_dist,
    end_dist                = end_dist,
    net_home_gain_m         = net_home_gain_m,
    total_path_m            = total_path_m,
    path_efficiency_global  = path_efficiency_global,
    progress_ratio_global   = NA_real_,  # nechávám volitelné, pokud ho počítáš jinde
    min_distance_to_home    = min_distance_to_home,
    t_first_hit             = t_first_hit,
    bearing_work            = bearing_work,
    corr_path_len_m         = corr_path_len_m,
    corr_home_gain_m        = corr_home_gain_m,
    corr_contrib_gain       = corr_contrib_gain,
    time_in_correction_ratio= time_in_correction_ratio
  )
}

dt <- last_offshore[navigation_event_id == 23]

res <- detect_donavigation_signed(
  dt,
  window_size = 2,
  min_delta = 5,
  min_duration = 10,
  min_points = 3,
  max_correction_duration = 600,
  min_home_gain = 5,         # m
  min_progress_ratio = 0.6   # ≥60 % kroků blíž k RP
)


res$corrections  # metriky donavigačních úseků
# res$dt         # původní data + bearing_err_* + correction_flag + correction_phase_id
score <- score_trajectory_global(dt, res = res, hit_radius = 20, bearing_window = 7)
score


stopifnot(!is.null(res))

dt_corr <- copy(res$dt)
corrs   <- res$corrections

# 1) Zapiš correction_phase_id do dt_corr podle indexů
dt_corr[, correction_phase_id := NA_integer_]
if (nrow(corrs)) {
  for (j in seq_len(nrow(corrs))) {
    dt_corr[corrs$start_idx[j]:corrs$end_idx[j],
            correction_phase_id := corrs$correction_phase_id[j]]
  }
}

# 2) Připrav start/end body
starts <- dt_corr[corrs$start_idx,
                  .(fishid, navigation_event_id, easting, northing,
                    time_from_start,
                    correction_phase_id = corrs$correction_phase_id)]
ends   <- dt_corr[corrs$end_idx,
                  .(fishid, navigation_event_id, easting, northing,
                    time_from_start,
                    correction_phase_id = corrs$correction_phase_id)]
library(ggplot2)
library(ggrepel)

ggplot() +
  # celá trajektorie šedě
  geom_path(data = dt_corr[order(time_from_start)],
            aes(easting, northing, group = 1), color = "grey70", linewidth = 0.6) +
  geom_point(data = dt_corr[order(time_from_start)],
             aes(easting, northing, group = 1), color = "grey70", size = 0.5) +
  geom_point(data = dt_corr[order(time_from_start)],
             aes(mean(x_sleep), mean(y_sleep), group = 1), color = "black", size = 2, shape=4) +
  # donavigační úseky barevně (podle ID korekce)
  geom_path(data = dt_corr[!is.na(correction_phase_id)][order(time_from_start)],
            aes(easting, northing, group = correction_phase_id, color = factor(correction_phase_id)),
            linewidth = 1.2) +
  
  # začátky a konce
  geom_point(data = starts, aes(easting, northing, color = factor(correction_phase_id)),
             size = 3, shape = 21, fill = "white", stroke = 1.2) +
  geom_point(data = ends,   aes(easting, northing, color = factor(correction_phase_id)),
             size = 3, shape = 16) +
  
  # popisky start/end (volitelné)
  geom_label_repel(data = starts,
                   aes(easting, northing, label = paste0("start#", correction_phase_id),
                       color = factor(correction_phase_id)),
                   size = 3, show.legend = FALSE) +
  geom_label_repel(data = ends,
                   aes(easting, northing, label = paste0("end#", correction_phase_id),
                       color = factor(correction_phase_id)),
                   size = 3, show.legend = FALSE) +
  
  coord_equal() +
  labs(color = "Correction ID",
       title = "Donavigační úseky (peak → min)",
       x = "Easting (m, UTM 33N)", y = "Northing (m, UTM 33N)") +
  theme_minimal(base_size = 13)


# orintace pri navratu z ruznych smeru ####
last_offshore

ggplot(last_offshore[cluster_id == 1], aes(x = easting, y = northing, group = navigation_event_id))+
  geom_path()+
  
  facet_wrap(~fishid)
  

last_offshore[, distance_from_previous := sqrt((easting - data.table::shift(easting))^2 + 
                                                (northing - data.table::shift(northing))^2), 
             by = .(fishid, navigation_event_id)]
last_offshore[!is.na(distance_from_previous) , cum_dist := cumsum(distance_from_previous), by = .(fishid, navigation_event_id)]
# 1. Přepočítej bearing_diff do rozsahu 0–360
last_offshore[, bearing_diff_360 := (bearing_movement - bearing_to_home) %% 360]

bearing_diff_summary <- last_offshore[
  cum_dist <= 100 & !is.na(bearing_diff_360),
  .(mean_bearing_diff = as.numeric(
    circular::mean.circular(
      circular::circular(bearing_diff_360, units = "degrees", template = "geographics")
    )
  )),
  by = .(fishid, navigation_event_id, cluster_id)
]

ggplot(bearing_diff_summary, aes(x = mean_bearing_diff))+geom_histogram()+facet_wrap(~fishid)


# ---- 1. Příprava dat (používáme unikátní název dt_8b) ----
dt_8b <- as.data.table(bearing_diff_summary)

# Odfiltrování případných NA
dt_8b <- dt_8b[!is.na(mean_bearing_diff) & cluster_id  == 1]

# Mapování ID ryb na labely použité v článku
dt_8b[, fishid_new := fifelse(fishid == "T449313_1", "#1",
                              fifelse(fishid == "T449215_1", "#2", fishid))]
dt_8b[, fishid_new := factor(fishid_new, levels = c("#1", "#2"))]

# ---- n do facetů ----
counts_8b <- dt_8b[, .N, by = fishid_new]
counts_8b[, lab := paste0(as.character(fishid_new), " (n = ", N, ")")]
counts_8b[, fishid_new_n :=
            factor(lab,
                   levels = lab[match(c("#1", "#2"), fishid_new)])]

# připojit label s n do hlavních dat
dt_8b <- merge(dt_8b, counts_8b[, .(fishid_new, fishid_new_n)], by = "fishid_new", all.x = TRUE)

# Úprava úhlu: z rozsahu [0, 360] na [-180, 180]
# 0 = direct home, + vpravo, - vlevo
dt_8b[, ang := mean_bearing_diff]
dt_8b[, ang := ((ang + 180) %% 360) - 180]

# Parametry pro graf
bw <- 15 # šířka binu

# ---- 2. Histogram / Růžicový podklad ----
hist_8b <- dt_8b[, .N, by = .(fishid_new, fishid_new_n, bin = floor(ang / bw) * bw)]
hist_8b[, prop := N / sum(N), by = fishid_new]

# ---- 3. Cirkulární statistika (Průměr a délka vektoru R) ----
means_8b <- dt_8b[, .(
  mean_rad = atan2(sum(sin(ang * pi/180)), sum(cos(ang * pi/180))),
  R = sqrt((sum(cos(ang * pi/180)))^2 + (sum(sin(ang * pi/180)))^2) / .N
), by = .(fishid_new, fishid_new_n)]

means_8b[, mean_deg := mean_rad * 180/pi]
means_8b[, mean_deg := ((mean_deg + 180) %% 360) - 180]

# ---- 4. Bootstrap 95% CI (používá tvou definovanou funkci boot_ci_one) ----
set.seed(1)
ci_8b <- dt_8b[, boot_ci_one(ang, B = 2000), by = .(fishid_new, fishid_new_n)]

# Spojení statistik do jedné tabulky
stats_8b <- merge(means_8b, ci_8b, by = c("fishid_new", "fishid_new_n"), all.x = TRUE)

# ---- 5. Generování oblouku CI pro graf (používá tvou funkci make_arc) ----
arc_8b <- rbindlist(lapply(1:nrow(stats_8b), function(i){
  f  <- stats_8b$fishid_new[i]
  fn <- stats_8b$fishid_new_n[i]
  make_arc(stats_8b$ci_low[i], stats_8b$ci_high[i], y = 0.52, n = 120)[, `:=`(fishid_new = f, fishid_new_n = fn)]
}))

# ---- 6. Vykreslení grafu 8b ----
p_8b <- ggplot() +
  # Růžicový histogram v pozadí
  geom_col(
    data = hist_8b,
    aes(x = bin, y = prop),
    fill = "blue",
    color = "grey70",
    width = bw,
    linewidth = 0.25, alpha = 0.2
  ) +
  # Surová data (body na vnějším okraji)
  geom_point(
    data = dt_8b,
    aes(x = ang, y = 0.50, color = fishid_new),
    alpha = 0.55,
    size = 2.2,
    position = position_jitter(width = 0, height = 0.015)
  ) +
  # Průměrný vektor (délka škálovaná podle R)
  geom_segment(
    data = stats_8b,
    aes(x = mean_deg, xend = mean_deg, y = 0, yend = 0.45 * R),
    linewidth = 1.2,
    color = "black",
    lineend = "round",
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  # Modrý oblouk konfidenčního intervalu
  geom_path(
    data = arc_8b,
    aes(x = x, y = y),
    linewidth = 2,
    color = "blue",
    lineend = "round"
  ) +
  # Referenční čára pro přímý směr domů
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  # Nastavení polárních koordinát
  coord_polar(theta = "x", start = pi, direction = 1) +
  scale_x_continuous(
    breaks = c(-90, 0, 90, 180),
    labels = c("left", "direct", "right", "opposite"),
    limits = c(-180, 180),
    expand = c(0, 0)
  ) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0, 0.55), expand = c(0, 0)) +
  facet_wrap(~fishid_new_n, nrow = 1) +
  guides(color = "none") +
  theme_minimal(base_size = 15) +
  theme(
    text = element_text(face = "bold", size = 15),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    clip = "off"
  )

print(p_8b)


# new version 
# ------------------------------------------------------------
# Combined orientation plot (rose + raw points + mean + 95% CI)
# Faceted for both focal individuals (#1 and #2)
# Uses the SAME angular convention as your old rose plots:
#   0 = direct (top), + = right, - = left, range [-180, 180]
# ------------------------------------------------------------

library(data.table)
library(ggplot2)
results_detour <- results_detour %>%
  mutate(
    empirical_bearing_360 = (empirical_bearing + 360) %% 360,
    lcp_bearing_360 = (lcp_bearing + 360) %% 360,
    bearing_offset = (empirical_bearing_360 - lcp_bearing_360 + 360) %% 360
  )


# ---- prepare data ----
dt <- as.data.table(results_detour)

# keep only the phase you want
dt <- dt[transition_phase == "after_last" & !is.na(bearing_diff_adj) &cluster_id == 1]

# map fish IDs -> labels used in the paper
dt[, fishid_new := fifelse(fishid == "T449313_1", "#1",
                           fifelse(fishid == "T449215_1", "#2", fishid))]
dt[, fishid_new := factor(fishid_new, levels = c("#1", "#2"))]

dt <- dt[!duplicated(dt, by = c("fishid","navigation_event_id"))]

# ---- add n for facet labels ----
counts_dt <- dt[, .N, by = fishid_new]
counts_dt[, lab := paste0(as.character(fishid_new), " (n = ", N, ")")]

# enforce facet order (#1 then #2)
counts_dt <- counts_dt[match(c("#1", "#2"), fishid_new), ]
counts_dt[, fishid_new_n := factor(lab, levels = lab)]

# join facet label back
dt <- merge(dt, counts_dt[, .(fishid_new, fishid_new_n)], by = "fishid_new", all.x = TRUE)

# signed angle in [-180,180] with your definition:
# 0 = direct; + right; - left
dt[, ang := bearing_offset]
dt[, ang := ((ang + 180) %% 360) - 180]   # wrap, kdyby náhodou byly mimo

bw <- 15

# histogram (rose background) – BINUJ ang!
hist_dt <- dt[, .N, by=.(fishid_new, fishid_new_n, bin = floor(ang / bw) * bw)]
hist_dt[, prop := N / sum(N), by=fishid_new]

# ---- circular mean + vector length (R) from RAW angles (no binning) ----
means_dt <- dt[, .(
  mean_rad = atan2(sum(sin(ang * pi/180)), sum(cos(ang * pi/180))),
  R = sqrt((sum(cos(ang * pi/180)))^2 + (sum(sin(ang * pi/180)))^2) / .N
), by = .(fishid_new, fishid_new_n)]
means_dt[, mean_deg := mean_rad * 180/pi]
means_dt[, mean_deg := ((mean_deg + 180) %% 360) - 180]

# ---- bootstrap 95% CI of mean direction (per fish) ----
set.seed(1)

boot_ci_one <- function(x_deg, B = 2000){
  x_rad <- x_deg * pi/180
  boot_mean <- replicate(B, {
    xr <- sample(x_rad, replace = TRUE)
    atan2(sum(sin(xr)), sum(cos(xr)))
  })
  boot_deg <- boot_mean * 180/pi
  boot_deg <- ((boot_deg + 180) %% 360) - 180
  
  mu <- atan2(sum(sin(x_rad)), sum(cos(x_rad))) * 180/pi
  mu <- ((mu + 180) %% 360) - 180
  
  d <- boot_deg - mu
  d <- ((d + 180) %% 360) - 180
  
  q <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
  ci_low  <- mu + q[[1]]
  ci_high <- mu + q[[2]]
  ci_low  <- ((ci_low  + 180) %% 360) - 180
  ci_high <- ((ci_high + 180) %% 360) - 180
  
  data.table(ci_low = ci_low, ci_high = ci_high)
}

ci_dt <- dt[, boot_ci_one(ang, B = 2000), by = .(fishid_new, fishid_new_n)]

stats_dt <- merge(means_dt, ci_dt, by = c("fishid_new", "fishid_new_n"), all.x = TRUE)

# ---- helper: generate CI arc points robustly even if CI crosses -180/180 ----
make_arc <- function(low, high, y, n = 120){
  low0  <- (low  + 360) %% 360
  high0 <- (high + 360) %% 360
  if (high0 < low0) high0 <- high0 + 360
  xs <- seq(low0, high0, length.out = n)
  xs2 <- ((xs + 180) %% 360) - 180
  data.table(x = xs2, y = y)
}

arc_dt <- rbindlist(lapply(1:nrow(stats_dt), function(i){
  f  <- stats_dt$fishid_new[i]
  fn <- stats_dt$fishid_new_n[i]
  make_arc(stats_dt$ci_low[i], stats_dt$ci_high[i], y = 0.52, n = 120)[, `:=`(fishid_new = f, fishid_new_n = fn)]
}))

# ---- plot ----
p <- ggplot() +
  geom_col(
    data = hist_dt,
    aes(x = bin, y = prop),
    fill = "blue",
    color = "grey70",
    width = bw,
    linewidth = 0.25, alpha = 0.2
  ) +
  geom_point(
    data = dt,
    aes(x = ang, y = 0.50, color = fishid_new),
    alpha = 0.55,
    size = 2.2,
    position = position_jitter(width = 0, height = 0.015)
  ) +
  geom_segment(
    data = stats_dt,
    aes(x = mean_deg, xend = mean_deg, y = 0, yend = 0.45 * R),
    linewidth = 1.2,
    color = "black",
    lineend = "round",
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_path(
    data = arc_dt,
    aes(x = x, y = y),
    linewidth = 2,
    color = "blue",
    lineend = "round"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_polar(theta = "x", start = pi, direction = 1) +
  scale_x_continuous(
    breaks = c(-90, 0, 90, 180),
    labels = c("left", "direct", "right", "opposite"),
    limits = c(-180, 180),
    expand = c(0, 0)
  ) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0, 0.55), expand = c(0, 0)) +
  facet_wrap(~fishid_new_n, nrow = 1) +
  guides(color = "none") +
  theme_minimal(base_size = 15) +
  theme(
    text = element_text(face = "bold", size = 15),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    clip = "off"
  )

p

stat_table_begR <- dt[,.(angle = bearing_offset, fishid_new, phase = "begR")]
stat_table_begR[, rad_ang := angle  * pi / 180]

RESPONSE_begR <- cbind(cos(stat_table_begR$rad_ang), sin(stat_table_begR$rad_ang))

model_intercept <- manova(RESPONSE_begR ~ 1)
summary(model_intercept, intercept = TRUE, test = "Pillai")


model_fish <- manova(RESPONSE_begR  ~ fishid_new, data = stat_table_begR)
summary(model_fish)

stat_table_vv <- dt_8b[,.(angle = ang, fishid_new, phase = "vv")]
stat_table_vv[, rad_ang := angle * pi / 180]

RESPONSE_vv <- cbind(cos(stat_table_vv$rad_ang), sin(stat_table_vv$rad_ang))

model_intercept_vv <- manova(RESPONSE_vv ~ 1)
summary(model_intercept_vv, intercept = TRUE, test = "Pillai")

model_fish_vv <- manova(RESPONSE_vv  ~ fishid_new, data = stat_table_vv)
summary(model_fish_vv)

stat_table <- rbind(stat_table_begR, stat_table_vv)
RESPONSE <- cbind(cos(stat_table$rad_ang), sin(stat_table$rad_ang))

model_comparison <- manova(RESPONSE ~ phase + fishid_new, data = stat_table)
summary(model_comparison, test = "Pillai")

model_comparison_2 <- manova(RESPONSE ~ phase:fishid_new, data = stat_table)
summary(model_comparison_2, test = "Pillai")

anova(model_comparison, model_comparison_2)

# 1. Získej rezidua
resids <- residuals(model_comparison)

# 2. Test multivariační normality
# install.packages("MVN")
library(MVN)
mvn_result <- MVN::mvn(resids, mvn_test = "mardia")
mvn_result$multivariate_normality
test_mvn <- mvn(resids, mvn_test  = "mardia")
plot(test_mvn)
# 3. Boxův M test (homogenita rozptylu)
# install.packages("biotools")
library(biotools)
box_res <- boxM(resids, stat_table$fishid_new)
print(box_res)

library(cowplot)

pa <- ggdraw(p) +
  draw_label("a)", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

pb <- ggdraw(p_8b) +
  draw_label("b)", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

 bearinf_return_gg_new <- plot_grid(pa, pb, ncol = 1)
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/bearinf_return_gg_new.png",
  plot =  bearinf_return_gg_new ,
  width = 9,      # šířka v palcích
  height = 10,     # výška v palcích
  dpi = 300
)

# vizualizace last offshore ####

last_offshore
last_offshore[fishid == "T449313_1", fishid_new := "#1"]
last_offshore[fishid == "T449215_1", fishid_new := "#2"]
last_offshore_tab <- last_offshore[,.(dist_m = max(cum_dist, na.rm = T)), by = .(fishid_new, navigation_event_id)]
last_offshore_med <- last_offshore_tab[,.(med_m = median(dist_m ), max_m = max(dist_m )), by = fishid_new]

last_off_dist_gg <- ggplot(last_offshore_tab, aes(x = fishid_new, y = dist_m, fill = fishid_new)) +
  geom_violin(trim = FALSE, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, aes(color = fishid_new)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c( "#00BFC4","#F8766D")) +
  scale_color_manual(values = c( "#00BFC4", "#F8766D")) +
  labs(
    y = "Length of last offshore segment (m)",
    title = "a)"
  ) +
  theme_minimal(base_size = 15)+
  theme(
    text = element_text(size = 15),
    #plot.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )
last_off_dist_gg

last_off_map_gg <- ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3, alpha = 0.3) +
  geom_path(data = last_offshore, aes(x = easting, y = northing, group = navigation_event_id), size = 0.5) +
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y),
             shape = 21, fill = "white", color = "black", size = 3) +
  scale_color_viridis_c(option = "plasma", name = "Dist to RP (m)") +
  labs(title = "b)") +
  facet_wrap(~fishid_new)+
  coord_sf(xlim = c(461200, 462800), ylim = c(5410100, 5411100), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3) +
  theme_minimal()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
   # plot.title = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.title = element_text(size = 15),
    legend.position = "right"
  )
last_off_map_gg 

last_offshore_gg <- last_off_dist_gg/last_off_map_gg 
last_offshore_gg
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/last_offshore_gg.png",
  plot =  last_offshore_gg ,
  width = 9,      # šířka v palcích
  height = 8,     # výška v palcích
  dpi = 300
)

# Výběr první pozice pro každý navigation_event_id
first_points <- last_offshore[
  cum_dist <= 100 & !is.na(bearing_movement),
  .SD[1],
  by = .(fishid, navigation_event_id, cluster_id)
]

# Data z prvních 100 m
data_first100 <- last_offshore[
  cum_dist <= 100 & !is.na(easting) & !is.na(northing)
]

# Seřadit body správně v čase
setorder(data_first100, fishid, navigation_event_id, timestamp)

# Spočítat rozdíly mezi body (použijeme shift)
data_first100[, `:=`(
  dx = data.table::shift(easting, type = "lead") - easting,
  dy = data.table::shift(northing, type = "lead") - northing
), by = .(fishid, navigation_event_id)]

# Vypočet absolutního směru pohybu v radiánech → pak převod na stupně
data_first100[, bearing_absolute := atan2(dx, dy) * 180 / pi]
data_first100[, bearing_absolute := (bearing_absolute + 360) %% 360]

# Vezmeme první bod každé události + směrový vektor z prvního segmentu
arrow_data_abs <- data_first100[
  , .SD[1],
  by = .(fishid, navigation_event_id, cluster_id)
]

vec_length <- 100
arrow_data_abs[, bearing_rad := bearing_absolute * pi / 180]
arrow_data_abs[, `:=`(
  x_end = easting + vec_length * sin(bearing_rad),  # sin, cos přehozeny kvůli směru os
  y_end = northing + vec_length * cos(bearing_rad)
)]

ggplot() +
  geom_sf(data = shape.rimov, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_segment(data = arrow_data_abs[cluster_id == 1],
               aes(x = easting, y = northing,
                   xend = x_end, yend = y_end,
                   color = fishid),
               arrow = arrow(length = unit(0.15, "cm")),
               linewidth = 0.8) +
  facet_wrap(~fishid)+
  coord_sf(xlim = c(461200, 462800), ylim = c(5410000, 5411000), expand = FALSE) +
  labs(
    title = "True movement direction (first 100 m)",
    x = "Easting", y = "Northing", color = "Fish ID"
  ) +
  theme_minimal(base_size = 13)


# 5. Vykreslení mapy se šipkami a polygonem nádrže
ggplot() +
  geom_sf(data = shape.rimov, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_segment(data = arrow_data[cluster_id ==1],
               aes(x = easting, y = northing,
                   xend = x_end, yend = y_end,
                   color = fishid),
               arrow = arrow(length = unit(0.15, "cm")),
               linewidth = 0.8) +
  coord_sf(xlim = c(461200, 462800), ylim = c(5410000, 5411000), expand = FALSE) +
  labs(
    title = "Initial 100 m movement vectors from home site",
    x = "Easting (m)", y = "Northing (m)", color = "Fish ID"
  ) +
  theme_minimal(base_size = 13)


# Fig. 3 # heat map ####
# Vytvoř si kopii dat pro práci
return_points <- copy(last_seg_sub[ !is.na(easting) & !is.na(northing) & cluster_id == 1])
return_points[fishid == "T449215_1", fishid_new := "#2"]
return_points[fishid == "T449313_1", fishid_new := "#1"]
return_points[, fishid_new := factor(fishid_new, levels = c("#1", "#2"))]
# Urči velikost gridu v metrech
grid_size <- 10

# Zaokrouhlení souřadnic na grid
return_points[, easting_bin := floor(easting / grid_size) * grid_size]
return_points[, northing_bin := floor(northing / grid_size) * grid_size]

# Spočítej počet výskytů v každé buňce
raster_counts <- return_points[, .N, by = .(easting_bin, northing_bin, fishid_new, transition_phase)]
raster_counts[, transition_phase := factor(transition_phase, levels = c("before_first", "after_last"))]
library("ggspatial")

sleep_spot_cluster1 <-  return_points[,.(x = mean(x_sleep), y= mean(y_sleep))]
ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  
  geom_raster(data = raster_counts, aes(x = easting_bin, y = northing_bin, fill = N)) +
  
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y),
             shape = 21, fill = "white", color = "black", size = 3) +
  
  #scale_fill_viridis_c(option = "magma", direction = -1, name = "Number of locations") +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt", name = "Number of locations")+
  coord_sf(xlim = c(461200, 462800), ylim = c(5408500, 5411000), expand = FALSE) +
  
  annotation_scale(location = "br", width_hint = 0.3) +
  
  facet_wrap(~fishid_new + transition_phase) +
  
  labs(
    title = "Heatmap of fish trajectories before and after foraging",
    x = "Easting [m]", y = "Northing [m]"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 10),
    legend.position = "right"
  )

# Parametr: velikost výřezu
buffer <- 1000
library(ggplot2)
library(patchwork)
library(ggplot2)
library(patchwork)
library(cowplot)  # pro get_legend

# Fixní rozsah pro výřezy
fixed_xlim <- c(461800, 462800)
fixed_ylim <- c(5410000, 5410950)

# Nové názvy pro fáze
phase_labels <- c(
  "before_first" = "Outbound",
  "after_last" = "Return"
)

# RP souřadnice – společné pro všechny výřezy
x0 <- sleep_spot_cluster1$x[1]
y0 <- sleep_spot_cluster1$y[1]

# Funkce pro výřez (bez legendy)
make_zoom_fixed <- function(fish_id, phase, x0, y0, xlim, ylim) {
  df <- raster_counts[fishid_new == fish_id & transition_phase == phase]
  df_zoom <- df[
    easting_bin >= xlim[1] & easting_bin <= xlim[2] &
      northing_bin >= ylim[1] & northing_bin <= ylim[2]
  ]
  
  ggplot() +
    geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
    geom_tile(data = df_zoom, aes(x = easting_bin, y = northing_bin, fill = N)) +
    geom_point(aes(x = x0, y = y0), shape = 21, fill = "white", color = "black", size = 2.5) +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 34),
                         name = "N loc.") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = paste(fish_id, phase)) +
    theme_minimal(base_size = 9) +
    theme(
      text= element_text(size = 15),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6, linetype = "dashed"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      legend.position = "none"
    )
}

# Výřezy pro každou rybu a fázi
plots <- list(
  bf5 = make_zoom_fixed("#2", "before_first", x0, y0, xlim = fixed_xlim, ylim = fixed_ylim),
  af5 = make_zoom_fixed("#2", "after_last",  x0, y0, xlim = fixed_xlim, ylim = fixed_ylim),
  bf3 = make_zoom_fixed("#1", "before_first", x0, y0, xlim = fixed_xlim, ylim = fixed_ylim),
  af3 = make_zoom_fixed("#1", "after_last",  x0, y0, xlim = fixed_xlim, ylim = fixed_ylim)
)

# Hlavní mapa uprostřed s legendou
main_plot <- ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  geom_tile(data = raster_counts, aes(x = easting_bin, y = northing_bin, fill = N)) +
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y),
             shape = 21, fill = "white", color = "black", size = 3) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 34),
                       name = "N Loc.") +
  annotate("rect",
           xmin = fixed_xlim[1], xmax = fixed_xlim[2],
           ymin = fixed_ylim[1], ymax = fixed_ylim[2],
           color = "grey70", fill = NA, linewidth = 0.6, linetype = "dashed")+
  facet_grid(fishid_new ~ transition_phase,
             labeller = labeller(
               transition_phase = as_labeller(phase_labels)
             )) +
  coord_sf(xlim = c(461200, 462800), ylim = c(5408700, 5411000), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  )

# Definuj pozice jednotlivých ploch
area_bf5  <- area(t = 1, l = 1, b = 1, r = 1)
area_af5  <- area(t = 1, l = 3, b = 1, r = 3)
area_bf3  <- area(t = 2, l = 1, b = 2, r = 1)
area_af3  <- area(t = 2, l = 3, b = 2, r = 3)
area_main <- area(t = 1, l = 2, b = 2, r = 2)  # hlavní mapa přes dva řádky

# Sestavení výsledného layoutu
final_plot <- wrap_plots(
  bf5 = plots$bf5,
  main = main_plot,  # ten už má legendu
  af5 = plots$af5,
  bf3 = plots$bf3,
  af3 = plots$af3,
  design = c(area_bf3, area_main, area_af3, area_bf5, area_af5)
)
# Výstup
print(final_plot)
# Uložení
# ggsave("central_layout_with_zoom_panels.pdf", layout_plot, width = 12, height = 8)

# rozdily vyuziti prostoru mezi odchodem a navratem ####
# ⬇️ Parametry
# Nastav velikost gridu

last_seg_sub[fishid == "T449215_1", fishid_new := "#2"]
last_seg_sub[fishid == "T449313_1", fishid_new := "#1"]

bin_size <- 50  # příklad
# Rasterizace
raster_counts_comp <- last_seg_sub %>%
  filter(cluster_id == 1)%>%
  mutate(
    easting_bin = floor(easting / bin_size) * bin_size,
    northing_bin = floor(northing / bin_size) * bin_size
  ) %>%
  group_by(fishid_new, transition_phase, easting_bin, northing_bin) %>%
  summarise(N = n(), .groups = "drop") %>%
  data.table()

# Získej fishid_new
fish_ids <- unique(raster_counts_comp$fishid_new)

# Výpočet overlapu a korelací
all_merged <- data.table()

for (fid in fish_ids) {
  before_data <- raster_counts_comp[fishid_new == fid & transition_phase == "before_first"]
  after_data  <- raster_counts_comp[fishid_new == fid & transition_phase == "after_last"]
  
  merged <- merge(
    before_data[, .(easting_bin, northing_bin, before_N = N)],
    after_data[, .(easting_bin, northing_bin, after_N = N)],
    by = c("easting_bin", "northing_bin"),
    all = TRUE
  )
  
  merged[is.na(before_N), before_N := 0]
  merged[is.na(after_N), after_N := 0]
  merged[, fishid := fid]
  
  all_merged <- rbind(all_merged, merged)
}

# Korelace
cor_summary <- all_merged[, .(
  rho = round(cor(before_N, after_N, method = "spearman"), 3)
), by = fishid]

# Spoj zpět
all_merged <- merge(all_merged, cor_summary, by = "fishid")

# Graf
ggplot(all_merged, aes(x = before_N, y = after_N, color = fishid)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~fishid) +
  labs(
    title = "Raster overlap between phases",
    subtitle = "Before foraging vs Return trajectory density",
    x = "Before foraging (counts per bin)",
    y = "Return (counts per bin)",
    color = "Fish ID"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  geom_text(
    data = cor_summary,
    aes(x = Inf, y = Inf, label = paste0("ρ = ", rho)),
    hjust = 1.1, vjust = 1.5,
    inherit.aes = FALSE,
    size = 4
  )
# rozdil pouze v oblasti blizko RP 
# Výpočet rozdílu: spojíš před a po
# Fixní rozsah
# 1. Nastav souřadnicový výřez (oblast kolem RP)
fixed_xlim <- c(461800, 462800)
fixed_ylim <- c(5410000, 5410950)

# 2. Vytvoř diff_data: rozdíl mezi after a before (už by měl existovat, ale pro jistotu):
diff_data <- merge(
  raster_counts_comp[transition_phase == "before_first", .(easting_bin, northing_bin, before_N = N, fishid_new)],
  raster_counts_comp[transition_phase == "after_last",  .(easting_bin, northing_bin, after_N = N, fishid_new)],
  by = c("easting_bin", "northing_bin", "fishid_new"),
  all = TRUE
)
diff_data[is.na(before_N), before_N := 0]
diff_data[is.na(after_N),  after_N  := 0]

# 3. Výřez do okolí RP
diff_zoom <- diff_data[
  easting_bin >= fixed_xlim[1] & easting_bin <= fixed_xlim[2] &
    northing_bin >= fixed_ylim[1] & northing_bin <= fixed_ylim[2]
]

# 4. Spočítej celkový počet detekcí ve výřezu podle ryby a fáze
totals_zoom <- diff_zoom[, .(
  total_before = sum(before_N),
  total_after  = sum(after_N)
), by = fishid_new]

# 5. Spoj zpět do diff_zoom a vypočítej hustoty
diff_zoom <- merge(diff_zoom, totals_zoom, by = "fishid_new")
diff_zoom[, `:=`(
  density_before = before_N / total_before,
  density_after  = after_N  / total_after
)]

# 6. Spočítej Spearmanovu korelaci mezi hustotami
rho_density_zoom <- diff_zoom[, .(
  rho_density = round(cor(density_before, density_after, method = "spearman"), 3)
), by = fishid_new]

# 7. Výstup
print(rho_density_zoom)

# Výpočet rozdílu hustoty
diff_zoom[, diff_density := density_after - density_before]
max_abs_diff <- max(abs(diff_zoom$diff_density), na.rm = TRUE)
# Mapa rozdílu hustoty
# Plot
density_overlap_gg <- ggplot() +
  geom_tile(data = diff_zoom, aes(x = easting_bin, y = northing_bin, fill = diff_density)) +
  geom_sf(data = shape.rimov, fill = NA, color = "black", size = 0.3) +  
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y), shape = 21, fill = "white", color = "black", size = 3)+
 # geom_point(aes(x = x0, y = y0), shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~fishid_new, ncol = 2) +
  coord_sf(xlim = fixed_xlim, ylim =fixed_ylim, expand = FALSE) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_diff, max_abs_diff),
    name = expression(atop(Delta~"density", (ret. - out.)))
  )+
  labs(
    title = "Difference in relative bin usage (density)",
    subtitle = "Positive = more used during return, Negative = more used during outbound movement",
    x = "Easting [m]", y = "Northing [m]"
  ) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")
  )+
  
  # >>> zde přidáš severku <<<
  annotation_north_arrow(
    location = "br",                # umístění: "tl", "tr", "bl", "br"
    which_north = "true",           # klasický sever
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_nautical 
  )
density_overlap_gg 

density_overlap_gg_f <- density_overlap_gg +
  plot_annotation(
    caption = "Positive = more used during return, Negative = more used during outbound movement",
    theme = theme(
      plot.caption = element_text(size = 15, hjust = 0.5)
    )
  )

heatmap_gg <- ggarrange(final_plot, density_overlap_gg_f, nrow = 2, heights = c(0.65, 0.35), labels = c("a)", "b)"),font.label = list(size = 15, color = "black", face = "plain", family = NULL))
heatmap_gg
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/heatmap_gg.png",
  plot = heatmap_gg,
  width = 10,      # šířka v palcích
  height = 10,     # výška v palcích
  dpi = 300
)

library(rKIN)
# selection of data 
  # #5
last_seg_sub[shore_50_merged == 1, shore_50_merged_name := "inshore" ]
last_seg_sub[shore_50_merged == 0, shore_50_merged_name := "offshore" ]

dt_5bf_dt <-  last_seg_sub[fishid == "T449215_1" & transition_phase == "before_first" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_5_bf <- dt_5bf_dt[, .(X = bottom_depth, Y = -depth_interp, species = shore_50_merged_name)]
  
dt_5al_dt <-  last_seg_sub[fishid == "T449215_1" & transition_phase == "after_last" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_5_al <- dt_5al_dt[, .(X = bottom_depth, Y = -depth_interp, species = shore_50_merged_name)]

dt_3bf_dt <-  last_seg_sub[fishid == "T449313_1" & transition_phase == "before_first" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_3_bf <- dt_3bf_dt[, .(X = bottom_depth, Y = -depth_interp, species = shore_50_merged_name)]

dt_3al_dt <-  last_seg_sub[fishid == "T449313_1" & transition_phase == "after_last" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_3_al <- dt_3al_dt[, .(X = bottom_depth, Y = -depth_interp, species = shore_50_merged_name)]

dt_bf_dt <-  last_seg_sub[transition_phase == "before_first" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_bf <- dt_bf_dt[, .(X = bottom_depth, Y = -depth_interp, species = fishid)]

dt_al_dt <-  last_seg_sub[transition_phase == "after_last" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_al <- dt_al_dt[, .(X = bottom_depth, Y = -depth_interp, species = fishid)]

dt_3_dt <-  last_seg_sub[fishid == "T449313_1" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_3 <- dt_3_dt[, .(X = bottom_depth, Y = -depth_interp, species = transition_phase)]

dt_5_dt <-  last_seg_sub[fishid == "T449215_1" & !is.na(shore_50_merged)& !is.na(bottom_depth)& !is.na(depth_interp)]
rkin_5 <- dt_5_dt[, .(X = bottom_depth, Y = -depth_interp, species = transition_phase)]

hval_set <- c(2,2)

# calculation of overlay
rkin_5_bf_kin<- estKIN(data=rkin_5_bf, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)
rkin_5_al_kin<- estKIN(data=rkin_5_al, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)
rkin_3_bf_kin<- estKIN(data=rkin_3_bf, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)
rkin_3_al_kin<- estKIN(data=rkin_3_al, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)

rkin_bf_kin<- estKIN(data=rkin_bf, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)
rkin_al_kin<- estKIN(data=rkin_al, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)

rkin_3_kin<- estKIN(data=rkin_3, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)
rkin_5_kin<- estKIN(data=rkin_5, x="X", y="Y", group="species", levels=c(25,50,95), hval = hval_set)

# Extract the area of each polygon
rkin_5_bf_line <- fortify(rkin_5_bf_kin$estObj)
rkin_5_al_line <- fortify(rkin_5_al_kin$estObj)
rkin_3_bf_line <- fortify(rkin_3_bf_kin$estObj)
rkin_3_al_line <- fortify(rkin_3_al_kin$estObj)

rkin_bf_line <- fortify(rkin_bf_kin$estObj)
rkin_al_line <- fortify(rkin_al_kin$estObj)
rkin_3_line <- fortify(rkin_3_kin$estObj)
rkin_5_line <- fortify(rkin_5_kin$estObj)

inshore_col <- "red"
offshore_col <- "blue"
pol_col <- "grey"
alpha_vec <- c(0.6,0.3,0.1,0.6,0.3,0.1 )
ratio_size <- 4

rkin_5_bf_gg <-    
  rkin_5_bf_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("c) #2 - Outbound")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  guides(fill = guide_legend(nrow = 2))+
  theme_bw()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    #plot.title = element_text(hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_blank()
  )

rkin_5_bf_gg 

rkin_5_al_gg <-    
  rkin_5_al_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
    geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
    geom_sf(aes(fill=  Group_ConfInt))+
    ggtitle("d) #2 - Return")+
    scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
    labs(x = "bottom depth (m)", y = "depth (m)")+
    scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                      na.value = "grey")+
    #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
    geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
    coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
    #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
    #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
    theme_bw()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    #plot.title = element_text(hjust = 1),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_blank()
  )

rkin_5_al_gg 

rkin_3_bf_gg <-    
  rkin_3_bf_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("a) #1 - Outbound")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    #plot.title = element_text(hjust = 1),
    legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
rkin_3_bf_gg 

rkin_3_al_gg <-    
  rkin_3_al_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  #ggtitle("Most")+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("b) #1 - Return")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    #plot.title = element_text(hjust = 1),
    legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
rkin_3_al_gg 
  
rkin_gg <- (rkin_3_bf_gg + rkin_3_al_gg )/(rkin_5_bf_gg + rkin_5_al_gg )
rkin_gg

# rkin picture ####
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/rkin_depth_gg.jpg",  # nebo .jpg
  plot = rkin_gg,
  width = 8,    # šířka v palcích
  height = 8,    # výška v palcích
  dpi = 300      # rozlišení pro tisk
)

rkin_bf_gg <-    
  rkin_bf_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  #ggtitle("Most")+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("#3 Before forage")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.y = element_blank()
  )
rkin_bf_gg 

rkin_al_gg <-    
  rkin_al_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  #ggtitle("Most")+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("#3 Return")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 1),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
rkin_al_gg 

rkin_3_gg <-    
  rkin_3_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  #ggtitle("Most")+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("#3 Before forage")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.y = element_blank()
  )
rkin_3_gg 

rkin_5_gg <-    
  rkin_5_line %>%  
  mutate(Group_ConfInt = paste(Group,"_",ConfInt, "%", sep = ""))%>%
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  #ggtitle("Most")+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill=  Group_ConfInt))+
  ggtitle("#3 Return")+
  scale_y_discrete(labels=c("30" = "30", "-25" = "25", "-20" = "20", "-15" = "15", "-10" = "10", "-5" = "5", "0" = "0"))+
  labs(x = "bottom depth (m)", y = "depth (m)")+
  scale_fill_manual(values = alpha(c( offshore_col,  offshore_col, offshore_col, inshore_col,inshore_col,inshore_col),alpha_vec), 
                    na.value = "grey")+
  #geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "black", linewidth  = 1, linetype = "dashed")+
  geom_polygon(data = data.table(x = c(0,0,35,0),y = c(0,-35,-35,0)), aes(x = x, y= y), fill = pol_col, col = "black")+
  coord_sf(ylim = c(2,-33.2), x= c(1.9, 40))+
  #geom_sf_text( x=5, y= 2, label= overlap_tab[OverlapID == "pike_95" & lake == "Rimov",]$wels_95, size = ratio_size, col = pike_col )+
  #geom_sf_text( x=35, y= 2, label= overlap_tab[OverlapID == "wels_95" & lake == "Rimov",]$pike_95, size = ratio_size, col = wels_col )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 1),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
rkin_5_gg 

# miss point ans surise[sunset]  ####
# jestli vzdalenost od RP je zavisla na to jestli byla pred nebo po 

miss_offshore <- merge(return_accuracy_cl1, last_offshore_tab, by =c("navigation_event_id"))
unique(miss_offshore$move_type)
miss_offshore[move_type == "only_shore"]
miss_offshore[move_type == "only offshore" ]
miss_offshore[move_type == "shore to offshore"]
miss_offshore[move_type == "shore to offshore", dist_to_rp := 0]
ggplot(miss_offshore[move_type != "only_shore"], aes( y = dist_to_rp, x = dist_m, col  = fishid_new))+ geom_point()+facet_wrap(~fishid_new)
# předpokládám, že miss_m = dist_to_rp
# a vzdálenost trasy = dist_m

miss_offshore[, angle_rad := atan(dist_to_rp / dist_m)]
miss_offshore[, angle_deg := atan(dist_to_rp / dist_m) * (180 / pi)]
# převod na stupně

# mean and CI
calc_circ <- function(a){
  m <- mean.circular(circular(a, units="degrees"))
  boot <- replicate(5000,
                    mean.circular(circular(sample(a, replace=TRUE), units="degrees")))
  ci <- quantile(boot, c(0.025,0.975))
  list(mean=as.numeric(m), lci=ci[1], uci=ci[2])
}

miss_offshore[move_type=="offshore to shore",
              calc_circ(angle_deg),
              by=fishid_new]
###

ggplot(miss_offshore[move_type == "offshore to shore" & dist_m > 0],
       aes(x = angle_deg, fill = fishid_new)) +
  geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
  labs(x = "Angular error (degrees)",
       y = "Count",
       fill = "Individual") +
  theme_minimal(base_size = 14)

ggplot(miss_offshore[move_type == "offshore to shore" & dist_m > 0],
       aes(y = angle_deg, x = dist_m, color = fishid_new)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(y = "Angular error (°)",
       x = "Return distance (m)",
       color = "Individual") +
  theme_minimal(base_size = 14)

miss_gg_dist_length <- ggplot(miss_offshore[move_type == "offshore to shore" & dist_m > 0],
       aes(y = dist_to_rp , x = dist_m, color = fishid_new)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(y = "Lenght of the last offhore segment (m)",
       x = "Distance to RP (m)",
       color = "Individual") +
  theme_bw(base_size = 15)
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/miss_gg_dist_lengt.jpg",  # nebo .jpg
  plot = miss_gg_dist_length,
  width = 7,    # šířka v palcích
  height = 5,    # výška v palcích
  dpi = 300      # rozlišení pro tisk
)

ggplot(miss_offshore[move_type == "offshore to shore" & dist_m > 0],
       aes( x = dist_m, color = fishid_new)) +
  geom_density(size = 2, alpha = 0.8) +
  labs(y = "Density",
       x = "Return distance (m)",
       color = "Individual") +
  theme_minimal(base_size = 14)

accuracy_light <- merge(
  return_accuracy_cl1,
  light_summary[, .(fishid, navigation_event_id, arrival_before_sunrise, arrival_sunrise)],
  by = c("fishid", "navigation_event_id"),
  all.x = TRUE
)

accuracy_light[, .(r = cor(arrival_sunrise, dist_to_rp, use = "complete.obs")), by = fishid]

accuracy_light[arrival_before_sunrise == T,arrival_before_sunrise_new := "Before sunset"]
accuracy_light[arrival_before_sunrise == F,arrival_before_sunrise_new := "After sunset"]
accuracy_light[,arrival_before_sunrise_new := factor(arrival_before_sunrise_new, levels = c("Before sunset", "After sunset"))]
accuracy_light[fishid == "T449313_1", fishid_new := "#1"]
accuracy_light[fishid == "T449215_1", fishid_new := "#2"]
accuracy_light[, fishid_new := factor(fishid_new, levels = c("#1", "#2"))]
accuracy_light[move_type == "shore to offshore", dist_to_rp := 0]

ggplot(accuracy_light[!is.na(arrival_sunrise)], 
       aes(x = arrival_sunrise, y = dist_to_rp)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~fishid)+
  labs(
    x = "Minutes relative to sunrise (negative = before)",
    y = "Distance from RP (m)",
    title = "Return accuracy vs. timing from sunrise"
  ) +
  theme_minimal(base_size = 13)
# miss point picture sunset x sunrise + map ####
nav_return_sf_cl2 <- merge(nav_return_sf_cl1, accuracy_light[!is.na(arrival_before_sunrise) & move_type != "only shore",.(navigation_event_id, move_type, arrival_before_sunrise_new)], by = c("navigation_event_id"))
med_miss_point <- accuracy_light[!is.na(arrival_before_sunrise) & move_type != "only shore", .(med_dist = median(dist_to_rp), na.rm = T), by =.(fishid_new, arrival_before_sunrise_new)] 
med_miss_point
miss_point_gg <- ggplot(accuracy_light[!is.na(arrival_before_sunrise) & move_type != "only shore"], 
       aes(x = arrival_before_sunrise_new, y = dist_to_rp, fill = fishid_new)) +
  geom_violin(trim = FALSE, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5, aes(color = fishid)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c( "#00BFC4","#F8766D")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  facet_wrap(~fishid_new)+
  
  labs(
    x = "Returned before sunrise",
    y = "Distance from RP (m)",
    title = "Return accuracy vs. timing relative to sunrise"
  ) +
  theme_minimal(base_size = 15)+
  theme(
    text = element_text(size = 15),
    plot.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )
miss_point_gg
# map
nav_return_sf_cl2 <- nav_return_sf_cl2 %>%
  mutate(
    fishid_new = case_when(
      fishid == "T449215_1" ~ "#2",
      fishid == "T449313_1" ~ "#1",
      TRUE ~ fishid
    ),
    fishid_new = factor(fishid_new, levels = c("#1", "#2"))
  )

miss_point_map_gg <- ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  geom_sf(data = nav_return_sf_cl2, aes(color = dist_to_rp), size = 2) +
  geom_point(data = sleep_spot_cluster1, aes(x = x, y = y),
             shape = 21, fill = "white", color = "black", size = 3) +
  scale_color_viridis_c(option = "plasma", name = "Dist to RP (m)") +
  labs(title = "Points from which fish navigated around srore to RP") +
  facet_wrap(~fishid_new)+
  coord_sf(xlim = c(461200, 462800), ylim = c(5410200, 5411100), expand = FALSE) +
  theme_minimal()+
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.title = element_text(size = 15),
    legend.position = "right"
  )

miss_point_gg_final <- ggarrange(miss_point_gg, miss_point_map_gg, ncol = 1, nrow = 2, labels = c("a)", "b)"),font.label = list(size = 15, color = "black", face = "plain", family = NULL))
miss_point_gg_final
ggsave(
  filename = "~/Teri/pikeperch_navigation/output/miss_point_gg.jpg",  # nebo .jpg
  plot = miss_point_gg_final,
  width = 8,    # šířka v palcích
  height = 8,    # výška v palcích
  dpi = 300      # rozlišení pro tisk
)

ggplot(accuracy_light[!is.na(arrival_before_sunrise)], 
       aes(x = move_type, y = dist_to_rp, fill = arrival_before_sunrise)) +
  geom_violin(scale = "width", trim = FALSE, color = "black",drop = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.3, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(width = 0.8), 
              shape = 21, size = 2, fill = "red", color = "black") +
  facet_wrap(~fishid)+
  labs(
    x = "Move type",
    y = "Distance from RP at return (m)",
    fill = "Before sunrise",
    title = "Return accuracy by move type and sunrise timing"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

accuracy_light[!is.na(arrival_before_sunrise), 
               .N, 
               by = .(move_type, arrival_before_sunrise,fishid)][
                 order(move_type, -arrival_before_sunrise,fishid)
               ]

# Agregace dat
count_df <- accuracy_light[!is.na(arrival_before_sunrise), 
                           .N, 
                           by = .(move_type, arrival_before_sunrise,fishid)
]

count_sum <- accuracy_light[, 
                           .N, 
                           by = .(move_type, fishid)
]
count_percent <- count_sum[
  , .(N = sum(N)), by = .(fishid, move_type)
][
  , .(move_type, percent = N / sum(N) * 100), by = fishid
]
# Plot
ggplot(count_df, aes(x = move_type, y = N, fill = arrival_before_sunrise)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    x = "Move type",
    y = "Number of return events",
    fill = "Before sunrise",
    title = "Count of return types by sunrise timing"
  ) +
  facet_wrap(~fishid)+
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(count_percent, aes(x = move_type, y = percent)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", fill = "steelblue") +
  labs(
    x = "Move type",
    y = "Proportion of return events (%)",
    title = "Return move types as percentage per fish"
  ) +
  facet_wrap(~fishid) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ted kategorie departure 
departure_data <- merge(
  movement_summary[transition_phase == "before_first"],
  light_summary[, .(fishid, navigation_event_id, departure_after_sunset)],
  by = c("fishid", "navigation_event_id"),
  all.x = TRUE
)

departure_counts <- departure_data[!is.na(departure_after_sunset),
                                   .N,
                                   by = .(move_type, departure_after_sunset)
][order(move_type, -departure_after_sunset)]

ggplot(departure_counts, aes(x = move_type, y = N, fill = departure_after_sunset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    x = "Move type (departure)",
    y = "Number of events",
    fill = "After sunset",
    title = "Departure move types vs. sunset timing"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

departure_data <- merge(
  movement_summary[transition_phase == "before_first"],
  light_summary[, .(fishid, navigation_event_id, departure_sunset)],
  by = c("fishid", "navigation_event_id"),
  all.x = TRUE
)

ggplot(departure_data[!is.na(departure_sunset)], 
       aes(x = move_type, y = departure_sunset)) +
  geom_violin(fill = "lightgray", scale = "width", trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  labs(
    x = "Move type (departure)",
    y = "Minutes relative to sunset (negative = before)",
    title = "Departure timing relative to sunset by move type"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# start and end of forage ####
last_seg_sub 

# předpoklad: data.table last_seg_sub, čas v `timestamp`
library(data.table)

# seřadit, aby "první/poslední" dávaly smysl
setorder(last_seg_sub, navigation_event_id, timestamp)

# a) BEFORE_FIRST: pro každé navigation_event_id 1. a poslední řádek
bef <- last_seg_sub[transition_phase == "before_first"
][, .SD[unique(c(1L, .N))], by = navigation_event_id]

# b) AFTER_LAST: pro každé navigation_event_id jen 1. řádek
aft <- last_seg_sub[transition_phase == "after_last"
][, .SD[1L], by = navigation_event_id]

# spojit
sel <- rbindlist(list(bef, aft), use.names = TRUE)
sel

ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  geom_point(data = sel[cluster_id == 1 &  fishid == "T449313_1"], aes(x = easting, y = northing, col = transition_phase), size = 2) +
  geom_point(aes(x = x0, y = y0), shape = 21, fill = "white", color = "black", size = 2.5) +
  #scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 34),
   #                    name = "N loc.") +
  coord_sf(xlim = c(range(sel$easting)), ylim = c(range(sel$northing)), expand = FALSE) +
  #labs(title = paste(fish_id, phase)) +
  theme_minimal(base_size = 9) +
  theme(
    text= element_text(size = 15),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6, linetype = "dashed"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggplot() +
  geom_sf(data = shape.rimov, fill = "darkgrey", color = "black", linewidth = 0.3) +
  geom_point(data = sel[cluster_id == 1 &  fishid == "T449215_1"], aes(x = easting, y = northing, col = transition_phase), size = 2) +
  geom_point(aes(x = x0, y = y0), shape = 21, fill = "white", color = "black", size = 2.5) +
  #scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 34),
  #                    name = "N loc.") +
  coord_sf(xlim = c(range(sel$easting)), ylim = c(range(sel$northing)), expand = FALSE) +
  #labs(title = paste(fish_id, phase)) +
  theme_minimal(base_size = 9) +
  theme(
    text= element_text(size = 15),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6, linetype = "dashed"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )
