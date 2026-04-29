##############################################
# Segmenatation of trajetories 
#load packages
library(data.table)
library(RPostgreSQL)
library(sp) 
library(lubridate)
library(geosphere)
library(rgdal)
library(raster)
library(ggplot2)
library(rpostgis)
library(sf)
library(terra)
library(suncalc)
library(tidyverse)
library(zoo)
library(gdistance)

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


# upload positions from database
positions <- positions[ up_validpos==TRUE & up_gamres < 75, ]
positions[, date := as.Date(up_timestamp_utc)]

# calculation of swimming distance
comp.dist <- function (easting, northing, depth, timestamp)
  {
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
positions[, week := week(up_timestamp_utc )]
positions[, tag_date := paste(fi_fishid, date, sep = "_")]


# getting shape of the Rimov
shape.rimov  <- st_read("~/Teri/shp_files/rimov_469/Rimov_pol_469m_UTM33.shp") 

# Definuj polohu Rimova
lat <- 48.85
lon <- 14.49

# Definuj datum/daty
dates <- seq(as.Date("2017-06-01"), as.Date("2017-12-07"), by = "day")

# Spočítej časy
sun_data <- data.table(getSunlightTimes(date = dates, lat = lat, lon = lon, keep = c("sunrise", "sunset")))

# light data 
light_int <- data.table(read_csv("~/Teri/data/Teri_Vor_Rimov_envi.csv"))
light_int[, timestamp := dmy_hms(Timestamp)]
light_int[, timestamp := as.POSIXct(round(as.numeric(timestamp)), origin = "1970-01-01", tz = "UTC")]

# Předpokládáme, že tvoje tabulka se jmenuje dt a timestamp je ve formátu POSIXct
setkey(light_int, timestamp)

# Vytvoříme nový vektor všech minut v rozsahu dat
minutes_all <- data.table(timestamp = seq(min(light_int$timestamp), max(light_int$timestamp), by = "1 sec"))

# Sloučíme s původní tabulkou
dt_full <- merge(minutes_all, light_int[, .(timestamp, light_up)], by = "timestamp", all.x = TRUE)

# Interpolace pomocí naapproxfun (příp. zoo::na.approx)
dt_full[, light_up_interp := approx(x = timestamp, y = light_up, xout = timestamp, rule = 2)$y]

# kalkulace vzdalenosti od brehu
# Převedeme na sf body
positions_sf <- st_as_sf(positions, coords = c("up_easting", "up_northing"), crs = st_crs(shape.rimov))
# Výpočet vzdáleností (v metrech, pokud CRS je metrický)


# select only days with detected resting place
daily_sleep_locations <- data.table(read_csv("~/Teri/pikeperch_navigation/data/daily_sleep_locations.csv"))
daily_sleep_locations[, tag_date := paste(fi_fishid, date, sep = "_")]
unique(daily_sleep_locations$fi_fishid)

daily_sleep_locations[
  fi_fishid == "T449215_1" & 
    date %in% as.Date(c("2017-07-17", "2017-07-18", "2017-07-19", 
                        "2017-07-20", "2017-07-21", "2017-07-24", 
                        "2017-07-31", "2017-08-03", "2017-08-04"))
]

# take only dates and fish that we detected resting place 
pos_sub <- merge(positions, daily_sleep_locations[,.(x_sleep, y_sleep, tag_date)], by = c("tag_date"))

# Předpoklad: máš nav_last_seg jako data.table
pos_sub[, up_timestamp_utc := as.POSIXct(round(as.numeric(up_timestamp_utc)), origin = "1970-01-01", tz = "UTC")]
pos_sub <- merge(pos_sub, dt_full[,.(timestamp, light_up_interp)], by.x = c("up_timestamp_utc"), by.y = c("timestamp"), all.x = T)
pos_sub[, light_lux := 120*light_up_interp]

data <- pos_sub[up_validpos==TRUE & up_gamres < 50,.(fishid = fi_fishid, up_gamres, up_validpos, timestamp = up_timestamp_utc,easting = up_easting,  northing = up_northing, bottom_depth = up_bottom_depth, distfromdam, egam = up_egam,ngam =  up_ngam, depth = up_depth, x_sleep, y_sleep, light_lux)]

# Výpočet vzdálenosti od domácí oblasti
data[, distance_to_home := sqrt((easting - x_sleep)^2 + (northing - y_sleep)^2)]

# home radius threshold 
home_radius <- 30

# 5. Určení jestli je doma
data[, is_home := distance_to_home <= home_radius]

#ggplot(data_clean_inlake, aes(x = easting, y =  northing, col = is_home))+geom_point()
# write data that we need to do everyrthing again 
#fwrite(data, "~/Teri/pikeperch_navigation/data/ready_positions.csv")

# 6. Identifikace návratových událostí
# 6. Detekce návratových událostí s podmínkou minimální délky mimo domov
min_departure_duration <- 20  # v minutách

data[, id := .I]  # přegenerovat indexy

library(data.table)

# Pokud není, převedeme na data.table
setDT(data)

# Pro jistotu: seřadit podle fishid a timestamp
setorder(data, fishid, timestamp)

# Předem si nastav minimální délku výletu (v minutách)
min_departure_duration <- 60

# Rozdělit data podle jedince
data_split <- split(data, by = "fishid")


navigation_summary_all <- rbindlist(lapply(data_split, function(subdata) {
  navigation_events <- list()
  leaving <- NULL
  
  for (i in 1:nrow(subdata)) {
    row <- subdata[i]
    
    if (!row$is_home && is.null(leaving)) {
      leaving <- list(
        start_idx = row$id,
        start_time = row$timestamp,
        start_point = c(row$easting, row$northing)
      )
      
    } else if (row$is_home && !is.null(leaving)) {
      duration_min <- as.numeric(difftime(row$timestamp, leaving$start_time, units = "mins"))
      
      if (duration_min >= min_departure_duration) {
        leaving$end_idx <- row$id
        leaving$end_time <- row$timestamp
        leaving$end_point <- c(row$easting, row$northing)
        
        navigation_events[[length(navigation_events) + 1]] <- leaving
      }
      leaving <- NULL
    }
  }
  
  if (length(navigation_events) == 0) return(NULL)
  
  # VLOŽÍME subdata do prostředí pomocí closure
  fishid_current <- subdata$fishid[1]
  
  rbindlist(lapply(navigation_events, function(ev) {
    segment <- subdata[id %between% c(ev$start_idx, ev$end_idx)]
    
    direct_distance <- sqrt((ev$start_point[1] - ev$end_point[1])^2 + 
                              (ev$start_point[2] - ev$end_point[2])^2)
    
    path_distance <- sum(sqrt(
      (segment$easting[-length(segment$easting)] - segment$easting[-1])^2 + 
        (segment$northing[-length(segment$northing)] - segment$northing[-1])^2
    ))
    
    straightness <- ifelse(path_distance > 0, direct_distance / path_distance, NA)
    
    duration_sec <- as.numeric(difftime(ev$end_time, ev$start_time, units="secs"))
    duration_min <- duration_sec / 60
    
    n_obs <- nrow(segment)
    ideal_n_obs <- ceiling(duration_sec / 15)
    coverage_ratio <- round(n_obs / ideal_n_obs, 3)
    
    data.table(
      fishid = fishid_current,
      start_time = ev$start_time,
      end_time = ev$end_time,
      duration_min = duration_min,
      direct_distance_m = direct_distance,
      path_length_m = path_distance,
      straightness_index = straightness,
      n_observations = n_obs,
      ideal_observations = ideal_n_obs,
      observation_coverage = coverage_ratio
    )
  }))
}))

navigation_summary_all[, date_start := as.Date(start_time)]
navigation_summary_all[, date_end := as.Date(end_time)]

# calculation of time from sunset/sunrise
navigation_summary_all <- merge(navigation_summary_all,  sun_data[,.(date_start = date,sunrise_sdate = sunrise, sunset_sdate = sunset )], by = c("date_start"))
navigation_summary_all <- merge(navigation_summary_all,  sun_data[,.(date_end = date,sunrise_edate = sunrise, sunset_edate = sunset )], by = c("date_end"))

navigation_summary_all[, time_from_sunrise := as.numeric(difftime(end_time,sunrise_edate, units = "mins"))]
navigation_summary_all[, time_from_sunset := as.numeric(difftime(start_time,sunset_sdate, units = "mins"))]


# adding index for each occasion
navigation_summary_all[, navigation_event_id := .I] 

# musime zjists zda spi na stejnem miste
setDT(daily_sleep_locations)
setDT(navigation_summary_all)

# Přidáme sloupec "date" pro porovnání
navigation_summary_all[, end_date := as.Date(end_time)]

# Spojíme navigační eventy s místem spánku, kde ryba spala v DEN návratu
navigation_with_sleep <- merge(navigation_summary_all, 
                               daily_sleep_locations[, .(fishid = fi_fishid, date, x_sleep, y_sleep)],
                               by.x = c("fishid", "end_date"),
                               by.y = c("fishid", "date"),
                               all.x = TRUE)

# Připravíme předešlá místa spánku pro každou rybu
daily_sleep_locations[, prev_date := date + 1]  # posuneme vpřed pro join na návrat
resting_previous <- daily_sleep_locations[, .(fishid = fi_fishid, prev_date, x_sleep_prev = x_sleep, y_sleep_prev = y_sleep)]

# Spojíme s `navigation_with_sleep` podle (fishid, end_date = prev_date)
navigation_with_sleep <- merge(navigation_with_sleep, resting_previous,
                               by.x = c("fishid", "end_date"),
                               by.y = c("fishid", "prev_date"),
                               all.x = TRUE)
# Vzdálenost mezi dnešním a včerejším resting place
navigation_with_sleep[, return_distance := sqrt((x_sleep - x_sleep_prev)^2 + (y_sleep - y_sleep_prev)^2)]


# zooming to standard
# tady vybirame jenom dobry subset kdyy mimimalni doba inaktivity je 3 hodin, maximalni 23 hodin 
# maximalnlni vzdalenost rest place je 50 m od predesle noci
# minimalni covarega je 0% - tj mame aspon 20% dat z te dane trajektorie
min_duration <- 60*3
max_duration <- 60*23
min_return_distance <- 50
min_observation_coverage <- 0

# finalni dataset eventu 
nav_sum_sub <- navigation_with_sleep[duration_min > min_duration & duration_min < max_duration &
                                       return_distance < min_return_distance & observation_coverage > min_observation_coverage]

n_events <- nav_sum_sub[, .(n_event = length(unique(navigation_event_id))), by = fishid]
n_events
fwrite(nav_sum_sub, "~/Teri/pikeperch_navigation/data/nav_sum_sub.csv")

nav_sum_sub <- data.table(read_csv("~/Teri/pikeperch_navigation/data/nav_sum_sub.csv"))
# testing       
navigation_summary_all[navigation_event_id == 10, ]
nav_sum_sub[fishid == "T449202_1"]
nav_sum_sub[fishid == "T449208_1" & date_start =="2017-07-04"]


data <- data.table(read_csv("~/Teri/pikeperch_navigation/data/ready_positions.csv"))

# --- Inicializace tabulky ---
# Smyčka pro rozdělení trajektorií na segmenty a klasifikaci jejich částí
# Každý navigační event je zpracován zvlášť podle fish_id a časového rozsahu
# Pro vybrané navigation_event_id lze specifikovat vlastní hodnotu bufferu (buffer_map), jinak se použije defaultní (default_buffer)
# Buffer určuje zónu kolem domovského místa (RP), kde nedetekujeme foraging eventy
# Výsledkem je klasifikace části trajektorie na foraging fáze (first, middle, last) a přechodové fáze (before_first, after_last, between)

# Před cyklem si připrav data_split:
data_split <- split(data, by = "fishid")

# OSTATNI JEDINCI
# pro ostatni pouzivame kombinacci speed_threshold <- 0.1 ;  no buffer  - default_buffer <- 0
# ostatni parametry jsou ty same

# Parametry
# Parametry
window_half <- 5
si_threshold <- 0.75
speed_threshold <- 0.2 

# --- Inicializace tabulky ---
navigation_segments <- data.table()

# --- Individuální buffery pro vybrané eventy ---
buffer_map <- c(
  "54" = 90,
  "514" = 80,
  "532" = 110
)
default_buffer <- 0 # original 250

for (i in seq_len(nrow(nav_sum_sub))) {
  
  fish_id <- nav_sum_sub$fishid[i]
  start_time <- nav_sum_sub$start_time[i]
  end_time <- nav_sum_sub$end_time[i]
  nav_event_id <- nav_sum_sub$navigation_event_id[i]
  
  data_fish <- data_split[[as.character(fish_id)]]
  
  segment <- data_fish[timestamp >= start_time & timestamp <= end_time]
  if (nrow(segment) < (window_half * 2 + 1)) next
  
  # Výpočet vzdálenosti a bufferu
  segment[, dist_home := sqrt((easting - x_sleep)^2 + (northing - y_sleep)^2)]
  
  this_buffer <- if (as.character(nav_event_id) %in% names(buffer_map)) {
    buffer_map[[as.character(nav_event_id)]]
  } else {
    default_buffer
  }
  segment[, near_home := dist_home < this_buffer]
  
  # Výpočet metrik
  segment[, `:=`(
    si_local = sapply(1:.N, function(j) {
      idx_start <- max(1, j - window_half)
      idx_end <- min(.N, j + window_half)
      if ((idx_end - idx_start) < 2) return(NA_real_)
      dx <- diff(egam[idx_start:idx_end])
      dy <- diff(ngam[idx_start:idx_end])
      d_direct <- sqrt((egam[idx_start] - egam[idx_end])^2 + (ngam[idx_start] - ngam[idx_end])^2)
      d_total <- sum(sqrt(dx^2 + dy^2))
      if (d_total == 0) return(NA_real_)
      return(d_direct / d_total)
    }),
    
    speed_local = sapply(1:.N, function(j) {
      idx_start <- max(1, j - window_half)
      idx_end <- min(.N, j + window_half)
      if (idx_end == idx_start) return(NA_real_)
      dt <- as.numeric(difftime(timestamp[idx_end], timestamp[idx_start], units = "secs"))
      d <- sqrt((egam[idx_start] - egam[idx_end])^2 + (ngam[idx_start] - ngam[idx_end])^2)
      if (dt == 0) return(NA_real_)
      return(d / dt)
    })
  )]
  
  # Detekce foraging eventu
  segment[, foraging_event := si_local < si_threshold | speed_local < speed_threshold]
  segment[near_home == TRUE, foraging_event := FALSE]
  segment[, foraging_group := rleid(foraging_event)]
  
  segment[, foraging_phase := "none"]
  segment[, transition_phase := "none"]
  
  forage_blocks <- segment[foraging_event == TRUE, .(
    start_idx = .I[1],
    end_idx = .I[.N]
  ), by = foraging_group]
  
  # --- Klasifikace fází ---
  if (nrow(forage_blocks) > 0) {
    first_block <- forage_blocks[1]
    last_block <- forage_blocks[.N]
    
    segment[first_block$start_idx:first_block$end_idx, foraging_phase := "first"]
    if (nrow(forage_blocks) > 2) {
      for (j in 2:(nrow(forage_blocks) - 1)) {
        blk <- forage_blocks[j]
        segment[blk$start_idx:blk$end_idx, foraging_phase := "middle"]
      }
    }
    if (nrow(forage_blocks) > 1) {
      segment[last_block$start_idx:last_block$end_idx, foraging_phase := "last"]
    }
    
    if (first_block$start_idx > 1) {
      segment[1:(first_block$start_idx - 1), transition_phase := "before_first"]
    }
    if (last_block$end_idx < nrow(segment)) {
      segment[(last_block$end_idx + 1):.N, transition_phase := "after_last"]
    }
    if (nrow(forage_blocks) > 1) {
      for (j in 1:(nrow(forage_blocks) - 1)) {
        blk1 <- forage_blocks[j]
        blk2 <- forage_blocks[j + 1]
        if ((blk2$start_idx - blk1$end_idx) > 1) {
          segment[(blk1$end_idx + 1):(blk2$start_idx - 1), transition_phase := "between"]
        }
      }
    }
  } else {
    # Pokud nejsou žádné foraging bloky, vše je přechod
    segment[, transition_phase := "whole"]
  }
  
  segment[, navigation_event_id := nav_event_id]
  navigation_segments <- rbind(navigation_segments, segment, fill = TRUE)
}

navigation_segments[, n_obs := .N, by = navigation_event_id]
navigation_segments[,  duration_sec := as.numeric(difftime(max(timestamp), min(timestamp),units = "secs")), by = navigation_event_id]
navigation_segments[, ideal_n_obs := as.numeric(ceiling(duration_sec / 15))]
navigation_segments[, coverage_ratio := round(n_obs / ideal_n_obs, 3)]

# only subset 
nav_last_seg_fin <- navigation_segments[navigation_event_id %in% unique(nav_sum_sub$navigation_event_id) & !is.na(x_sleep)]
unique(nav_last_seg_fin[is.na(x_sleep)]$fishid)

#fwrite(navigation_segments, "~/Teri/pikeperch_navigation/data/navigation_segments.csv")
fwrite(navigation_segments, "~/Teri/pikeperch_navigation/data/navigation_segments_ss01_bz0.csv")

