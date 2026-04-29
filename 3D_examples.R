library(terra)

# Load necessary libraries
# If you don't have them installed, run install.packages("terra"), install.packages("ggplot2"), install.packages("plotly")
library(terra)
library(ggplot2)
library(plotly)
library(tidyverse)
library(RPostgreSQL)
library(rpostgis)
library(data.table)
library(zoo)


data_3d <- data.table(read_csv("~/Teri/pikeperch_navigation/data/navigation_segments.csv"))
full_data <- data.table(read_csv("~/Teri/pikeperch_navigation/data/ready_positions.csv"))

nav_ev_id <- c(87, 72, 103, 127, 123,137,152)

full_data_merge <- merge(full_data, data_3d[,.(fishid, timestamp,transition_phase, navigation_event_id )], by = c("fishid", "timestamp"), all.x =T)

library(data.table)

# Předpoklad: full_data_merge je data.table
setkey(full_data_merge, fishid, timestamp)

setkey(full_data_merge, fishid, timestamp)

# 1. Najdi všechny relevantní události
events <- full_data_merge[
  navigation_event_id %in% c(100,37,92,347,73,60,47),
  .(event_time = timestamp, nav_id = navigation_event_id),
  by = fishid
]

# 2. Pro každou událost vyber záznamy do 10 minut po ní a připoj navigation_event_id
post_events <- events[
  , full_data_merge[
    fishid == .BY$fishid &
      timestamp > event_time &
      timestamp <= event_time + 120
  ][
    , nav_source := nav_id
  ],
  by = .(fishid, event_time, nav_id)
]

# 3. Volitelně: filtruj jen ty s transition_phase == "after_last" nebo NA
final_selection<- post_events[transition_phase == "after_last" | is.na(transition_phase)]
final_selection[, depth_interp := na.approx(depth, x = timestamp, na.rm = FALSE), by = .(fishid, navigation_event_id)]

setDT(final_selection)
setnames(final_selection, make.unique(names(final_selection)))  # z fishid udělá "fishid" a "fishid.1"
setkey(final_selection, fishid, nav_source, timestamp)

unique(final_selection$nav_source)
ggplot() +
  geom_raster(data= hloubka_df , aes(x = x, y = y, fill =lyr.1), na.rm = T) +
  geom_path(data = final_selection, aes (x = easting, y = northing, group = as.factor(nav_source)), col = "red")+
 # geom_point(data = final_selection, aes (x = easting, y = northing, group = nav_source), col = "black")+
  coord_sf(crs = "EPSG:32633", xlim = c(461600, 462800), ylim =c(5410000, 5411000)) + # Použijte coord_sf pro správné souřadnice
  
  geom_point(data = home_points, aes (x = home_x, y = home_y), col = "red", size = 2)

# r <- pgGetRast(con, c("teri", "lakerastertiles"), rast = "rast", clauses = "WHERE lr_id = 1")
# crs(r) <- "EPSG:32633"  # nebo jiný kód podle skutečné projekce
# plot(r)
# 
# hladina <- r
# # Řekněme známá kóta hladiny (měřeno třeba při plném stavu):
# hladina <- 469.0  # metry nad mořem, uprav podle reality
# 
# # Raster je nadmořská výška dna:
# dem <- r  # tvoje data
# 
# # Výpočet hloubky
# hloubka <- hladina - dem
# 
# # Mimo nádrž nastavíme NA (negativní hodnoty nejsou hloubka)
# hloubka[hloubka <= 0] <- NA
# 
# # Plot hloubkové mapy
# plot(hloubka, main = "Hloubka vody v nádrži (m)")
# 
# # Uložení do GeoTIFF
# writeRaster(hloubka, "~/Teri/pikeperch_navigation/data/Rimov_depth_raster_33N.tif", overwrite = TRUE)
library(raster)
hloubka <- raster("~/Teri/pikeperch_navigation/data/Rimov_depth_raster_33N.tif")
plot(hloubka)
# Původní raster je r, zmenšíme ho třeba 5x
hloubka_aggr <- terra::aggregate(hloubka, fact = 7, fun = mean, na.rm = TRUE)

# Převod do data.frame pro ggplot:
hloubka_df <- as.data.frame(hloubka_aggr, xy = TRUE, na.rm = TRUE)
#fwrite(hloubka_df, "~/Teri/pikeperch_navigation/data/Rimov_depth_raster_agg.csv")

home_points <- last_seg_sub[last_seg_sub$navigation_event_id == nav_ev_id &
                              last_seg_sub$transition_phase == "after_last", .(
  home_x = mean(x_sleep, na.rm = TRUE),
  home_y = mean(y_sleep, na.rm = TRUE),
  home_depth = 0
)]  # nebo navigation_event_id, pokud chceš podle události


nav_ev_id <- c(87, 72, 103, 127, 123,137,152)
#nav_ev_id <- c(100,37,92,347,73,60,47)
# 7. VIZUALIZACE S GGPLOT2
ggplot() +
  geom_path(data = last_seg_sub[last_seg_sub$navigation_event_id == nav_ev_id &
                                  last_seg_sub$transition_phase == "after_last"  ], aes (x = easting, y = northing, group = navigation_event_id, col = as.factor(shore_50_merged)))+
  geom_point(data = last_seg_sub[last_seg_sub$navigation_event_id == nav_ev_id &
                                  last_seg_sub$transition_phase == "after_last"  ], aes (x = easting, y = northing, group = navigation_event_id, col = as.factor(shore_50_merged)))+
  
  geom_point(data = home_points, aes (x = home_x, y = home_y), col = "red", size = 2)+
  scale_colour_manual(values = c("orange", "black"))+
  scale_fill_viridis_c(
    option = "D",             # Můžete zkusit jiné (A, B, C, D, E)
    name = "Hloubka [m]",
    na.value = "transparent", # Důležité pro odstranění šedého/žlutého pozadí
    direction = -1,           # Pro hloubku často chcete, aby tmavé bylo hluboké
    limits = c(0, 45),        # Omezte rozsah legendy na smysluplné hodnoty (max. hloubka Římova je kolem 43m)
    breaks = c(0, 10, 20, 30, 40) # Pěkné zlomy v legendě
  ) +
  coord_sf(crs = "EPSG:32633", xlim = c(461600, 462800), ylim =c(5410000, 5411000)) + # Použijte coord_sf pro správné souřadnice
  labs(
    #title = "Hloubka vody v nádrži Římov (agregováno) [m]",
    x = "",
    y = "",
    col = "breh < 50m"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray95", linetype = "dotted"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

unique(final_selection$navigation_event_id)
ggplot() +
  geom_raster(data= hloubka_df , aes(x = x, y = y, fill =lyr.1), na.rm = T) +
  geom_path(data = final_selection[navigation_event_id %in% c(37, 73)], aes (x = easting, y = northing, group = navigation_event_id, col = as.factor(navigation_event_id)), linewidth = 1.5)+
  #geom_point(data = final_selection, aes (x = easting, y = northing, group = navigation_event_id))+
  
  geom_point(data = home_points, aes (x = home_x, y = home_y), col = "red1", size = 2)+
  scale_colour_manual(values = c("orangered", "black"))+
  scale_fill_viridis_c(
    option = "D",             # Můžete zkusit jiné (A, B, C, D, E)
    name = "Hloubka [m]",
    na.value = "transparent", # Důležité pro odstranění šedého/žlutého pozadí
    direction = -1,           # Pro hloubku často chcete, aby tmavé bylo hluboké
    limits = c(0, 45),        # Omezte rozsah legendy na smysluplné hodnoty (max. hloubka Římova je kolem 43m)
    breaks = c(0, 10, 20, 30, 40) # Pěkné zlomy v legendě
  ) +
  coord_sf(crs = "EPSG:32633", xlim = c(461600, 462800), ylim =c(5409000, 5411000)) + # Použijte coord_sf pro správné souřadnice
  labs(
    #title = "Hloubka vody v nádrži Římov (agregováno) [m]",
    x = "",
    y = "",
    col = "id"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray95", linetype = "dotted"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggplot() +
  geom_raster(data= hloubka_df , aes(x = x, y = y, fill =lyr.1), na.rm = T) +
  geom_path(data = final_selection[navigation_event_id %in% c(347)], aes (x = easting, y = northing, group = navigation_event_id, col = as.factor(navigation_event_id)), linewidth = 1.5)+
  #geom_point(data = final_selection, aes (x = easting, y = northing, group = navigation_event_id))+
  
  geom_point(data = home_points, aes (x = home_x, y = home_y), col = "red1", size = 2)+
  scale_colour_manual(values = c("orangered", "black"))+
  scale_fill_viridis_c(
    option = "D",             # Můžete zkusit jiné (A, B, C, D, E)
    name = "Hloubka [m]",
    na.value = "transparent", # Důležité pro odstranění šedého/žlutého pozadí
    direction = -1,           # Pro hloubku často chcete, aby tmavé bylo hluboké
    limits = c(0, 45),        # Omezte rozsah legendy na smysluplné hodnoty (max. hloubka Římova je kolem 43m)
    breaks = c(0, 10, 20, 30, 40) # Pěkné zlomy v legendě
  ) +
  coord_sf(crs = "EPSG:32633", xlim = c(461600, 462800), ylim =c(5410100, 5411000)) + # Použijte coord_sf pro správné souřadnice
  labs(
    #title = "Hloubka vody v nádrži Římov (agregováno) [m]",
    x = "",
    y = "",
    col = "id"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray95", linetype = "dotted"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# plotly 

library(plotly)
library(terra)
# Výběr trajektorie
traj_data <- last_seg_sub %>%
  filter(navigation_event_id == nav_ev_id, transition_phase == "after_last")

traj_data <- final_selection
# Začátek trajektorie - y souřadnice
y_start <- traj_data$northing[1]

# Převod na terra objekt (pokud už není)
hloubka_terra <- rast(hloubka_aggr)

# Získání hranic pomocí terra syntaxe
e <- ext(hloubka_terra)
xmin <- e$xmin
xmax <- e$xmax
ymax <- e$ymax
ymin <- y_start - 50

# Oříznutí (vše v rámci balíčku terra)
hloubka_crop <- crop(hloubka_terra, ext(xmin, xmax, ymin, ymax))

# Nepovinné: Vyhlazení rastru pro hladší okraje
#hloubka_crop <- disagg(hloubka_crop, fact = 3, method = "bilinear")

# Převod na matrix a korekce osy Y
z_matrix <- as.matrix(hloubka_crop, wide = TRUE)
z_matrix <- z_matrix[nrow(z_matrix):1, ]

x_vals <- terra::xFromCol(hloubka_crop, 1:ncol(hloubka_crop))
y_vals <- terra::yFromRow(hloubka_crop, 1:nrow(hloubka_crop))
y_vals <- rev(y_vals)  # nutné převrátit osy, aby seděla matice a prostor

# Převod nadmořské výšky dna
hladina <- 0  # nastav podle hladiny v m n.m.
z_matrix_invert <- z_matrix - hladina

# Pokud máš absolutní nadmořskou výšku dna v trajektorii, dopočítáš hloubku ryby:
# traj_data$depth_bottom_elev <- ...  # Pokud máš v datech
# traj_data$depth_utm <- traj_data$depth_bottom_elev - traj_data$depth

# Pokud ryba plave podle hloubky vůči hladině:
traj_data$depth_utm <- traj_data$depth_interp  # Negativní, aby hloubka šla dolů

camera_fixed <- list(
  eye = list(x = 2, y = 5, z = -1),  # směr odkud se koukáš
  center = list(x = 0, y = 0, z = 0)  # střed pohledu
)

# seřadit a vyčistit (pro jistotu)
# seřazení a odstranění NA
traj_fix <- traj_data[order(traj_data$timestamp), ]
traj_fix <- subset(traj_fix, !is.na(easting) & !is.na(northing) )

#traj_fix <- traj_fix[order(traj_fix$timestamp), ]
#traj_fix <- subset(traj_fix, !is.na(easting) & !is.na(northing) & !is.na(depth_utm))

# Výpočet mean pozice pro každý fishid / trajektorii
home_points <- traj_fix[, .(
  home_x = mean(x_sleep, na.rm = TRUE),
  home_y = mean(y_sleep, na.rm = TRUE),
  home_depth = 0
)]  # nebo navigation_event_id, pokud chceš podle události


# barvy (ideálně hex)
traj_colors <- c(
  "100" = "#006400",  # chartreuse4
  "37"  = "#7FFF00",  # chartreuse
  "92"  = "#CD3333",  # firebrick3
  "347" = "#FF3030",  # firebrick1
  "73"  = "#00BFFF",  # deepskyblue
  "60"  = "#00B2FF",  # deepskyblue2 approx
  "47"  = "#0064A6"   # deepskyblue4 approx
)

fallback_col <- "#000000"  # černá pro případy, kdy barva chybí

p <- plot_ly() %>%
  add_surface(
    x = x_vals, y = y_vals, z = -z_matrix_invert,
    colorscale = "Cividis",
    cmin = min(-z_matrix, na.rm = TRUE),
    cmax = max(-z_matrix, na.rm = TRUE),
    showscale = F,
    colorbar = list(title = list(text = "Bottom depth [m]")),
    opacity = 1,
    showlegend = FALSE
  )

# přidej trajektorie – deterministicky a s jistou barvou
ids <- sort(unique(traj_fix$nav_source))

for (id in ids) {
  col <- traj_colors[as.character(id)]
  if (is.na(col)) col <- fallback_col
  col <- unname(col)             # <<< klíčové: zahoď jméno vektoru
  
  df <- traj_fix[nav_source == id][order(timestamp)]
  
  p <- p %>%
    add_trace(
      data = df,
      x = ~easting, y = ~northing, z = ~(-depth_interp),
      type = "scatter3d", mode = "lines",
      line = list(width = 7, color = col),
      name = paste("nav", id),
      legendgroup = paste0("nav_", id),
      showlegend = F
    )
}

# home body
p <- p %>%
  add_trace(
    data = home_points,
    x = ~home_x, y = ~home_y, z = 0,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3, color = "black", symbol = "x"),
    name = "home", legendgroup = "home", showlegend = TRUE
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "Easting (m)"),
      yaxis = list(title = "Northing (m)"),
      zaxis = list(title = "Depth (m)"),
      bgcolor = "white",
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.5),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 0.6),
                    center = list(x = 0, y = 0, z = -0.1))
    ),
    showlegend = F
  )

p

library(htmlwidgets)
library(htmlwidgets)
# p je tvůj plotly objekt
p <- p %>% 
  htmlwidgets::prependContent(
    htmltools::tags$h3("Interactive 3D bathymetric model with pikeperch return trajectories"),
    htmltools::tags$p("Users can rotate, zoom and explore the bathymetry and recorded return paths to the resting platform (RP)."),
    htmltools::tags$ul(
      htmltools::tags$li(htmltools::tags$span(style="color:green;", "Green – final navigation phase along the shoreline")),
      htmltools::tags$li(htmltools::tags$span(style="color:blue;", "Blue – final offshore phase followed by nearshore navigation to the RP")),
      htmltools::tags$li(htmltools::tags$span(style="color:red;", "Red – precise final approach from offshore to the RP"))
    )
  )

saveWidget(p,
           "~/Teri/pikeperch_navigation/trajektorie_graf.html",
           selfcontained = TRUE)


