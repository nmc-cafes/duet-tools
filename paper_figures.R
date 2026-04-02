library(here)
library(tidyverse)
library(terra)
library(tidyterra)
library(scico)
library(patchwork)

arr_to_rst <- function(dir, name, crop_ext = NULL, n_row=908, n_col=1208){
  arr <- read.table(paste0(dir,"/",name,".txt"))
  mat <- as.matrix(arr)
  rst <- rast(nrow = n_row,
              ncol = n_col,
              ext = ext(c(0,n_col,0,n_row)),
              vals = mat,
              names = name)
  rst <- flip(rst, direction="vertical")
  if(!is.null(crop_ext)){
    rst <- crop(rst, crop_ext)
  }
  return(rst)
}

rst_to_df <- function(rst){
  df <- as.data.frame(rst, xy=T)
  name <- names(rst)[1]
  before_after <- strsplit(name, "_")[[1]][1]
  fuel_type <- strsplit(name, "_")[[1]][2]
  df$before_after <- before_after
  df$fuel_type <- fuel_type
  names(df)[3] <- "loading"
  return(df)
}

crop_ext <- ext(c(0,100,25,125)) #crop to make heterogeneity easier to see

original_loading <- arr_to_rst(here("figures-data","Arrays"),"original_loading",crop_ext)
original_grass <- arr_to_rst(here("figures-data","Arrays"),"original_grass",crop_ext)
original_litter <- arr_to_rst(here("figures-data","Arrays"),"original_litter",crop_ext)
original_deciduous <- arr_to_rst(here("figures-data","Arrays"),"original_deciduous",crop_ext)
original_coniferous <- arr_to_rst(here("figures-data","Arrays"),"original_coniferous",crop_ext)

calibrated_loading <- arr_to_rst(here("figures-data","Arrays"),"calibrated_loading",crop_ext)
calibrated_grass <- arr_to_rst(here("figures-data","Arrays"),"calibrated_grass",crop_ext)
calibrated_litter <- arr_to_rst(here("figures-data","Arrays"),"calibrated_litter",crop_ext)
calibrated_deciduous <- arr_to_rst(here("figures-data","Arrays"),"calibrated_deciduous",crop_ext)
calibrated_coniferous <- arr_to_rst(here("figures-data","Arrays"),"calibrated_coniferous",crop_ext)

original_stack <- c(original_loading, original_grass, original_litter, original_deciduous, original_coniferous)
calibrated_stack <-  c(calibrated_loading, calibrated_grass, calibrated_litter, calibrated_deciduous, calibrated_coniferous)

loading_stack <- c(original_loading, calibrated_loading)
grass_stack <- c(original_grass, calibrated_grass)
litter_stack <- c(original_litter, calibrated_litter)
deciduous_stack <- c(original_deciduous, calibrated_deciduous)
coniferous_stack <- c(original_coniferous, calibrated_coniferous)

###########
og_loading_df <- rst_to_df(original_loading)
og_grass_df <- rst_to_df(original_grass)
og_litter_df <- rst_to_df(original_litter)
og_deciduous_df <- rst_to_df(original_deciduous)
og_coniferous_df <- rst_to_df(original_coniferous)

cal_loading_df <- rst_to_df(calibrated_loading)
cal_grass_df <- rst_to_df(calibrated_grass)
cal_litter_df <- rst_to_df(calibrated_litter)
cal_deciduous_df <- rst_to_df(calibrated_deciduous)
cal_coniferous_df <- rst_to_df(calibrated_coniferous)

duet_df <- bind_rows(og_loading_df,
                     og_grass_df,
                     og_litter_df,
                     og_deciduous_df,
                     og_coniferous_df,
                     cal_loading_df,
                     cal_grass_df,
                     cal_litter_df,
                     cal_deciduous_df,
                     cal_coniferous_df) %>%
  mutate(before_after = factor(before_after,
                               levels = c("original","calibrated"),
                               labels = c("Original", "Calibrated")),
         fuel_type = factor(fuel_type,
                            levels = c("loading","grass","coniferous","deciduous"),
                            labels = c("Total","Grass/\nHerbaceous","Coniferous\nLitter","Deciduous\nLitter"))) %>%
  mutate(y = y-25)

#### Plot
total <- duet_df %>%
  filter(fuel_type == "Total") %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  facet_grid(before_after~fuel_type) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,5),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  scale_y_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA,
                                        color = NA),
        legend.position = 'bottom',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

grass <- duet_df %>%
  filter(fuel_type == "Grass/\nHerbaceous") %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  facet_grid(before_after~fuel_type) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   # limits = c(0,4),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  scale_y_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA,
                                        color = NA),
        legend.position = 'bottom',
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

coniferous <- duet_df %>%
  filter(fuel_type == "Coniferous\nLitter") %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  facet_grid(before_after~fuel_type) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,5),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  scale_y_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA,
                                        color = NA),
        legend.position = 'bottom',
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.title = element_blank())

deciduous <- duet_df %>%
  filter(fuel_type == "Deciduous\nLitter") %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  facet_grid(before_after~fuel_type) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   # limits = c(0,4),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  scale_y_continuous(expand=c(0,0), breaks = c(25,50,75,100)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA,
                                        color = NA),legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

figure <- coniferous + deciduous + grass + total + plot_layout(nrow = 1) +
  plot_annotation(caption = bquote('Fine Fuel Loading (kg m'^-2*")"),
                  theme = theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))))
figure
ggsave("figure1.jpg", figure, path = here("figures-data","Plots"), width=9, height=9*0.6)

###############
## Landfire

# Import calibrated arrays
landfire_calibrated <- arr_to_rst(here("figures-data","Arrays"), "calibrated_loading_landfire", crop_ext)
landfire_litter <- arr_to_rst(here("figures-data","Arrays"), "calibrated_litter_landfire", crop_ext)
landfire_grass <- arr_to_rst(here("figures-data","Arrays"),"calibrated_grass_landfire", crop_ext)

cal_loading_lf_df <- rst_to_df(landfire_calibrated) %>% mutate(y = y-25)
cal_litter_lf_df <- rst_to_df(landfire_litter) %>% mutate(y = y-25)
cal_grass_lf_df <- rst_to_df(landfire_grass) %>% mutate(y = y-25)

# Import landfire arrays that calibration targets were derived from
landfire_loading <- arr_to_rst(here("figures-data","Arrays"), "landfire_loading", crop_ext=NULL, n_row=45, n_col=47)
landfire_fueltype <- arr_to_rst(here("figures-data","Arrays"), "landfire_fueltype", crop_ext=NULL, n_row=45, n_col=47)

landfire_loading_df <- rst_to_df(landfire_loading)
landfire_fueltype_df <- rst_to_df(landfire_fueltype) %>%
  select(-fuel_type) %>%
  rename(fuel_type = loading) %>%
  mutate(fuel_type = factor(fuel_type,
                            levels = c(1,-1,0),
                            labels = c("Grass","Litter","Both")))

# Plot
original_duet <- og_loading_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,5),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.title = element_blank())

calibrated_loading <- cal_loading_lf_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,5),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

calibrated_litter <- cal_litter_lf_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,2),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

calibrated_grass <- cal_grass_lf_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,2),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.title = element_blank())

landfire_targets <- landfire_loading_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  scale_fill_scico(palette="managua",
                   direction = -1,
                   limits = c(0,2),
                   na.value = scico(palette = "managua",
                                    direction = -1,
                                    n=2)[2]) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

landfire_fueltypes <- landfire_fueltype_df %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=fuel_type)) +
  scale_fill_manual(values = c("forestgreen","orange3","gray")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(fill = "Fuel Type") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# Patchwork them together

bottom_left <- original_duet
bottom_right <- calibrated_loading
bottom_panel <- bottom_left + bottom_right + plot_layout(guides = "collect")

top_left <- (landfire_fueltypes + landfire_targets)
top_right <- (calibrated_grass + calibrated_litter) + plot_layout(guides = "collect") & theme(legend.position = "right")
top_panel <- (top_left | top_right)

final_plot <- top_panel / bottom_panel + plot_layout(heights = c(1,3)) +
  plot_annotation(caption = bquote('Fine Fuel Loading (kg m'^-2*")"),
                  theme = theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))) &
  theme(plot.tag.position = c(0,1))

# Save to file and fix plot spacing externally
ggsave("figure2_raw.jpg", final_plot, path = here("figures-data","Plots"), width=9, height=9)
