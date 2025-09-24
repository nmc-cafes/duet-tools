library(here)
library(tidyverse)
library(terra)
library(tidyterra)
library(scico)
library(patchwork)

arr_to_rst <- function(dir, name, crop_ext = NULL){
  arr <- read.table(paste0(dir,"/",name,".txt"))
  mat <- as.matrix(arr)
  rst <- rast(nrow = 908,
              ncol = 1208,
              ext = ext(c(0,1208,0,908)),
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

# ggplot() +
#   geom_spatraster(data = loading_stack) +
#   facet_wrap(~lyr) +
#   scale_fill_scico(palette = "managua", 
#                    direction = -1,
#                    limits = c(0,4),
#                    na.value = "yellow") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_fixed() +
#   theme_bw()
# 
# ggplot() +
#   geom_spatraster(data = grass_stack) +
#   facet_wrap(~lyr) +
#   scale_fill_scico(palette = "managua", 
#                    direction = -1,
#                    limits = c(0,2),
#                    na.value = "yellow") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_fixed() +
#   theme_bw()
# 
# ggplot() +
#   geom_spatraster(data = deciduous_stack) +
#   facet_wrap(~lyr) +
#   scale_fill_scico(palette = "managua", 
#                    direction = -1,
#                    limits = c(0,0.4),
#                    na.value = "yellow") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_fixed() +
#   theme_bw()
# 
# ggplot() +
#   geom_spatraster(data = coniferous_stack) +
#   facet_wrap(~lyr) +
#   scale_fill_scico(palette = "managua", 
#                    direction = -1,
#                    limits = c(0,2),
#                    na.value = "yellow") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_fixed() +
#   theme_bw()

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
                            labels = c("Total\n","Grass/\nHerbaceous","Coniferous\nLitter","Deciduous\nLitter"))) %>%
  mutate(y = y-25)

#### Plot
total <- duet_df %>%
  filter(fuel_type == "Total\n") %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=loading)) +
  facet_grid(before_after~fuel_type) +
  scale_fill_scico(palette="managua",
                   direction = 1,
                   # limits = c(0,4),
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
  theme(legend.position = 'bottom',
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
                   direction = 1,
                   # limits = c(0,4),
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
  theme(legend.position = 'bottom',
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
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'bottom',
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
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed() +
  labs(y = "Y (m)",
       x = "X (m)",
       fill = bquote('Fine Fuel\nLoading (kg m'^-2*")")) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

figure <- coniferous + deciduous + grass + total + plot_layout(nrow = 1) + 
  plot_annotation(title = bquote('Fine Fuel\nLoading (kg m'^-2*")"),
                  tag_levels = 'A')
figure
ggsave("paper_figure.jpg", figure, path = here("figures-data","Plots"), width=9, height=9*0.6)
