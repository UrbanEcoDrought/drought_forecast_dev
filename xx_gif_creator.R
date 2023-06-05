library(magick)
library(ggplot2)
library(dplyr)
library(tidyr)


dir_out <- "figures/test_forecast"

imgs <- list.files(dir_out, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
# img_animated

## save to disk
image_write(image = img_animated,
            path = "figures/test_forecast/ranom_walk.gif")
