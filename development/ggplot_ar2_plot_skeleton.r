library(tidyverse)
library(cyclar)

ggplot() + 
  scale_x_continuous(breaks = seq(-2, 2, 0.25), lim = c(-2, 2)) + 
  scale_y_continuous(breaks = seq(-1, 1, 0.25), lim = c(-1, 1)) + 
  # triangle
  geom_polygon(data = data.frame(x = c(-2, 0, 2, -2), y = c(-1, 1, -1, -1)),
            aes(x, y), fill = NA, color = "black") + 
  # parabola
  geom_function(fun = function(x) -0.25 * x^2) + 
  # vertical line at 0
  geom_segment(data = data.frame(x1 = 0, x2 = 0, y1 = -1, y2 = 1),
               aes(x = x1, xend = x2, y = y1, yend = y2))
  
