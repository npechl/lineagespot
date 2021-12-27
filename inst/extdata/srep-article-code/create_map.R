

# load libraries -----------------------------

# citation("maps")
# citation("mapdata")
# citation("dplyr")
# citation("ggplot2")
# citation("ggrepel")

library(maps)
library(mapdata)

library(dplyr)

library(ggplot2)
library(ggrepel)

# Get the world polygon and extract Greece ------------------------

gr <- map_data("world") %>% filter(region == "Greece")

data <- world.cities %>% 
    filter(country.etc == "Greece") %>% 
    filter(name == "Thessaloniki")


ggplot() +
    
  geom_polygon(data = gr, aes(x = long, y = lat, group = group), 
               fill = "grey", alpha = 0.8) +
    
  geom_text_repel(data = data, aes(x = long, y = lat, label = name), 
                  size = 5) +
    
  geom_point(data = data , aes(x = long, y = lat), 
             color = "red", size=3) +
    
  theme_void() + 
  coord_map() 


