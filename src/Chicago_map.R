library("ggmap")

# get_stamenmap() does NOT require a google API. 
chicago_bbox = c(left = -87.9, bottom = 41.6, right = -87.5, top = 42.06)
chicago <- get_stamenmap(bbox = chicago_bbox, maptype = "toner-background", color = "bw") 
ggmap(chicago)
save(chicago, file = "./data/chicago_map.RData")
save(chicago_bbox, file = "./data/chicago_bbox.RData")
