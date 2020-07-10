#This function makes state plots for k-means clustering maps.
#input should be a single dataframe with three columns, long, lat, layer
kmeans_map =
  function(map_df){
    colour = c("10" = "white", "3"="springgreen3", "9"="slateblue2", "1"="tomato2", "5"="goldenrod1", "6"="gray50", "8"="gray17", "4"="brown", "2"="darkslategray4", "7"="deeppink3")
    mapplot =
      ggplot(map_df, aes(x, y))+
        geom_raster(aes(fill = layer))+
        #scale_fill_manual("", values=redo_color, labels = region_names)+
        geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "black", fill = NA)+
        theme_map()+
      scale_fill_manual("", values=colour)
    return(mapplot)
    }
