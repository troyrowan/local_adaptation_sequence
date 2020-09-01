#This function makes state plots for k-means clustering maps.
#input should be a single dataframe with three columns, long, lat, layer
kmeans_map =
  function(map_df, colors = NULL, region = layer){
    mapplot =
      ggplot(map_df, aes(x, y))+
        geom_raster(aes(fill = region))+
        #scale_fill_manual("", values=redo_color, labels = region_names)+
        geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "black", fill = NA)+
        theme_map()+
      scale_fill_manual("", values=colors)
    return(mapplot)
    }


mapplotting2 <- function(data, title = "", legend_title = "", vcolor = NULL, max = 10, min = 0.1, alpha = 0.6, facet = NULL){
  colour = c("10" = "white", "3"="springgreen3", "9"="slateblue2", "1"="tomato2", "5"="goldenrod1", "6"="gray50", "8"="gray17", "4"="brown", "2"="darkslategray4", "7"="deeppink3")
  region_names = c("1" = "Fescue Belt", "2" = "Southeast", "3" = "Forested Mountains", "4" = "Desert", "5" = "Arid Prairie", "6" = "Corn Belt", "7" = "Upper Midwest", "8" = "Rainforest", "9" = "High Plains", "10" = "HOT PLACE")
  all_states<-map_data("state")
  mapplot = ggplot()+
    geom_count(aes(x = data$long, y = data$lat, colour=factor(data$zone)), alpha = alpha)+
    scale_size(range = c(min, max), guide = FALSE)+ #This changes max dot size for altering figure size
    scale_color_manual("Climate Zone", values=colour, labels = region_names[data$zone])+
    #scale_color_manual("Climate Zone", values=colour, labels = paste("Zone ", levels(data$zone), " : ", table(data$zone), sep = ""))+
    ggtitle(title)+
    geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "black", fill = NA)+
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme_map()
  return(mapplot)
}

mapplotting <- function(data, title = "", legend_title = "", vcolor = NULL, max = 10, min = 0.1, alpha = 0.6, facet = NULL){
  colorkey = c("HighPlains"="springgreen3", "UpperMidwest"="slateblue2", "Desert"="tomato2", "AridPrairie"="goldenrod1", "FescueBelt"="gray17", "Rainforest" = "brown", "CornBelt"="brown", "Southeast"="darkslategray4", "ForestedMountains"="deeppink3", "Foothills" = "gray50", "ExtremeDesert" = "gray50")
  all_states<-map_data("state")
  mapplot = ggplot()+
    geom_count(aes(x = data$long, y = data$lat, color=factor(data$region)), alpha = alpha)+
    scale_size(range = c(min, max), guide = FALSE)+ #This changes max dot size for altering figure size
    scale_color_manual("Climate Zone", values=colorkey[data$region])+
    #scale_color_manual("Climate Zone", values=colour, labels = paste("Zone ", levels(data$zone), " : ", table(data$zone), sep = ""))+
    ggtitle(title)+
    geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "black", fill = NA)+
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme_map()
  return(mapplot)
}
