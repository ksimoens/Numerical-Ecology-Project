g1 <- ggmap::ggmap(map,extent="panel", fill="grey")
g1 <- g1 + geom_point(aes(x=coord$lon,y=coord$lat,fill=MEM4), data=geo_dbMEM, shape=21,size=3) + scale_fill_gradient2(low="white", mid="grey", high="black") +
  theme(axis.text.x=element_text(colour="black"))+
   theme(axis.text.y=element_text(colour="black"))+
 labs(x="longitude") +
 labs(x="latitude") +
 labs(x="longitude") +
 labs(y="latitude") +
 theme(panel.border = element_rect(colour="black", fill=NA, size=1),
         axis.title=element_text(size=12,colour="black",family="Helvetica"),
         legend.title=element_text(size=10,colour="black",family="Helvetica"))
