library(ggplot2)

data <- read.csv("Input/economy.csv",header=T)

data.lan <- data[data$type=="landings",]
data.lan.UK <- data.lan[data.lan$country=="UK",]
data.lan.EU <- data.lan[data.lan$country=="EU",]

data.ear <- data[data$type=="earnings",]
data.ear$value <- data.ear$value/1000
data.ear.UK <- data.ear[data.ear$country=="UK",]
data.ear.EU <- data.ear[data.ear$country=="EU",]

p <- ggplot() + 
		geom_polygon(aes(x=c(data.lan.EU$year,rev(data.lan.EU$year)),y=c(data.lan.EU$value,rep(0,length(data.lan.EU$value))),fill='European Union'),col='black',alpha=0,size=1) +
		geom_polygon(aes(x=c(data.lan.UK$year,rev(data.lan.UK$year)),y=c(data.lan.UK$value,rep(0,length(data.lan.UK$value))),fill='United Kingdom'),col='black',alpha=0.5,size=1) +
		scale_fill_manual("", breaks=c("European Union","United Kingdom"), values=c("white","black")) +
		ylab("landings (tons)") +
		scale_x_continuous(breaks=2007:2019, labels=2007:2019) +
		theme_bw() + theme(panel.grid.minor=element_blank(), axis.title.x=element_blank(), legend.position='top',legend.text=element_text(size=15))

ggsave("economy_landings.png",p,device='png',width=20,height=10,units='cm')

p <- ggplot() + 
		geom_polygon(aes(x=c(data.ear.EU$year,rev(data.ear.EU$year)),y=c(data.ear.EU$value,rep(0,length(data.ear.EU$value))),fill='European Union'),col='black',alpha=0,size=1) +
		geom_polygon(aes(x=c(data.ear.UK$year,rev(data.ear.UK$year)),y=c(data.ear.UK$value,rep(0,length(data.ear.UK$value))),fill='United Kingdom'),col='black',alpha=0.5,size=1) +
		scale_fill_manual("", breaks=c("European Union","United Kingdom"), values=c("white","black")) +
		ylab("landing profit (million EUR)") +
		scale_x_continuous(breaks=2007:2019, labels=2007:2019) +
		theme_bw() + theme(panel.grid.minor=element_blank(), axis.title.x=element_blank(), legend.position='top',legend.text=element_text(size=15))

ggsave("economy_earnings.png",p,device='png',width=20,height=10,units='cm')
		
