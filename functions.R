library(hierfstat)

PCoA <- function(Fst){
	A <- -0.5 * Fst * Fst

	identity <- diag(ncol(A))
	ones <- rep(1,ncol(A))
	side <- identity - ones %*% t(ones) / ncol(A)

	Delta <- side %*% A %*% side

	Evalues <- eigen(Delta)$values
	print(Evalues)

	c <- abs(Evalues[length(Evalues)])
	print(c)
	c.mat <- matrix(rep(c,ncol(A)*ncol(A)),ncol=ncol(A))
	c.mat <- c.mat - c*diag(ncol(A))

	Fst.corr <- sqrt(Fst * Fst + 2*c.mat)
	A <- -0.5 * Fst.corr * Fst.corr
	Delta <- side %*% A %*% side

	Evalues <- eigen(Delta)$values
	Evectors <- eigen(Delta)$vectors
	names <- paste0("PCo",1:ncol(Evectors))
	colnames(Evectors) <- names
	rownames(Evectors) <- rownames(Fst)

  Evalues <- Evalues[which(Evalues > 1e-10)]
  Evectors <- Evectors[,which(Evalues > 1e-10)]

  for(i in 1:ncol(Evectors)){
    Evectors[,i] <- Evectors[,i]*Evalues[i]
  }

	return(list(values=Evalues,vectors=Evectors))
}

countryList <- function(populations){
  reference <- data.frame(
                  pop = c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul",
                          "Jer",
                          "Idr",
                          "Hel",
                          "Ale","Sky","The","Tor",
                          "Cor","Hoo","Kil","Mul","Ven",
                          "Laz","Tar","Sar",
                          "Oos",
                          "Ber","Flo","Sin","Tro",
                          "Vig",
                          "Gul","Kav","Lys"),
                  country = c(rep("GRB",14),"CHA","FRA","DEU",rep("HEL",4),rep("IRL",5),rep("ITA",3),"NDL",rep("NOR",4),"ESP",rep("SVE",3))
                )
  for(i in 1:nrow(reference)){
    populations[populations == reference$pop[i]] <- reference$country[i]
  }

  return(populations)
}

regionList <- function(populations){
  reference <- data.frame(
                  pop = c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul",
                          "Jer",
                          "Idr",
                          "Hel",
                          "Ale","Sky","The","Tor",
                          "Cor","Hoo","Kil","Mul","Ven",
                          "Laz","Tar","Sar",
                          "Oos",
                          "Ber","Tro",
                          "Flo","Sin",
                          "Vig",
                          "Gul","Kav","Lys"),
                  region = c(rep("ATL",17),rep("MED",4),rep("ATL",5),rep("MED",3),rep("ATL",3),rep("SKA",2),"ATL",rep("SKA",3))
                )
  for(i in 1:nrow(reference)){
    populations[populations == reference$pop[i]] <- reference$region[i]
  }

  return(populations)
}

makePlotDF <- function(Evectors,names){
  df.plot <- Evectors[,1:2] %>% as.data.frame()
  rownames(df.plot) <- names
  df.plot$country <- rownames(df.plot)
  df.plot$country <- countryList(df.plot$country)
  df.plot$region <- rownames(df.plot)
  df.plot$region <- regionList(df.plot$region)

  print(df.plot)
  return(df.plot)
}

getColours <- function(pop){
	col_key <- data.frame(
			population = c("CHA","DEU","ESP","FRA","GRB","HEL","IRL","ITA","NDL","NOR","SVE","ATL","MED","SKA"),
			colour = c(rep(c("#C7E020FF","#24868EFF","#440154FF"),3),c("#C7E020FF","#440154FF"),c("#C7E020FF","#24868EFF","#440154FF")),
			shape = c(rep(15,3),rep(17,3),rep(19,3),18,18,15,17,19))

	col_res <- col_key[col_key$population %in% unique(pop),]$colour
	sha_res <- col_key[col_key$population %in% unique(pop),]$shape
	return(list(colour = col_res,shape = sha_res))
}

calcFst <- function(gendata){
	diffstats <- pairwise.WCfst(gendata)

	print(diffstats)

	rownames(diffstats) <- substr(rownames(diffstats),1,3)
	colnames(diffstats) <- substr(colnames(diffstats),1,3)
	diffstats[diffstats < 0] <- 0
	diffstats[is.na(diffstats)] <- 0

	return(diffstats)
}

makeRDA <- function(rowname){
	str_res <- strsplit(rowname, " +")[[1]]
	print(str_res)
	str_res <- str_res[str_res!="+"]
	n <- length(str_res)
	index <- which(str_res=="|")
	if(length(index) == 0){
		index <- n+1
	} else {
		index <- index[1]
	}
	variable <- as.data.frame(key[names(key)==str_res[1]])
	if(index > 2){
		for(i in 2:(index-1)){
			var <- str_res[i]
			variable <- cbind(variable,key[names(key)==var])
		}
	}
	print(names(variable))
	if(index != n+1){
		conditional <- as.data.frame(key[names(key)==str_res[index+1]])
		if(index + 1 < n){
			for(i in (index+2):(n)){
				var <- str_res[i]
				conditional <- cbind(conditional,key[names(key)==var])
			}
		}
		rda.out <- rda(freq.atl, variable, conditional)
		print(names(conditional))
	} else {
		rda.out <- rda(freq.atl, variable)
	}
	return(rda.out)
}