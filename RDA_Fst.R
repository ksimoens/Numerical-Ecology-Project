library(vegan)
library(adespatial)

Fst.atl <- read.csv("Output/PCo/PCo.sel.atl.csv",header=T,row.names=1)

env <- read.table("Output/EnvMatrix.txt",header=T,row.names=1)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

Fst.atl <- Fst.atl[match(rownames(AEM),rownames(Fst.atl)),]
env.atl <- env.atl[match(rownames(AEM),rownames(env.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

env.atl.norm <- scale(env.atl, center=T, scale=T)

rda.env <- rda(Fst.atl, env.atl.norm)
R2adj.env <- RsquareAdj(rda.env)
print(anova(rda.env,permutations = how(nperm=9999)))
sel.env <- forward.sel(Fst.atl, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
print(sel.env)

env.atl.norm.sel <- env.atl.norm[,sel.env$order]

rda.linear <- rda(Fst.atl,linear)
R2adj.linear <- RsquareAdj(rda.linear)
print(anova(rda.linear,permutations = how(nperm=9999)))

rda.dbMEM <- rda(Fst.atl, dbMEM)
R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
print(anova(rda.dbMEM,permutations = how(nperm=9999)))
sel.dbMEM <- forward.sel(Fst.atl, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
print(sel.dbMEM)

dbMEM.sel <- dbMEM[,sel.dbMEM$order]

rda.AEM <- rda(Fst.atl, AEM)
R2adj.AEM <- RsquareAdj(rda.AEM)
print(anova(rda.AEM,permutations = how(nperm=9999)))
sel.AEM <- forward.sel(Fst.atl, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
print(sel.AEM)

fractions <- vector()
R2adj <- vector()
degrees <- vector()
Fvalues <- vector()
pvalues <- vector()

rda.part <- varpart(Fst.atl, env.atl.norm.sel, linear, dbMEM.sel)
rda.part.df <- as.data.frame(rbind(rda.part$part$fract,rda.part$part$indfract,rda.part$part$contr1))
rda.part.df <- rda.part.df[rda.part.df$Testable==T,]
rownames(rda.part.df) <- c("[env]","[lin]","[MEM]","[env] + [lin]","[env] + [MEM]","[lin] + [MEM]",
	"[env] + [lin] + [MEM]","[env] | [lin] + [MEM]","[lin] | [env] + [MEM]","[MEM] | [env] + [lin]",
	"[env] | [MEM]","[env] | [lin]","[lin] | [MEM]","[lin] | [env]","[MEM] | [env]","[MEM] | [lin]")

key <- list("[env]"=env.atl.norm.sel, "[lin]"=linear, "[MEM]"=dbMEM.sel)

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
		rda.out <- rda(Fst.atl, variable, conditional)
		print(names(conditional))
	} else {
		rda.out <- rda(Fst.atl, variable)
	}
	return(rda.out)
}

for(i in 1:nrow(rda.part.df)){
	fractions <- c(fractions, rownames(rda.part.df)[i])
	R2adj <- c(R2adj, rda.part.df$Adj.R.square[i])

	rda <- makeRDA(rownames(rda.part.df)[i])
	print(summary(rda))
	an <- anova(rda, permutations = how(nperm=9999))
	deg <- paste0(an$Df[1],"(",an$Df[1]+an$Df[2],")")
	degrees <- c(degrees, deg)
	Fvalues <- c(Fvalues, an$F[1])
	pvalues <- c(pvalues, an$"Pr(>F)"[1])
}

df.parts <- data.frame(fraction=fractions, RsquareAdj=R2adj, dof=degrees, F=Fvalues, "p-value"=pvalues)
print(df.parts)