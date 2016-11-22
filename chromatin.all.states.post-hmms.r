# Chromatin state processing script (the all-states bit)
# Copyright © 2014-15, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

# version 0.9.6

# NB: The viterbi function will only work when R is run with the --max-ppsize=500000 command-line option (no rstudio for you!)

#x11()
library(ape)
load(".RData");


source("~/Dropbox/R/chromatin.functions.iii.r")
source("~/Dropbox/R/genexp.modal.prop.ii.r")
source("~/Dropbox/R/gene.fdr.r")

### Assign modal states to genes
st.end.max <- 0
for (state in st.sta:st.end) {
  if ( is.null(chrom.states[[paste(name,".na.vit.",state,sep="")]])) {
    next;
  } else {
    st.end.max <- state
  }
  
  chrom.states[[paste(name,".na.states.",state,sep="")]] <- make.states.gff(chrom.states[[paste(name,".na.vit.",state,sep="")]]$states,data.na.ex)
}

### Plot state means for all states
curr.path <- "state.mean.plots"
dir.create(curr.path)

### Plot state means for all states using matrix mode

for (target.state in c(st.sta:st.end.max)) {
  state.mean.hmap <- data.frame(name=c("polII","Brm","H1","HP1","Pc"))
  for (state in 1:target.state) {
    state.mean.hmap[,paste("s",state)] <- chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$distribution$mean[[state]]
  }
  
  state.mean.sc <- t(state.mean.hmap[,2:(target.state+1)] )
  state.mean.hmap.sc <- t( scale(state.mean.sc, center = FALSE, scale = apply(state.mean.sc,2, sd, na.rm = TRUE) ) )
  
  ### State transitions plotting, clustering and trees
  pdf(paste("transMat.heatmap",name,target.state,"pdf",sep="."))
  heatmap(chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$transMat,col=colorRampPalette(c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1)))(256))
  dev.off()
  
  hc <- hclust(dist(chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$transMat))
  
  pdf(paste("hclust.transMat",name,target.state,"pdf",sep="."))
  #plot( as.phylo(hc), type="radial")
  plot( as.phylo(hc), type="unrooted", lab4ut="axial")
  dev.off()
  
  pdf(paste("hclust.stateMeans.scaled",name,target.state,"pdf",sep="."))
  state.hc <- hclust(dist(t(state.mean.hmap.sc)))
  #plot( as.phylo(state.hc), type="radial")
  plot( as.phylo(state.hc), type="unrooted", lab4ut="axial")
  dev.off()

  pdf(paste(curr.path,"/","State.",target.state,".heatmap.means.sc.pdf",sep=""),height=5,width=0.8*target.state)
  hmap(as.matrix(state.mean.hmap.sc[c(3,5,4,2,1),]), st.labels=c(1:target.state))
  #dev.copy(png,paste(curr.path,"/","State.",target.state,".heatmap.means.sc.png",sep=""),height=500,width=80*target.state);dev.off()
  dev.off()
  
  state.mean.hmap.order <- data.frame(name=c("polII","Brm","H1","HP1","Pc"))
  for (state in hc$order) {
    state.mean.hmap.order[,paste("s",state)] <- chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$distribution$mean[[state]]
  }
  
  pdf(paste(curr.path,"/","State.",target.state,".heatmap.means.order-transMat.sc.pdf",sep=""),height=5,width=0.8*target.state)
  state.mean.order.sc <- t(state.mean.hmap.order[,2:(target.state+1)] )
  state.mean.order.hmap.sc <- t( scale(state.mean.order.sc, center = FALSE, scale = apply(state.mean.order.sc,2, sd, na.rm = TRUE) ) )
  hmap(as.matrix(state.mean.order.hmap.sc[c(3,5,4,2,1),]), st.labels=hc$order)
  #dev.copy(png,paste(curr.path,"/","State.",target.state,".heatmap.means.order-transMat.sc.png",sep=""),height=500,width=80*target.state);dev.off()
  dev.off()
  
  
  state.mean.hmap.order <- data.frame(name=c("polII","Brm","H1","HP1","Pc"))
  for (state in state.hc$order) {
    state.mean.hmap.order[,paste("s",state)] <- chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$distribution$mean[[state]]
  }
  
  pdf(paste(curr.path,"/","State.",target.state,".heatmap.means.order-hclust-state.sc.pdf",sep=""),height=5,width=0.8*target.state)
  state.mean.order.sc <- t(state.mean.hmap.order[,2:(target.state+1)] )
  state.mean.order.hmap.sc <- t( scale(state.mean.order.sc, center = FALSE, scale = apply(state.mean.order.sc,2, sd, na.rm = TRUE) ) )
  hmap(as.matrix(state.mean.order.hmap.sc[c(3,5,4,2,1),]), st.labels=state.hc$order)
  #dev.copy(png,paste(curr.path,"/","State.",target.state,".heatmap.means.order-hclust-state.sc.png",sep=""),height=500,width=80*target.state);dev.off()
  dev.off()
}


single.state.template <- paste("
## Enter the target state and sample name
target.state <-",st.end.max,"
name <-\"",name,"\"

### assign names and colours to states, and simplify
# (if not simplifying, then just identify all states)

# state names
simp.names[[name]] <- c(
	\"TrxG permissive\",
	\"Null permissive\",
	\"HP1 repressive\", 
	\"PcG mixed\", 
	\"PcG repressive\", 
	\"TrxG repressive\", 
	\"Null repressive\"
)

# state colours (generally leave as is)
simp.cols[[name]] <- c(
	\"#dc3912\", # TrxG p
	\"#FF9900\", # Null p
	\"#4e9a06\", # Hp1 r
	\"#3366cc\", # PcG mix
	\"#2050aa\", # PcG r
	\"#996666\", # TrxG r
	\"#666666\"  # Null r
)

# list of simplified states, with each element being a vector of the states in the simplified state
simp.states[[name]] <- list(  
	c(),	# Trxg permissive
	c(),	# Null permissive
	c(),    # HP1 repressive
	c(),    # PcG mixed
	c(),	# PcG repressive
	c(),    # Trxg repressive
	c() 	# Null repressive
)
",sep="")

write.table(single.state.template,paste("single.state",st.end.max,"r",sep="."),quote=F,row.names=F,col.names=F)

### redo gene expression (may not have been completed if performing partial analysis)
chrom.states[[paste(name,".na.states.",target.state,sep="")]] <- make.states.gff(chrom.states[[paste(name,".na.vit.",target.state,sep="")]]$states,data.na.ex)
all.modal <- gene.exp.modal(chrom.states[[paste(name,".na.states.",target.state,sep="")]], genes.file="/mnt/data/Genomes/dmel_release/DmR6/DmR6.genes.gff")
chrom.states[[paste(name,".na.genes.",target.state,sep="")]] <- all.modal[["exp"]]
gene.count <- all.modal[["count"]]
exp.genes <- gene.exp(data.na.ex, genes.file="/mnt/data/Genomes/dmel_release/DmR6/DmR6.genes.gff")

# expressed genes in states
state.table.all.states <- output.state.table(states=chrom.states[[paste(name,".na.genes.",target.state,sep="")]],state.names=c(1:target.state),exp.genes=exp.genes,state.order=c(1:target.state),fdr=0.01)
write.table(state.table.all.states,"state.table.all.states.txt", row.names=F, col.names=T, quote=F)
capture.output( print(state.table.all.states, print.gap=3), file="state.table.all.states.capture.txt")


### save data
write.table(chrom.states[[paste(name,".na.genes.",target.state,sep="")]],paste(curr.path,"/","single.state.dat",sep=""))

# all states
for (i in c(1:target.state)) {
  out.names <- chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$name[ chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$state == i ]
  write.table(out.names,paste(curr.path,"/","state.",i,".genes.txt",sep=""),row.names=F,col.names=F,quote=F)
}

save.image("all.states.RData")
cat("All done.\n")