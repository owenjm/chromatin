# Chromatin state processing script (the single-state bit)
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

# version 0.9.5

# This script should be run after the "chromatin.all.states" R script.
# After analysing the data from the above script, a single state is now chosen
# to analyse further.

# This script should be run after and inherit all variables from the chromatin.all.states script.

library(ape)

# X11
x11()

# Other source files
source("~/Dropbox/R/gene.fdr.r")
source("~/Dropbox/R/genexp.modal.prop.ii.r")
source("~/Dropbox/R/chromatin.functions.iii.r")

### Read CLI options
input.args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

read.ops <- function (x) {
  for (op in x) {
    op <- gsub("^--","",op)
    y <- unlist(strsplit(op,"="))
    
    if (y[1] == "help") {
      cat("Options:\n")
      for (n in names(single.op.args)) {
        cat(paste("  ",n,"=",single.op.args[[n]],"\n",sep=""))
      }
      q()
    }
    
    if (!is.null(single.op.args[[ y[1] ]])) {
	  print(y[2])
	  if ( grepl("^[[:digit:]]*$", y[2])) {
		single.op.args[[ y[1] ]] <<- as.integer(y[2])
	  } else {
		single.op.args[[ y[1] ]] <<- y[2]
	  }
    } else {
      cat("Error: Option",y[1],"not recognised ...\n")
      q()
    }
  }
}

write.ops <- function () {
  out.df <- data.frame()
  for (n in names(single.op.args)) {
	v <<- single.op.args[[n]]
	df.line <- data.frame(
	  option=n,
	  value=v
	)
	out.df <- rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.single.txt",row.names=F)
}

single.op.args <- list(
  "input" = Sys.glob("single*.r"),
  "data" = "all.states.RData",
  "genes.file" = "/mnt/data/Genomes/dmel_release/DmR6/DmR6.genes.gff",
  "tfs" = "/home/owen/dros.transcription.factors.tfs-update.txt",
  "pad" = 250
)

read.ops(input.args)

if (is.na(single.op.args[["input"]])) {
  cat("Error: input file must be set (use --input=[single.state.r file])\n\n")
  q()
}

write.ops()

simp.names <- list()
simp.cols <- list()
simp.states <- list()

data.file <- tryCatch({
    load(single.op.args[["data"]])
  },  warning = function (e) {
    cat("Error: unable to read .RData file (use --data=[.RData file] to set)\n\n")
    q()
  }
)

# input file is read from the commandline option:
source(single.op.args[["input"]])
print(single.op.args)

# Date -- for working directory
adate <- format(Sys.time(), "%Y-%m-%d_%H:%M%P")
wd <- paste(name,adate,sep=".")
dir.create(wd)
setwd(wd)

# read TFs data
all.tfs <- read.table(single.op.args[["tfs"]],quote="")


###########################################
### redo gene expression (may not have been completed if performing partial analysis)
chrom.states[[paste(name,".na.states.",target.state,sep="")]] <- make.states.gff(chrom.states[[paste(name,".na.vit.",target.state,sep="")]]$states,data.na.ex)
all.modal <- gene.exp.modal(chrom.states[[paste(name,".na.states.",target.state,sep="")]], genes.file=single.op.args[["genes.file"]], pad=single.op.args[["pad"]])
chrom.states[[paste(name,".na.genes.",target.state,sep="")]] <- all.modal[["exp"]]
gene.count <- all.modal[["count"]]

print(gene.count)

### save data
write.table(chrom.states[[paste(name,".na.genes.",target.state,sep="")]],paste("single.state.dat"))


###########################################
### Heatmaps
state.mean.hmap <- data.frame(name=c("polII","Brm","H1","HP1","Pc"))
for (state in 1:target.state) {
  state.mean.hmap[,paste("s",state,sep="")] <- chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$distribution$mean[[state]]
}

state.mean.sc <- t(state.mean.hmap[,2:(target.state+1)] )
state.mean.hmap.sc <- t( scale(state.mean.sc, center = FALSE, scale = apply(state.mean.sc,2, sd, na.rm = TRUE) ) )

hmap(as.matrix(state.mean.hmap[c(3,5,4,2,1),unlist(simp.states[[name]])+1]))
dev.copy(png,paste("State.",target.state,".heatmap.means.png",sep=""),height=500,width=80*target.state);dev.off()
dev.copy(pdf,paste("State.",target.state,".heatmap.means.pdf",sep=""),height=5,width=0.8*target.state);dev.off()

hmap(as.matrix(state.mean.hmap.sc[c(3,5,4,2,1),unlist(simp.states[[name]])]))
dev.copy(png,paste("State.",target.state,".heatmap.means.sc.png",sep=""),height=500,width=80*target.state);dev.off()
dev.copy(pdf,paste("State.",target.state,".heatmap.means.sc.pdf",sep=""),height=5,width=0.8*target.state);dev.off()

###########################################
### Pie charts
# Simplified states
gene.pie.simplified(
                    data.frame(
                              name=chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$name,
                              state=chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$state
                              ),
                    cols=simp.cols[[name]],
                    names=simp.names[[name]],
                    simplified.states=simp.states[[name]],
                    state.order=c(1:length(simp.states[[name]])),
                    num.genes=gene.count
                  )
dev.copy(png,paste("simplified.state.pie.png",sep=""),width=600,height=600);dev.off()
  dev.copy(pdf,paste("simplified.state.pie.pdf",sep=""));dev.off()

# All states pie
full.cols <- vector()
for (i in 1:length(simp.states[[name]])) {
  for (j in unlist(simp.states[[name]][i])) {
    full.cols[j] <- simp.cols[[name]][i]
  }
}

gene.pie.simplified(
                    chrom.states[[paste(name,".na.genes.",target.state,sep="")]],
                    names=c(1:target.state),
                    simplified.states=c(1:target.state),
                    state.order=unlist(simp.states[[name]]),
                    cols=full.cols,
                    num.genes=gene.count
                    )
dev.copy(png,paste("all.states.pie.png",sep=""),width=600,height=600);dev.off()
  dev.copy(pdf,paste("all.states.pie.pdf",sep=""));dev.off()



### Gene lists
# simplified states
for ( i in c(1:length(simp.states[[name]])) ) {
  out.names <- chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$name[ chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$state %in% simp.states[[name]][[i]] ]
  write.table(out.names,paste(simp.names[[name]][i],".genes.txt",sep=""),row.names=F,col.names=F,quote=F)
}

# all states
for (i in c(1:target.state)) {
  out.names <- chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$name[ chrom.states[[paste(name,".na.genes.",target.state,sep="")]]$state == i ]
  write.table(out.names,paste("state",i,"genes.txt",sep="."),row.names=F,col.names=F,quote=F)
}


###########################################
### gene expression
exp.genes <- gene.exp(data.na.ex, genes.file=single.op.args[["genes.file"]])
write.table(exp.genes, paste(name,"gene.expression.dat",sep="."))
write.table(exp.genes$name[exp.genes$FDR < 0.01], paste(name,"expressed.genes.txt",sep="."),row.names=F,col.names=F,quote=F)


###########################################
#### simplified states viterbi -- for output and comparison
states.reverse <- vector()
for (i in 1:length(simp.states[[name]])) {
  for (j in unlist(simp.states[[name]][i])) {
    states.reverse[j] <- i
  }
}
simp.vit.gff <- chrom.states[[paste(name,".na.states.",target.state,sep="")]]
simp.vit.gff$state <- states.reverse[simp.vit.gff$state]

write.table(data.frame(chrom.states[[paste(name,".na.states.",target.state,sep="")]]$chr,".",".",chrom.states[[paste(name,".na.states.",target.state,sep="")]]$start,chrom.states[[paste(name,".na.states.",target.state,sep="")]]$end,chrom.states[[paste(name,".na.states.",target.state,sep="")]]$state),paste("all.states.gff",sep=""), col.names=F, row.names=F, quote=F)

write.table(data.frame(simp.vit.gff$chr,".",".",simp.vit.gff$start,simp.vit.gff$end,simp.vit.gff$state),paste("simplified.states.gff",sep=""), col.names=F, row.names=F, quote=F)

simp.vit.genes <- gene.exp.modal(simp.vit.gff, genes.file=single.op.args[["genes.file"]], pad=single.op.args[["pad"]])[["exp"]]

gene.pie.simplified(
                    data.frame(
                              name=simp.vit.genes$name,
                              state=simp.vit.genes$state
                              ),
                    cols=simp.cols[[name]],
                    names=simp.names[[name]],
                    simplified.states=c(1:length(simp.names[[name]])),
                    state.order=c(1:length(simp.states[[name]])),
                    num.genes=gene.count
                  )
dev.copy(png,paste("simp.vit.pie.png",sep=""),width=600,height=600);dev.off()
  dev.copy(pdf,paste("simp.vit.pie.pdf",sep=""));dev.off()

all.states.simp.genes <- chrom.states[[paste(name,".na.genes.",target.state,sep="")]]
all.states.simp.genes$state <- states.reverse[all.states.simp.genes$state]

gene.comp <- data.frame(name = simp.vit.genes$name, state.simp = simp.vit.genes$state, state.all = all.states.simp.genes$state, prop.simp = simp.vit.genes$prop, prop.all = all.states.simp.genes$prop, gene.size = simp.vit.genes$gene.size)
gene.mismatch <- subset(gene.comp, state.simp != state.all)

gene.mis.hist <- vector()
gene.mis.names <- list()
for (i in 1:length(simp.states[[name]])) {
  gene.mis.hist[i] <- length(gene.mismatch$name[gene.mismatch$state.all == i])
  gene.mis.names[[ simp.names[[name]][i] ]] <- gene.mismatch$name[gene.mismatch$state.all == i]
}

plot(gene.mis.hist)
dev.copy(pdf,paste("gene.mis.hist.pdf",sep=""));dev.off()

###########################################
### state tables
# expressed genes in states
state.table.all.states <- output.state.table(states=chrom.states[[paste(name,".na.genes.",target.state,sep="")]],state.names=simp.names[[name]][states.reverse],exp.genes=exp.genes,state.order=unlist(simp.states[[name]]),fdr=0.01)
write.table(state.table.all.states,"state.table.all.states.txt", row.names=F, col.names=T, quote=F)
capture.output( print(state.table.all.states, print.gap=3), file="state.table.all.states.capture.txt")

# expressed genes in simp.vit states
state.table.simp.vit <- output.state.table(states=simp.vit.genes,state.names=simp.names[[name]],exp.genes=exp.genes,fdr=0.01)
write.table(state.table.simp.vit,"state.table.simp.vit.states.txt", row.names=F, col.names=T, quote=F)
capture.output( print(state.table.simp.vit, print.gap=3), file="state.table.simp.vit.states.capture.txt")


# histograms of proportional coverage
for ( i in c( 1:length(simp.names[[name]]) ) ) {
  simp.vit.state <- subset(simp.vit.genes, state == i)
  hist(simp.vit.state$prop, col=simp.cols[[name]][i], main=simp.names[[name]][i], xlim=c(0,1), xlab="Proportion of gene body covered by modal state")
  dev.copy(pdf,paste("simp.vit.state.coverage.",simp.names[[name]][i],".pdf",sep=""));dev.off()
    dev.copy(png,paste("simp.vit.state.coverage.",simp.names[[name]][i],".png",sep=""));dev.off()
}

hist(simp.vit.genes$prop, col="#666666", main="All states", xlim=c(0,1), xlab="Proportion of gene body covered by modal state")
dev.copy(pdf,paste("simp.vit.state.coverage.all.pdf",sep=""));dev.off()
  dev.copy(png,paste("simp.vit.state.coverage.all.png",sep=""));dev.off()


###########################################
#full gene lists for all simplified states (for GO term analysis, etc)
for (i in 1:length(simp.states[[name]])) {
  out.names <- simp.vit.genes$name[ simp.vit.genes$state == i ]
  write.table(out.names,paste("simp.vit.genes.state",simp.names[[name]][i],"genes.txt",sep="."),row.names=F,col.names=F,quote=F)
}

### Boxplots
boxplot(simp.vit.genes$gene.size ~ simp.vit.genes$state,outline=F)
dev.copy(pdf,paste("simp.vit.genes.states.boxplot.pdf",sep=""));dev.off()


### Gene state boxplots
for (i in c(1:target.state)) {
  data.cut <- data.na.ex[ chrom.states[[paste(name,".na.vit.",target.state,sep="")]]$states==i, c(4:8)]
  title <- paste("State ",i)
  boxplot(data.cut, col=c("orange","darkred","grey50","darkgreen","darkblue"),ylim=c(-5,5),main=title,names=c("polII","Brm","H1","HP1","Pc"))
  abline(0,0)
  dev.copy(pdf, paste("detailed.boxplot.state",i,"pdf",sep=".")); dev.off()
}


### reverse lookup lists for colours
all.states <- list()
all.cols <- list()
all.name <- list()
real.cols <- list()
real.names <- list()
for (i in 1:length(simp.states[[name]])) {
  all.states[[name]] <- append(all.states[[name]], simp.states[[name]][[i]])
  for (j in 1: length(simp.states[[name]][[i]])) {
    all.cols[[name]] <- append(all.cols[[name]], simp.cols[[name]][i])
    all.name[[name]] <- append(all.name[[name]], simp.names[[name]][i])
  }
}

for (i in 1:length(all.states[[name]])) {
  real.cols[[name]] <- append(real.cols[[name]], all.cols[[name]][ all.states[[name]] == i ])
  real.names[[name]] <- append(real.names[[name]], all.name[[name]][ all.states[[name]] == i ] )
}



### State transitions plotting, clustering and trees
heatmap(chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$transMat,col=colorRampPalette(c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1)))(256))
dev.copy(pdf, paste("transMat.heatmap",name,target.state,"pdf",sep=".")); dev.off()

hc <- hclust(dist(chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$transMat))

plot( as.phylo(hc), type="radial", tip.color=real.cols[[name]])
dev.copy(pdf, paste("hclust.transMat",name,target.state,"pdf",sep=".")); dev.off()

state.hc <- hclust(dist(t(state.mean.hmap.sc)))
plot( as.phylo(state.hc), type="radial", tip.color=real.cols[[name]])
dev.copy(pdf, paste("hclust.stateMeans.scaled",name,target.state,"pdf",sep=".")); dev.off()



### protein binding boxplots -- unordered
for (i in c(4:8)) {
  prot.name <- names(data.na.ex)[i]
  prot.col=c("orange","darkred","grey50","darkgreen","darkblue")[i-3]
  binding <- list()
  for (j in c(1:target.state)) {
	binding[[j]] <- data.na.ex[ chrom.states[[paste(name,".na.vit.",target.state,sep="")]]$states == j, i]
  }
  title <- prot.name
  boxplot(binding,ylim=c(-5,5),col=prot.col, main=title,names=c(1:target.state))
  abline(0,0)
  dev.copy(pdf, paste("prot.binding.boxplot",prot.name,"pdf",sep=".")); dev.off()
}

### protein binding boxplots -- ordered
for (i in c(4:8)) {
  prot.name <- names(data.na.ex)[i]
  prot.col=c("orange","darkred","grey50","darkgreen","darkblue")[i-3]
  binding <- list()
  k <- 1;
  for (j in all.states[[name]]) {
	binding[[k]] <- data.na.ex[ chrom.states[[paste(name,".na.vit.",target.state,sep="")]]$states == j, i]
	k <- k+1
  }
  title <- prot.name
  boxplot(binding,ylim=c(-5,5),col=prot.col, main=title,names=all.states[[name]])
  abline(0,0)
  dev.copy(pdf, paste("prot.binding.boxplot",prot.name,"ordered.pdf",sep=".")); dev.off()
}


### write scaled states data
state.mean.hmap.name <- data.frame(name=c("polII","Brm","H1","HP1","Pc"))
for (state in 1:target.state) {
  state.mean.hmap.name[,paste(name,"s",state,sep=".")] <- chrom.states[[paste(name,".na.",target.state,sep="")]]$HMM$distribution$mean[[state]]
}

state.mean.hmap.name.sc <- t( scale(t(state.mean.hmap.name[,2:(target.state+1)] ), center=F) )
write.table(state.mean.hmap.name.sc,"state.mean.hmap.sc.named.ordered.dat")
write.table(real.cols[[name]],"state.mean.hmap.sc.cols.dat")


### barplots
# proportion of genes expressed
barplot(state.table.all.states$pc.of.state,col=all.cols[[name]],main=paste("Proportion of",name,"genes expressed by state"),names.arg=state.table.all.states$state)
dev.copy(pdf,paste("barpolot.prop.of",name,"genes.expr.by.state.pdf",sep="."),height=5,width=0.5*target.state);dev.off()

par(las=1)
par(mar=c(5,8,4,2))
barplot(rev(state.table.simp.vit$pc.of.state),col=rev(simp.cols[[name]]),main=paste("Proportion of",name,"genes expressed by state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T)
dev.copy(pdf,paste("barpolot.prop.of",name,"simp.vit.genes.expr.by.state.pdf",sep="."),height=5,width=5);dev.off()

# Enrichment of TFs
par(las=1)
par(mar=c(5,8,4,2))
barplot(rev(-log(state.table.simp.vit$p.value)),col=rev(simp.cols[[name]]),main=paste("Enrichment of TFs in state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T,xlab="-log(P value) Fischer's exact test"); abline(v=-log(0.05/length(simp.cols[[name]])));
dev.copy(pdf,paste("barpolot.enrichment.of.tfs.log-p-value",name,"simp.vit.pdf",sep="."),height=5,width=0.5*target.state);dev.off()

# Percentage of TFs
par(las=1)
par(mar=c(5,8,4,2))
barplot(rev(state.table.simp.vit$pc.tfs.in.state),col=rev(simp.cols[[name]]),main=paste("Proportion of TFs in state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T,xlab="Percentage of TFs")
dev.copy(pdf,paste("barpolot.percentage.of.tfs",name,"simp.vit.pdf",sep="."),height=5,width=0.5*target.state);dev.off()


### save data
write.table(simp.vit.genes,paste("single.state.simp.dat"))
write.table(simp.vit.genes$name,"all.analysed.genes.txt",row.names=F,col.names=F,quote=F)

### coloured IGV tracks
convert.hex.rgb <- function (x) {
  hex.r <- substr(x,2,3)
  hex.g <- substr(x,4,5)
  hex.b <- substr(x,6,7)
  dec.r <- strtoi(hex.r,16L)
  dec.g <- strtoi(hex.g,16L)
  dec.b <- strtoi(hex.b,16L)
  out <- paste(dec.r,",",dec.g,",",dec.b,sep="")
  return(out)
}


condense.gff <- function (input.df) {
	merged <- data.frame()
	
	chr.last <- ""
	state.start <- NA
	state.end <- NA
	state.current <- NA
	for (i in c(1:nrow(input.df))) {
	  chr <- input.df[i,1]
	  start <- input.df[i,2]
	  end <- input.df[i,3]
	  state <- input.df[i,4]
	  
	  if (chr == chr.last) {
		if (state == state.current) {
		  state.end <- end
		} else {
		  # save state block
		  merged <- rbind(merged, data.frame(chr=chr, start=state.start, end=state.end, state=state.current))
		  
		  state.start <- start
		  state.end <- end
		  state.current <- state
		}
	  } else {
		# new chromosome
		state.start <- start
		state.end <- end
		state.current <- state
	  }
	  
	  chr.last <- chr
	}
	
	return(merged)
}

make.states.bed <- function (datf, simp) {
  # expect "states gff" input: chr start end state
  datf[,"zero"] <- 0
  datf[,"plus"] <- "+"
  datf$cols <- sapply(simp.vit.gff$state, function(x) convert.hex.rgb(simp.cols[[name]][x]) )
  
  states.bed <- data.frame(datf$chr,datf$start,datf$end,datf$state,datf$zero,datf$plus,datf$start,datf$end,datf$cols)
  return(states.bed)
}

write('track name="" description="" visibility=2 itemRgb="On"', paste(name,"igv.cols.bed",sep="."))
write.table(make.states.bed(simp.vit.gff, simp.cols), paste(name,"igv.cols.bed",sep="."), append=T, quote=F, row.names=F,col.names=F)



save.image("single.Rdata")

#cat("Running goterm analysis ...")
#system("perl ~/Dropbox/perl/flymine.goterms.pl *.genes.txt")

setwd("../")
cat("Copying single.state.r file to working directory ...")
file.copy(single.op.args[["input"]],wd)

cat("All done.\n")



