# Chromatin state processing script (comparisons between chromatin states)
# Copyright © 2014, Owen Marshall

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

# version 0.9

# This script should be run after the "chromatin.all.states.r" and "chromatin.single.state.r"
# R scripts.  Paths to the two analysis directories to be compared need to be supplied.
# An optional list of genes can also be supplied to subset the data

### Read CLI options
input.args <- commandArgs(trailingOnly = TRUE)

read.ops <- function (x) {
  for (op in x) {
    op <- gsub("^--","",op)
    y <- unlist(strsplit(op,"="))
    
    if (y[1] == "help") {
      cat("Options:\n")
      for (n in names(op.args)) {
        cat(paste("  ",n,"=",op.args[[n]],"\n",sep=""))
      }
      q()
    }
    
    if (!is.null(op.args[[ y[1] ]])) {
      op.args[[ y[1] ]] <<- y[2]
    } else {
      cat("Error: Option",y[1],"not recognised ...\n")
      q()
    }
  }
}


op.args <- list(
  "n1" = NA,
  "n2" = NA,
  "n3" = NA,
  "black.name" = "black",
  "genes" = NA,
  "all.states" = F,
  "arrow.width" = 15,
  "tfs" = "/home/owen/dros.transcription.factors.tfs-update.txt"
)

read.ops(input.args)

n1.path <- op.args[["n1"]]
n2.path <- op.args[["n2"]]
n3.path <- op.args[["n3"]]
genes.file <- op.args[["genes"]]


### Load libs
#x11()

library(Gmisc)
library(grid)


### Read genes list
gene.list <- vector()
if (!is.na(genes.file)) {
  gene.list <- read.table(genes.file, quote="")
  cat("Gene list found!\n")
  print(gene.list)
}


### Read tfs
all.tfs <- read.table(op.args[["tfs"]],quote="")
print(head(all.tfs[[1]]))

### Read input data
simp.names <- list()
simp.cols <- list()
simp.states <- list()
simp.order <- list()
expr <- list()
expr.all <- list()

cat(paste(n1.path,"\n",n2.path,"\n"))

n1.files <- Sys.glob(paste(n1.path,"single*.r",sep="/"))
source(n1.files[1])
n1 <- name

n2.files <- Sys.glob(paste(n2.path,"single*.r",sep="/"))
source(n2.files[1])
n2 <- name

n3.files <- Sys.glob(paste(n3.path,"single*.r",sep="/"))
source(n3.files[1])
n3 <- name


# read in gene and state data files
expr[[n1]] <- read.table(paste(n1.path,"single.state.simp.dat",sep="/"))
expr[[n2]] <- read.table(paste(n2.path,"single.state.simp.dat",sep="/"))
expr[[n3]] <- read.table(paste(n3.path,"single.state.simp.dat",sep="/"))

# exp[[n1]] and exp[[n2]] are read in via the source script
names(expr[[n1]]) <- c("name",paste(n1,"state",sep="."))
names(expr[[n2]]) <- c("name",paste(n2,"state",sep="."))
names(expr[[n3]]) <- c("name",paste(n3,"state",sep="."))

print(length(expr.all[[n1]]$name))
print(length(expr.all[[n2]]$name))
print(length(expr.all[[n3]]$name))

st.df.2 <- merge(expr[[n1]], expr[[n2]], by=c("name"))
st.df.3 <- merge(expr[[n2]], expr[[n3]], by=c("name"))

transMtrx <- function (n1,n2,st.df) {
  ### reverse lookup lists for colours
  # I think I tied myself in knots here, but it seems to work ...
  all.states <- list()
  all.cols <- list()
  all.name <- list()
  real.cols <- list()
  real.names <- list()
  for (i in 1:length(simp.states[[n1]])) {
    all.states[[n1]] <- append(all.states[[n1]], simp.states[[n1]][[i]])
    for (j in 1: length(simp.states[[n1]][[i]])) {
      all.cols[[n1]] <- append(all.cols[[n1]], simp.cols[[n1]][i])
      all.name[[n1]] <- append(all.name[[n1]], simp.names[[n1]][i])
    }
  }
  
  for (i in 1:length(all.states[[n1]])) {
    real.cols[[n1]] <- append(real.cols[[n1]], all.cols[[n1]][ all.states[[n1]] == i ])
    real.names[[n1]] <- append(real.names[[n1]], all.name[[n1]][ all.states[[n1]] == i ] )
  }
  
  for (i in 1:length(simp.states[[n2]])) {
    all.states[[n2]] <- append(all.states[[n2]], simp.states[[n2]][[i]])
    for (j in 1: length(simp.states[[n2]][[i]])) {
      all.cols[[n2]] <- append(all.cols[[n2]], simp.cols[[n2]][i])
      all.name[[n2]] <- append(all.name[[n2]], simp.names[[n2]][i])
    }
  }
  
  for (i in 1:length(all.states[[n2]])) {
    real.cols[[n2]] <- append(real.cols[[n2]], all.cols[[n2]][ all.states[[n2]] == i ])
    real.names[[n2]] <- append(real.names[[n2]], all.name[[n2]][ all.states[[n2]] == i ] )
  }
  
  ### Transitions with vit.simple
  i1 <- n1
  i2 <- n2
  
  # keep to a list of provided genes if required
  if (length(gene.list)>0) {
    st.df <- subset(st.df, name %in% gene.list[[1]])
  }
  
  genes.total <- nrow(st.df)
  write.table(genes.total,file="genes.total.txt",row.names=F,col.names=F,quote=F)
  
  n1.test <- with(st.df,get(paste(i1,"state",sep=".")))
  n2.test <- with(st.df,get(paste(i2,"state",sep=".")))
  
  n1.tab <- factor(n1.test, labels=simp.names[[i1]][sort(unique(n1.test))])
  n2.tab <- factor(n2.test, labels=simp.names[[i2]][sort(unique(n2.test))])
  
  transition_mtrx <- table(n1.tab, n2.tab)
  write.table(transition_mtrx, paste(n1,n2,"test.transition.mtrx.dat",sep="."))
  
  return(transition_mtrx)
}



plotTrans <- function (n1,n2,st.df,max.wd) {
  ### reverse lookup lists for colours
  # I think I tied myself in knots here, but it seems to work ...
  all.states <- list()
  all.cols <- list()
  all.name <- list()
  real.cols <- list()
  real.names <- list()
  for (i in 1:length(simp.states[[n1]])) {
    all.states[[n1]] <- append(all.states[[n1]], simp.states[[n1]][[i]])
    for (j in 1: length(simp.states[[n1]][[i]])) {
      all.cols[[n1]] <- append(all.cols[[n1]], simp.cols[[n1]][i])
      all.name[[n1]] <- append(all.name[[n1]], simp.names[[n1]][i])
    }
  }
  
  for (i in 1:length(all.states[[n1]])) {
    real.cols[[n1]] <- append(real.cols[[n1]], all.cols[[n1]][ all.states[[n1]] == i ])
    real.names[[n1]] <- append(real.names[[n1]], all.name[[n1]][ all.states[[n1]] == i ] )
  }
  
  for (i in 1:length(simp.states[[n2]])) {
    all.states[[n2]] <- append(all.states[[n2]], simp.states[[n2]][[i]])
    for (j in 1: length(simp.states[[n2]][[i]])) {
      all.cols[[n2]] <- append(all.cols[[n2]], simp.cols[[n2]][i])
      all.name[[n2]] <- append(all.name[[n2]], simp.names[[n2]][i])
    }
  }
  
  for (i in 1:length(all.states[[n2]])) {
    real.cols[[n2]] <- append(real.cols[[n2]], all.cols[[n2]][ all.states[[n2]] == i ])
    real.names[[n2]] <- append(real.names[[n2]], all.name[[n2]][ all.states[[n2]] == i ] )
  }
  
  ### Transitions with vit.simple
  i1 <- n1
  i2 <- n2
  
  # keep to a list of provided genes if required
  if (length(gene.list)>0) {
    st.df <- subset(st.df, name %in% gene.list[[1]])
  }
  
  n1.test <- with(st.df,get(paste(i1,"state",sep=".")))
  n2.test <- with(st.df,get(paste(i2,"state",sep=".")))
  
  n1.tab <- factor(n1.test, labels=simp.names[[i1]][sort(unique(n1.test))])
  n2.tab <- factor(n2.test, labels=simp.names[[i2]][sort(unique(n2.test))])
  
  n1.num.st <- length(levels(n1.tab))
  n2.num.st <- length(levels(n2.tab))
  max.num.states <- max(n1.num.st, n2.num.st)
  
  arrow.trans.cols <- simp.cols[[n1]][as.integer(levels(factor(n1.test)))]
  print(arrow.trans.cols)
  box.names <- levels(n1.tab)
  if (n1.num.st < max.num.states) {
    for (i in c((n1.num.st+1):max.num.states)) {
      levels(n1.tab) <- c(levels(n1.tab),paste("null",i,sep="."))
    }
    
    if (length(arrow.trans.cols) < max.num.states) {
      for (i in c((length(arrow.trans.cols)+1):max.num.states)) {
        print(i)
        arrow.trans.cols <- c(arrow.trans.cols,"#000000")
        print(arrow.trans.cols)
      }
    }
    
    box.names <- levels(n2.tab)
  }
  
  if (n2.num.st < max.num.states) {
    for (i in c((n2.num.st+1):max.num.states)) {
      levels(n2.tab) <- c(levels(n2.tab),paste("null",i,sep="."))
    }
    #arrow.trans.cols <<- simp.cols[[n1]]
    box.names <- levels(n1.tab)
  }
  
  transition_mtrx <- table(n1.tab, n2.tab)
  
  arrow.trans <- vector()
  for (n1st in arrow.trans.cols) {
    arrow.trans <- c(arrow.trans, c(rep(n1st, times=length(levels(n2.tab)))))
  }
  
  # fix box text
  
  box.names <- gsub("TrxG permissive","TrxG",box.names, perl=T, ignore.case=T)
  box.names <- gsub("Null permissive","non-TrxG",box.names, perl=T, ignore.case=T)
  box.names <- gsub("HP1.*","HP1",box.names, perl=T, ignore.case=T)
  box.names <- gsub("PcG repressive","PcG R",box.names, perl=T, ignore.case=T)
  box.names <- gsub("TrxG repressive","TrxG R",box.names, perl=T, ignore.case=T)
  box.names <- gsub("Null repressive",op.args[["black.name"]],box.names, perl=T, ignore.case=T)
  box.names <- gsub("Basal",op.args[["black.name"]],box.names, perl=T, ignore.case=T)
  
  # Simp.vit Transition plot
  pdf(paste("transitions.from",i1,"to",i2,"width",op.args[["arrow.width"]],"pdf",sep="."))
  transitionPlot(
                 transition_mtrx,
                 type_of_arrow = "simple",
                 #box_width = 1/3,
                 min_lwd = unit(1, "mm"),
                 max_lwd = unit(max.wd, "mm"),
                 fill_start_box=simp.cols[[n1]][as.integer(levels(factor(n1.test)))],
                 box_txt=box.names,
                 fill_end_box=simp.cols[[n2]][as.integer(levels(factor(n2.test)))],
                 tot_spacing=0,
                 arrow_clr=arrow.trans,
                 new_page=F
                 )
  dev.off()
  
  
  
  ### Transition tables
  curr.path=paste("./transitions.from",i1,"to",i2,sep=".")
  dir.create(curr.path)
  trans.no <- data.frame()
  
  for (state in 1:length(simp.names[[i1]])) {
    
    total.in.state <- length(st.df$name[ with(st.df,get(paste(i2,"state",sep="."))) == state ])
    
    for (i in 1:length(simp.names[[i2]])) {
      # write tables:
      temp <- st.df[ with(st.df,get(paste(i1,"state",sep="."))) == state & with(st.df,get(paste(i2,"state",sep="."))) == i ,]
      
      if (nrow(temp) == 0) {next}
      
      trans.no <- rbind(trans.no,data.frame(transition=paste(simp.names[[i1]][state]," --> ", simp.names[[i2]][i]), genes=nrow(temp)))
      
      temp.names <- temp$name
      write.table(temp.names,file=paste(curr.path,paste("genes",i1,"state",simp.names[[i1]][state],"to",i2,"state",simp.names[[i2]][i],"txt",sep="."),sep="/"),row.names=F,col.names=F,quote=F)
      
      tf.names <- subset(temp, name %in% all.tfs[[1]])
      
      if (length(tf.names$name) > 0) {
        write.table(tf.names$name,file=paste(curr.path,paste("tfs",i1,"state",simp.names[[i1]][state],"to",i2,"state",simp.names[[i2]][i],"txt",sep="."),sep="/"),row.names=F,col.names=F,quote=F)
      }
    }
  }
  print(head(trans.no))
  
  write.table(trans.no[order(-trans.no$genes),],file=paste(i1,"to",i2,"summary.transitions.table.csv",sep="."), row.names=F, quote=F, sep="\t")
}

mat1 <- transMtrx(n1,n2,st.df.2)
mat2 <- transMtrx(n2,n3,st.df.3)

print(mat1)
print(mat2)

max.m1 <- max(mat1)
max.m2 <- max(mat2)

max.wd.m1 <- op.args[["arrow.width"]]
max.wd.m2 <- op.args[["arrow.width"]]

if (max.m1 > max.m2) {
  max.wd.m2 <- (max.m2/max.m1) * as.numeric(op.args[["arrow.width"]]) 
} else if (max.m2 > max.m1) {
  max.wd.m1 <- (max.m1/max.m2) * as.numeric(op.args[["arrow.width"]])
}

plotTrans(n1,n2,st.df.2,max.wd.m1)
plotTrans(n2,n3,st.df.3,max.wd.m2)


