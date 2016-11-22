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

# version 0.9.9

# NB: The viterbi function will only work when R is run with the --max-ppsize=500000 command-line option (no rstudio for you!)

### save random seed for future reproducibility
runif(1)
my.seed  <- .Random.seed
write.table(my.seed,".randomseed")

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

write.ops <- function () {
  out.df <- data.frame()
  for (n in names(op.args)) {
	v <<- op.args[[n]]
	df.line <- data.frame(
	  option=n,
	  value=v
	)
	out.df <- rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.txt",row.names=F)
}

op.args <- list(
  "name" = NA,
  "st.start" = 6,
  "st.end" = 20,
  "rhmm.iter" = 250,
  "rhmm.iter.init" = 10,
  "gff.path" = NA
)

read.ops(input.args)

if (is.na(op.args[["name"]])) {
  cat("Error: name must be set (use --name=[name])\n")
  q()
}

if (is.na(op.args[["gff.path"]])) {
  cat("Error: path to chromatin GFFs must be set (use --gff.path=[name])\n")
  q()
}

write.ops()


### Load libs
library(RHmm)
source('~/Dropbox/R/read.gff.r')

### read input data
# Data is read from a file called "input.variables.txt" in the source directory

name <- op.args[["name"]] # must conform to R variable naming conventions -- no spaces, etc
st.sta <- as.integer(op.args[["st.start"]]) # minimum number of HMM states to fit
st.end <- as.integer(op.args[["st.end"]]) # maximum number of HMM states to fit
rhmm.n.iter <- as.integer(op.args[["rhmm.iter"]])
rhmm.n.iter.init <- as.integer(op.args[["rhmm.iter.init"]])
base.path <- op.args[["gff.path"]]

polii.path <- Sys.glob(paste(base.path,"*olII*",sep="/"))
h1.path <- Sys.glob(paste(base.path,"*H1*",sep="/"))
hp1.path <- Sys.glob(paste(base.path,"*HP1*",sep="/"))
pc.path <- Sys.glob(paste(base.path,"*Pc*",sep="/"))
brm.path <- Sys.glob(paste(base.path,"*rm*",sep="/"))

log.file <- vector()
log.file <- paste()

cat(sprintf("Processing sample: %s\t(fitting from %d to %d states)\n",name,st.sta,st.end))

### Make master dataframes ...
cat("Building dataframes ...\n")
polii.dat <- read.gff(polii.path,paste("polii.",name,sep=""))
brm.dat <- read.gff(brm.path,paste("brm.",name,sep=""))
h1.dat <- read.gff(h1.path,paste("h1.",name,sep=""))
hp1.dat <- read.gff(hp1.path,paste("hp1.",name,sep=""))
pc.dat <- read.gff(pc.path,paste("pc.",name,sep=""))

data.all <- merge(polii.dat,brm.dat,by=c("chr","start","end"))
data.all <- merge(data.all,h1.dat,by=c("chr","start","end"))
data.all <- merge(data.all,hp1.dat,by=c("chr","start","end"))
data.all <- merge(data.all,pc.dat,by=c("chr","start","end"))

data.sorted <- data.all[order(data.all$chr,data.all$start),]
data.na <- data.sorted
data.na[data.na == 0] <- NA
data.na.ex <- na.omit(data.na)

# global list variables
chrom.states <- list()

### fit HMMs to data
cat("Fitting HMMs ...\n")
for (state in st.sta:st.end) {
  cat(paste("  Fitting",state,"states ...\n"))
  chrom.states[[paste(name,".na.",state,sep="")]] <- HMMFit(na.exclude(data.matrix(data.na[,4:8])),nStates=state,control=list(verbose=1,nInit=rhmm.n.iter,nIterInit=rhmm.n.iter.init))
  
  # Can only use viterbi with R --max-ppsize=500000
  chrom.states[[paste(name,".na.vit.",state,sep="")]] <- viterbi(chrom.states[[paste(name,".na.",state,sep="")]],data.matrix(na.exclude(data.na[,4:8])))
  
  save.image()
}

save.image()
cat("All done.\n")

