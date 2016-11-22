Mode <- function(x) {
  ux <- unique(x)
  tx <- tabulate(match(x, ux))
  mx <- ux[which.max(tx)]
  cx <- tx[which.max(tx)]
  return(c(mx,cx,length(x)))
}

gene.exp.modal <- function (input.df, genes.file="/mnt/data/Genomes/dmel_release/DmR6/DmR6.genes.gff", pad=0, subset.genes=T) {
  
  print(genes.file)
  
  cat("Reading genes data file ...\n")
  genes <- read.table(genes.file, comment.char="#", sep="\t", quote="", stringsAsFactors=F)
  
  names(genes) <-  c('chr','source','type','start','end','score','strand','c','details')
  
  print(subset.genes)
  
  # only subset if there is a type termed "gene"
  if ( (subset.genes) & (any(genes$type == 'gene'))) {
    genes <- subset(genes, type=='gene')
  }
  
  total <- length(genes$chr)
  
  genes$name <- sapply(genes$details, FUN = function (x) {regmatches(x,gregexpr("(?<=Name=).*?(?=;)", x, perl=T))} )
  genes <- genes[,c('chr','start','end','strand','name')]
  
  avg.exp <- data.frame()
  avg <- vector(length=1)
  avg.exp <- avg.exp[0,]
  
  count <- 0
    
  # unroll chromosomes for speed:
  for (chromo in unique(genes$chr)) {
    input.chr <- subset(input.df, chr==chromo)
        
    if (is.null(input.chr)) { next }
    
    genes.chr <- subset(genes, chr==chromo)
    
    for (i in c(1:length(genes.chr$name))) {
      expr <- data.frame(input.chr[ (input.chr$start <= genes.chr[i,"end"]) 
                                      & (input.chr$end >= genes.chr[i,"start"]), ] )
      if (length(expr[,1]) == 0) {next}
      
      # trim to gene boundaries ...
      expr$start[1] <- genes.chr[i,"start"]-pad
      expr$end[length(expr[,1])] <- genes.chr[i,"end"]+pad
      
      # gene length (not required here)
      gene.len <- genes.chr[i,"end"]-genes.chr[i,"start"]
      
      # roll through each row
      avg <- vector()
      for (j in c(1:nrow(expr))) {
        avg <- append(avg,rep(expr[j,4],(expr$end[j]-expr$start[j])))
      }
      
      if (is.null(avg)) { cat("Error!",chromo,genes.chr[i,"start"],genes.chr[i,"end"],avg,"\n");next }
      
      m <- Mode(avg)
      modeavg <- m[1]
      mprop <- m[2]/m[3]
      
      new.df <- data.frame(name=as.character(genes.chr[i,"name"]), state=modeavg, prop=mprop, gene.size=gene.len)
      
      avg.exp <- rbind(avg.exp,new.df);
      count <- count+1
      
      if (count %% 50 == 0) {cat(paste(count,"/",total,"genes averaged ...\r"))}
    }
  }
  
  cat("\nAll done.\n\n");
  return(list("exp" = avg.exp, "count" = total));
}
