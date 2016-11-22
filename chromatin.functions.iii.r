
#### simplified pie
gene.pie.simplified <- function (genes.df,cols=NULL,names=NULL,simplified.states=NULL,state.order=c(1:length(simplified.states)),num.genes=17718) {
    state.len <- vector()
    props <- vector()
    j <- 0
    out.cols<-vector()
    n.states <- length(simplified.states)
    cat("n.states = ",n.states,"\n")
    for (i in state.order) {
        cat("state ",i,":",simplified.states[[i]],"\n")
        j <- j+1
        cat(i,"\n")
        out.cols[j] <- cols[i]
        stlen <- length(genes.df$name[genes.df$state %in% simplified.states[[i]] ])
        state.len[j]=stlen
        fnum=format(round((stlen/num.genes*100),2),nsmall=2)
        props[j]=paste(names[i],"\n",fnum,"%")
    }
    rem <- num.genes-length(genes.df$name)
    state.len[length(state.len)+1] <- rem
    props[length(props)+1] <- paste("N/A\n",format(round((rem/num.genes*100),2),nsmall=2),"%")
    out.cols[length(out.cols)+1] <- "white"
    
    cat(out.cols,"\n")
    
    pie(state.len,col=out.cols,labels=props,border="white")
    return(state.len)
}

### chromatin factors heatmap
hmap <- function (x, st.labels=unlist(simp.states[[name]])) {
  image(c(1:length(x[1,])),c(1:length(x[,1])),t(x),col=colorRampPalette(c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1)))(256),xlab="State",ylab="protein",zlim=c(-max(x),max(x)), axes=F)
  axis(2, 1:length(x[,1]), labels = rev(c("polII","Brm","HP1","Pc","H1")), las = 2, tick = 0)
  axis(1, 1:length(x[1,]), labels = st.labels, las = 1, tick = 0)
}


### output data tables
output.state.table <- function (states,state.names,exp.genes, state.order=c(1:length(state.names)), fdr=0.01) {    
    all.tfs <- read.table("/home/owen/dros.transcription.factors.tfs-update.txt",quote="")
    out.df <- data.frame()
    
    for (i in state.order) {
        gene.state.names <- subset(states, state %in% i, name)[[1]]
        genes.expressed <- gene.state.names[gene.state.names %in% exp.genes$name[ exp.genes$FDR < fdr]]
        
        all.tfs.in.state <- gene.state.names[gene.state.names %in% all.tfs[[1]] ]
        tfs.present.in.state <- length(gene.state.names[gene.state.names %in% all.tfs[[1]] ])
        tfs.expressed.in.state <- length(genes.expressed[ genes.expressed %in% all.tfs[[1]] ])
        
        num.genes.in.state <- length(gene.state.names)
        
        #write.table(genes.expressed[ genes.expressed %in% all.tfs[[1]] ], paste("tfs.expressed.in.state",i,"txt",sep="."),col.names=F,row.names=F,quote=F)
        #write.table(all.tfs.in.state, paste("all.tfs.present.in.state",i,"txt",sep="."),col.names=F,row.names=F,quote=F)
        
        expression.value <- mean(exp.genes[exp.genes$name %in% genes.expressed,2])
        expression.state <- mean(exp.genes[exp.genes$name %in% gene.state.names,2])
      
        pc.exp.st <- length(genes.expressed)/length(gene.state.names) * 100
        pc.exp.tot <- length(genes.expressed)/length(exp.genes[,1]) * 100
        pc.tot.tfs <- length(all.tfs.in.state)/length(all.tfs[[1]]) * 100
        pc.tfs.in.state <- length(all.tfs.in.state)/length(gene.state.names) * 100
        
        num.tfs <- length(all.tfs[,1])
        num.genes <- length(exp.genes[,1])
        
        fisher.p <- fisher.test(matrix(c(
                    tfs.present.in.state,
                    num.tfs - tfs.present.in.state,
                    num.genes.in.state - tfs.present.in.state,
                    num.genes - num.tfs - num.genes.in.state + tfs.present.in.state
                ),2,2),alternative="greater")["p.value"]
        
        df.line <- data.frame(
            state = i,
            state.name = state.names[i],
            num.genes = length(gene.state.names),
            num.exp = length(genes.expressed),
            pc.of.state = pc.exp.st,
            pc.of.total = pc.exp.tot,
            mean.exp.exp = expression.value,
            mean.exp.all = expression.state,
            tfs.in.state = tfs.present.in.state,
            pc.tfs.in.state = pc.tfs.in.state,
            pc.all.tfs = pc.tot.tfs,
            fisher.exact.test = fisher.p,
            tfs.expressed = tfs.expressed.in.state
        )
        
        out.df <- rbind(out.df, df.line)
    }
    
    #names(out.df) <- c(
    #    "State",
    #    "State name",
    #    "Number of genes",
    #    "Number expressed",
    #    "% of state",
    #    "% of total",
    #    "Mean exp (expressed)",
    #    "Mean exp (all genes in state)",
    #    "TFs in state",
    #    "% TFs in state",
    #    "% of all TFs",
    #    "Fisher exact P-value",
    #    "TFs expressed"
    #)
    # print(out.df)
    return(out.df)
}


### make states
make.states.gff <- function (states,df) {
    # expect df to have first three columns: chr, start, end
    gff <- data.frame(chr=df$chr, start=df$start, end=df$end, state=states)
    return(gff)
}
