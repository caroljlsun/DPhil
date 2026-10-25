
#Name
#Description
#Date
#User


# Libraries ---------------------------------------------------------------

library(magrittr)
library(Gviz)
library(GenomicRanges)
library(data.table)
library(biomaRt)
library(Repitools)
library(Rsamtools)
library(GenomicAlignments)
library(cowplot)
library(data.table)

# Options -----------------------------------------------------------------
options(ucscChromosomeNames=FALSE)
axis <- GenomeAxisTrack(col="black", fontcolor="black")

# Import data -------------------------------------------------------------

# run multiple times till it works, funky connection
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")



# Function: get genes -----------------------------------------------------


getgenes <- function(chr,from,to,offset,AC009505.2=FALSE) {
    genes1 <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                     mart=e75.genemart,
                                     collapseTranscripts="meta",
                                     filters=list(with_ox_refseq_mrna=T),
                                     chromosome=chr,from=from,to=to)
    genes2 <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                     mart=e75.genemart,
                                     collapseTranscripts="meta",
                                     filters=list(with_ox_refseq_ncrna=TRUE),
                                     chromosome=chr,from=from,to=to)

    genes.df <- rbind(as.data.frame(genes1@range),as.data.frame(genes2@range))
    exons.dt <- data.table::as.data.table(genes.df)
    exons.dt <- exons.dt[ start > from & end < to,]
    exons.dt[,start:=start-offset]
    exons.dt[,end:=end-offset]
    exons.dt[,y:=ifelse(strand=="+",0.2,-0.2)]
    genes.dt <- exons.dt[,.(start=min(start),end=max(end),strand=unique(strand),y=unique(y)),
                         by=c("gene","symbol")]
    genes.dt[strand=="-",c("start","end"):=list(end,start)]
    if(AC009505.2)
        genes.dt <- rbind(genes.dt,
                          data.table(gene="ENSG00000235522",
                                     symbol = "AC009505.2", 
                                     start = 773.5,
                                     end = -9869.5, strand = "-",
                                     y=-0.2))
    return(list(exons=exons.dt,genes=genes.dt))
}


## genes <- getgenes(CHR,FROM,TO,OFFSET)

genes <- getgenes(11,76680005,76935036,100)



# Plotting figure 5 of peaky paper ----------------------------------------


## 1mm = 2.83465pt
pt <- 1/2.83465

baits <- 886743

plotone <- function(xb,zoom=NULL,pt=0.1) {
    if(!is.null(zoom))
        xb <- xb[dist>=zoom[[1]] & dist <= zoom[[2]], ]
    from <- min(xb$dist) 
    to <- max(xb$dist)
    xb[c(1,nrow(xb)),residual:=0]
    ## peaks.this <- peaks[end-offset > from & start-offset < to,]
    ## rna.this <- rna[end-offset > from & start-offset < to,]
    ## rna.this[,RNA.strand:=ifelse(y>0,"+","-")]
    #tracks.this <- lapply(tracks, function(x) x[end-offset > from & start-offset < to, ])
    ## raw 
    yv <- with(genes$genes[ start>from & end<to,],
               rep(seq(0.03,0.23,by=0.05),length.out=length(y)))
    plots <- list(ggplot(xb,aes(x=dist,y=N+1)) + geom_point(size=pt,col="black",fill="grey30") +
                      ##geom_path(aes(y=B+T)) +
                      scale_y_log10() + ggtitle("Counts (log10 scale)") + ylim(0,max(xb$N+1)),
                  ### resid
                  ggplot(xb, aes(x=(start.prey-offset), xend=(end.prey-offset),yend=residual,#x=dist,
                                y=residual))+
                      geom_point(col="black",fill="grey30",size=pt)+
                      ggtitle("NB residual")+
                      geom_hline(yintercept=0,col="grey"),
                  ### CHiCAGO
                  # ggplot(xb,aes(x=dist,
                  #               #x=(start.prey-offset),xend=(end.prey-offset),yend=chicago,
                  #               ymin=0,ymax=chicago,
                  #               y=chicago)) +
                  #     #geom_ribbon() +
                  #     geom_area(col="black",fill="grey30") +
                  #     geom_area(aes(y=pmin(5,chicago)),col="grey",fill="grey") +
                  #     #                                        geom_point(size=pt) +
                  #     ggtitle("CHiCAGO score (asinh scale)") + geom_hline(yintercept=5,col="grey") +
                  #     scale_y_continuous(trans=asinh_trans(),breaks=c(0,5,20,40),limits=c(0,max(xb$chicago))),
                  ### Credset
                  ## ggplot(xb,aes(xmin=(start.prey-offset),xmax=(end.prey-offset),ymin=0,ymax=h)) + geom_rect(fill="grey20",col="grey20") + ggtitle("MPPC finemap"),
                  ### MPPC
                  
                  # ggplot() +
                  #     geom_rect(aes(xmin=(start-offset),xmax=(end-offset)),
                  #               ymin=0,ymax=max(xb$mppc)^0.5,fill="violet",col="violet",
                  #               data=chromhmm) +
                  #     ## geom_point(aes(x=(mid.prey-offset), y=mppc),
                  #     ##            size=pt,data=xb) +
                  #     geom_area(aes(x=(mid.prey-offset), y=mppc),
                  #               col="black",fill="grey30",
                  #               data=xb) +
                  #     ggtitle("MPPC (sqrt scale)") +
                  #     scale_y_sqrt()
                  # ,
                  ### bamfiles
                  # ac27=ggplot(tracks.this$ac27,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("K3K27ac"),
                  # me1=ggplot(tracks.this$me1,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("H3K4me1"),
                  # me3=ggplot(tracks.this$me3,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("H3K4me3"),
                  ### CHROMHMM - B&W
                  ##ggplot(chromhmm, aes(xmin=(start-offset),xmax=(end-offset),ymin=0,ymax=1,fill=state,col=state)) + geom_rect()  +
                  ## chromhmm=ggplot(chromhmm, aes(xmin=(start-offset),xmax=(end-offset),ymin=0,ymax=1,fill=Chromatin)) + geom_rect(fill="grey20",col="grey20")  + ggtitle("active chromatin"),
                  ### RNAseq
                  
                  ## p.rna=ggplot(rna.this,aes(xmin=start-offset,xmax=end-offset,
                  ##                           ymax=sign(y) * log2(abs(y)),ymin=0)) +
                  ##                           ## fill=RNA.strand,col=RNA.strand)) +
                  ##   geom_rect() +
                  ##   theme(legend.position="none") +
                  ##   ggtitle("RNAseq (log2 counts)"), #+
                  ## scale_colour_manual("Genes (strand)",values=c("grey20","darkslateblue")) +
                  ## guides(colour=FALSE) +
                  ## scale_fill_manual("Genes (strand)",values=c("grey20","darkslateblue")) +
                  ## guides(fill=FALSE)
                  
                  ### genes
                  
                  g=ggplot(genes$genes[ start>from & end<to,],
                           aes(x=start,xend=end,col=strand))+
                      geom_segment(aes(y=sign(y) * yv, #rep(seq(0,0.15,by=0.05),length.out=length(y)),
                                       yend=sign(y) * yv), #rep(seq(0,0.15,by=0.05),length.out=length(y))),
                                   arrow = arrow(length = unit(0.1,"cm")),size=1) +
                      geom_text(aes(x=start,y=sign(y) * yv + sign(y) * 0.02, #rep(seq(0.05,0.2,by=0.05),length.out=length(y)),
                                    label=symbol),hjust=1,size=12 * pt) +
                      ## geom_rect(aes(xmin=start,xmax=end,ymin=y-0.03, ymax=y+0.03,col=strand),
                      ##           data=genes$exons[ start>from & end<to,]) +
                      scale_y_continuous(limits=c(-0.3,0.3)) +
                      scale_colour_manual("Genes (strand)",values=rev(c("grey20","darkslateblue"))) + guides(colour=FALSE)
    )
    ## add bait locations
    bloc <- xb[preyID %in% baits,]$dist %>% unique()%>%setdiff(.,0)
    plots <- lapply(plots, function(p) {
        p + #geom_vline(xintercept=bloc,col="grey",linetype="dashed",alpha=0.5) +
            #background_grid() +
            geom_vline(xintercept=0,col="Firebrick1",linetype="dashed") +
            xlim(from,to) + 
            theme(#legend.position="top",
                axis.text.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()
            )
    })
    plots[-c(1:4)] %<>% lapply(., function(p) {
        p + theme(axis.text.y=element_blank(), #scale_y_continuous(breaks=NULL) 
                  axis.line.y = element_blank(), 
                  axis.ticks.y = element_blank())
    })
    ## turn x axis back on for genes
    plots[[length(plots)]] <- plots[[length(plots)]] +
        scale_x_continuous("Distance from bait (bp)"## , trans="kb"
                           ,limits=c(from,to)) + theme(axis.text.x = element_text(),axis.title.x=element_text()) 
    ## make plot
    plot_grid(plotlist=plots,ncol=1,align="v",
              rel_heights=c(rep(1,length(plots)-1),2))
}
## if(interactive())
##     plotone(xb,zoom=TODO[[g]]$zoom[[i]])

#zoom = P[abs(P$dist)<1e6]


b2g <- data.table::fread("supplementary_files/bait2gene.csv")
head(b2g)
baits <- b2g$baitID

hind <- data.table::fread("supplementary_files/hind-position.csv")
hind[,hindID:=as.integer(hindID)]


xp <- P
xp <- addstartend(xp)


plotone(xp)


plotone(P,zoom=zoom[[g]]$zoom[[i]])
