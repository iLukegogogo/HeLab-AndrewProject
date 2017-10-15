load('RData/alignment.RData')
require(plyr)
require(dplyr)
require(foreach)
require(reshape2)
require(gplots)


get.zscore.2 <- function(m,a,n1,n2){ #n1 for IP, n2 for input
  condition.mean <- 0
  tmp            <- 2^(a-log2(sqrt(n1)) - log2(sqrt(n2)))    
  condition.var  <- (4*(1-tmp)* (1/log(2)) * (1/log(2))) / ((n1+n2)*tmp)
  z.score        <- (m-condition.mean) / sqrt(condition.var)
  z.score
}

get.zscore <- get.zscore.2


#########################################################################################################
# Analyze gao dataset
# Well, I will  correct the batch bias  with lowess
#########################################################################################################


chip.stage.vec <- c('2cell','4cell','8cell','ctrl','ESC','ICM','MIIOocyte','morula','oeKdm5b','siKdm5b','TE','TSC')
tag.cnt.matrix <- tag.cnt.matrix[,grepl(x=colnames(tag.cnt.matrix),pattern='all.rmdup.bam')]
tag.cnt.matrix <- tag.cnt.matrix[1:902,]


gao.enrichment.df <- foreach(chip.stage = chip.stage.vec,.combine='rbind') %do% {
    input.exp    <- sprintf("%s.input",  chip.stage)
    ip.exp       <- sprintf("%s.H3K4me3",chip.stage)
    input.flag   <- grepl(x=colnames(tag.cnt.matrix),pattern=input.exp)
    ip.flag      <- grepl(x=colnames(tag.cnt.matrix),pattern=ip.exp)
    input.matrix <- tag.cnt.matrix[,input.flag]
    ip.matrix    <- tag.cnt.matrix[,ip.flag]
    if(sum(input.flag) >= 2){
        input.tag.cnt <- apply(input.matrix,MARGIN = 1,sum) 
    }else{
        input.tag.cnt <- input.matrix  
    }
    if(sum(ip.flag) >= 2){
      ip.tag.cnt      <- apply(ip.matrix,MARGIN = 1,sum) 
    }else{
      ip.tag.cnt      <- ip.matrix  
    }
    
    
    M                   <-  log2(ip.tag.cnt) - log2(input.tag.cnt)
    A                   <-  (log2(ip.tag.cnt) + log2(input.tag.cnt))/2
    ip.total.tag.cnt    <-  sum(exp.meta.df[colnames(tag.cnt.matrix)[ip.flag]    %>% as.character,'total.tag.cnt'])
    input.total.tag.cnt <-  sum(exp.meta.df[colnames(tag.cnt.matrix)[input.flag] %>% as.character,'total.tag.cnt'])
    bias                <-  log2(ip.total.tag.cnt) - log2(input.total.tag.cnt)
    M                   <-  M - bias
    
    flag                <- ip.tag.cnt>0 & input.tag.cnt>0
    M                   <- M[flag]
    A                   <- A[flag]

    #   lowess fit to correct batch bias for each experiment  
    fit                 <-  loess(data = data.frame(M=M,A=A),M ~ A,span=0.3,degree=2)    
    adjusted.M          <-  residuals(fit) 

    plot(x=A,y=M,main=chip.stage)
    lowess(x=A,y=M) %>% lines
    plot(x=A,  y=adjusted.M,main=chip.stage)
    lowess(x=A,y=adjusted.M) %>% lines
    
    M                   <-  adjusted.M
    df                  <-  data.frame(id = rownames(tag.cnt.matrix)[flag],M=M,A=A,stage = chip.stage,n1=ip.total.tag.cnt,n2=input.total.tag.cnt)
    df                  <-  df[complete.cases(df),]
    df
}

gao.enrichment.df$stage             <- as.character(gao.enrichment.df$stage)
gao.enrichment.df$stage[gao.enrichment.df$stage == 'MIIOocyte'] <- 'oocyte'
gao.enrichment.df                   <- gao.enrichment.df[gao.enrichment.df$stage %in% c('oocyte','2cell','4cell','8cell','morula'),]
gao.enrichment.df$zscore            <- get.zscore(m = gao.enrichment.df$M,  gao.enrichment.df$A,gao.enrichment.df$n1,gao.enrichment.df$n2)
gao.enrichment.df$p.value           <- pnorm(abs(gao.enrichment.df$zscore),lower.tail = FALSE) * 2
tmp                                 <- ddply(gao.enrichment.df,.(id),nrow)
tmp                                 <- tmp[tmp$V1 == 5,]
gao.enrichment.df                   <- gao.enrichment.df[gao.enrichment.df$id %in% as.character(tmp$id),]

#Store chip-seq results into csv files 
foreach(stage = c('oocyte','2cell','4cell','8cell','morula')) %do% {
    flag     <- gao.enrichment.df$stage == stage
    df       <- gao.enrichment.df[flag,]
    df       <- arrange(df,desc(zscore))
    df$fdr   <- p.adjust(df$p.value,method='fdr')
    rs.file  <- sprintf('results/gao.%s.csv',stage)
    df$n1 <- NULL
    df$n2 <- NULL
    write.csv(df,file=rs.file,quote=FALSE)
}




####################################################
# Get gao.zscore.matrix
####################################################
tmp                                 <- gao.enrichment.df[,c('id','stage','zscore')]
gao.zscore.df                       <- dcast(tmp,formula = id~stage,value.var = 'zscore')
rownames(gao.zscore.df)             <- gao.zscore.df$id
gao.zscore.df$id                    <- NULL
gao.zscore.df                       <- gao.zscore.df[,c('oocyte','2cell','4cell','8cell','morula')]
colnames(gao.zscore.df)[1]          <- 'oocyte'
gao.zscore.df                       <- gao.zscore.df[complete.cases(gao.zscore.df),]
gao.zscore.matrix                   <- as.matrix(gao.zscore.df)
rt.var                              <- apply(gao.zscore.matrix,MARGIN=1,var)
rt.var                              <- sort(rt.var,decreasing = TRUE)
gao.zscore.matrix                   <- gao.zscore.matrix[names(rt.var),] # rank all RT family according to variation across 5 stages
write.csv(gao.zscore.matrix,file='results/gao.zscore.matrix.csv',quote=FALSE)


####################################################
#Get gao.M.matrix
####################################################
tmp                            <- gao.enrichment.df[,c('id','stage','M')]
gao.M.df                       <- dcast(tmp,formula = id~stage,value.var = 'M')
rownames(gao.M.df)             <- gao.M.df$id
gao.M.df$id                    <- NULL
gao.M.df                       <- gao.M.df[,c('oocyte','2cell','4cell','8cell','morula')]
colnames(gao.M.df)[1]          <- 'oocyte'
gao.M.df                       <- gao.M.df[complete.cases(gao.M.df),]
gao.M.matrix                   <- as.matrix(gao.M.df)
rt.var                         <- apply(gao.zscore.matrix,MARGIN=1,var)
rt.var                         <- sort(rt.var,decreasing = TRUE)
gao.M.matrix                   <- gao.M.matrix[names(rt.var),] # rank all RT family according to variation across 5 stages
write.csv(gao.M.matrix,file='results/gao.M.matrix.csv',quote=FALSE)



####################################################
#Generate heatmap by heatmap.2
####################################################
scaled.gao.zscore.matrix            <- apply(gao.zscore.matrix,MARGIN=1,scale) %>% t
colnames(scaled.dgao.zscore.matrix) <- c('oocyte','2cell','4cell','8cell','morula')

rt.family.and.type        <- read.table("./annotation/rt.family.and.type.txt", quote="\"", stringsAsFactors=FALSE)
rt.family.and.type$family <- 'LTR'
rt.family.and.type[grepl('SINE',rt.family.and.type$V2),'family'] <- 'SINE'
rt.family.and.type[grepl('LINE',rt.family.and.type$V2),'family'] <- 'LINE'
rt.family.and.type$V2                                            <-  NULL
colnames(rt.family.and.type)[1]                                  <- 'rt.id'
rownames(rt.family.and.type)                                     <- rt.family.and.type$rt.id

rt.family.and.type$color                                         <- '#E41A1C'
rt.family.and.type$color [rt.family.and.type$family == 'SINE']   <- '#4DAF4A'
rt.family.and.type$color [rt.family.and.type$family == 'LTR']    <- '#377EB8'



library(RColorBrewer)
my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
TT <- scaled.gao.zscore.matrix [1:70,]
rs <- heatmap.2( 
            x=as.matrix(TT),
            RowSideColors = rt.family.and.type[rownames(TT) %>% as.character,'color'],
            sepcolor=1,
            Colv = FALSE,
            Rowv = TRUE,
            tracecol = 'NA',
            cexCol = 2,
            cexRow = 0.65,
            dendrogram = 'row',
            margins = c(6,20),
            keysize=1,
            density.info='none',
            scale='none',
            col=my_palette,
            distfun=function(x) as.dist((1-cor(t(x)))/2)
)
legend("bottomright", legend=c('LINE','LTR','SINE') , fill=c('#E41A1C','#377EB8','#4DAF4A'),cex=2)



tmp                        <- rs$rowDendrogram[[1]][[2]]
rs$rowDendrogram[[1]][[2]] <- rs$rowDendrogram[[1]][[1]]
rs$rowDendrogram[[1]][[1]] <- tmp

tmp                             <- rs$rowDendrogram[[1]][[1]][[2]]
rs$rowDendrogram[[1]][[1]][[2]] <- rs$rowDendrogram[[1]][[1]][[1]]
rs$rowDendrogram[[1]][[1]][[1]] <- tmp

tmp                             <- rs$rowDendrogram[[1]][[1]][[2]]
rs$rowDendrogram[[1]][[1]][[2]] <- rs$rowDendrogram[[1]][[1]][[1]]
rs$rowDendrogram[[1]][[1]][[1]] <- tmp

tmp                                  <- rs$rowDendrogram[[1]][[1]][[2]][[1]]
rs$rowDendrogram[[1]][[1]][[2]][[1]] <- rs$rowDendrogram[[1]][[1]][[2]][[2]]
rs$rowDendrogram[[1]][[1]][[2]][[2]] <- tmp


#This is the heatmap andrew wants, clustering by heatmap.2
rs.1 <- heatmap.2( x=TT,
                 Colv = FALSE,
                 RowSideColors = rt.family.and.type[rownames(TT) %>% as.character,'color'],
                 Rowv=rs$rowDendrogram,
                 tracecol = 'NA',
                 cexCol = 2,
                 cexRow = 0.65,
                 dendrogram = 'row',
                 margins = c(6,20),
                 keysize=1,
                 density.info='none',
                 scale='none',
                 col=my_palette,
                 distfun=function(x) as.dist((1-cor(t(x)))/2)
)


##################################################
#Generate heatmap by heatmap.2 (clustered by PCA)
##################################################
pc1 <- prcomp(gao.zscore.matrix)$x[,1] %>% sort
heatmap.2( x=gao.zscore.matrix[names(pc1),],
           Colv = FALSE,
           Rowv = FALSE,
           tracecol = 'NA',
           cexCol = 2,
           cexRow = 0.8,
           dendrogram = 'none',
           margins = c(10,20),
           keysize=1,
           density.info='none',
           scale='none',#
           col=my_palette           
)



#######################
#The interesting RT-gene pairs were picked out by Andrew, they are
# CDK2ap1 - MT2B2
# GTDC1   - MTA_Mm
# Guca1a  - MT2_Mm
# Mtf2    - Mteb
# Pemt    - MT2_Mm
# Rbm25   - MT2_Mm
#
#I generated H3K4me3 modification figure for Pemt and CDK2AP1 and store them in the results/Figure.pptx
#
#######################


