library(Rsubread)
require(dplyr)
require(foreach)
require(parallel)
require(doParallel)
registerDoParallel(40)

setwd('alignment')

single.end.exp    <- c('4cell.input','ESC.H3K4me3.1','ESC.H3K4me3.2','ESC.input')
output.list       <- system('ls | grep all.bam|grep -v H3K27me3',intern=TRUE)
exp.name          <- gsub(output.list,pattern=".all.bam",replacement="")
paired.end.exp    <- exp.name[ !(exp.name %in% single.end.exp)]          
suffix.vec        <- c('.all.rmdup.bam','.unique.rmdup.bam')

single.end.bam.file <- foreach(x =single.end.exp,.combine='c' ) %do%{
    paste(x,suffix.vec,sep="")
}
paired.end.bam.file <- foreach(x =paired.end.exp,.combine='c' ) %do%{
    paste(x,suffix.vec,sep="")
}


rt.paired.rs <- featureCounts(files      = paired.end.bam.file,
                              annot.ext  = "../annotation/mm10.repbase.gtf",
                              isGTFAnnotationFile = TRUE,
                              GTF.featureType = 'exon',
                              GTF.attrType = 'repeat',
                              useMetaFeatures=TRUE,
                              allowMultiOverlap=TRUE,
                              countMultiMappingReads=TRUE,
                              isPairedEnd=TRUE,
                              requireBothEndsMapped=TRUE,
                              countChimericFragments=FALSE,
                              reportReads = FALSE,
                              nthreads = 16
                             )

pr.paired.rs <- featureCounts(files = paired.end.bam.file,
                              annot.ext="../annotation/mm10.protein.coding.promoter.gtf",
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType = 'exon',
                              GTF.attrType='gene_id',
                              useMetaFeatures=TRUE,
                              allowMultiOverlap=TRUE,
                              countMultiMappingReads=TRUE,
                              isPairedEnd=TRUE,
                              requireBothEndsMapped=TRUE,
                              countChimericFragments=FALSE,
                              reportReads = FALSE,
                              nthreads = 16
                             )
                             
                             
rt.single.rs <- featureCounts(files      = single.end.bam.file,
                              annot.ext  = "../annotation/mm10.repbase.gtf",
                              isGTFAnnotationFile = TRUE,
                              GTF.featureType = 'exon',
                              GTF.attrType = 'repeat',
                              useMetaFeatures=TRUE,
                              allowMultiOverlap=TRUE,
                              countMultiMappingReads=TRUE,
                              isPairedEnd=FALSE,
                              reportReads = FALSE,
                              nthreads = 16
                             )

pr.single.rs <- featureCounts(files = single.end.bam.file,
                              annot.ext="../annotation/mm10.protein.coding.promoter.gtf",
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType = 'exon',
                              GTF.attrType='gene_id',
                              useMetaFeatures=TRUE,
                              allowMultiOverlap=TRUE,
                              countMultiMappingReads=TRUE,
                              isPairedEnd=FALSE,
                              reportReads = FALSE,
                              nthreads = 16
                             )


rt.paired.tag.cnt.matrix <- rt.paired.rs [['counts']]
rt.single.tag.cnt.matrix <- rt.single.rs [['counts']]
rt.paired.tag.cnt.matrix <- rt.paired.tag.cnt.matrix[rownames(rt.single.tag.cnt.matrix),]

pr.paired.tag.cnt.matrix <- pr.paired.rs [['counts']]
pr.single.tag.cnt.matrix <- pr.single.rs [['counts']]
pr.paired.tag.cnt.matrix <- pr.paired.tag.cnt.matrix[rownames(pr.single.tag.cnt.matrix),]



rt.tag.cnt.matrix <- cbind(rt.paired.tag.cnt.matrix,rt.single.tag.cnt.matrix)
pr.tag.cnt.matrix <- cbind(pr.paired.tag.cnt.matrix,pr.single.tag.cnt.matrix)





rownames(pr.tag.cnt.matrix) <- sapply(rownames(pr.tag.cnt.matrix), 
                                      function(x)  { l <- strsplit(x,split = '\\.') %>% unlist  
                                                     l[1]
                                                    }
                                      )

tag.cnt.matrix           <- rbind(rt.tag.cnt.matrix,pr.tag.cnt.matrix)
colnames(tag.cnt.matrix) <- gsub(colnames(tag.cnt.matrix),pattern="X",replacement="")
tag.cnt.matrix           <- tag.cnt.matrix[,sort(colnames(tag.cnt.matrix))]


####################

# rt.paired.rs <- featureCounts(files      = paired.end.bam.file,
#                               annot.ext  = "../annotation/mm10.repbase.gtf",
#                               isGTFAnnotationFile = TRUE,
#                               GTF.featureType = 'exon',
#                               GTF.attrType = 'gene_id',
#                               useMetaFeatures=TRUE,
#                               allowMultiOverlap=TRUE,
#                               countMultiMappingReads=TRUE,
#                               isPairedEnd=TRUE,
#                               requireBothEndsMapped=TRUE,
#                               countChimericFragments=FALSE,
#                               reportReads = FALSE,
#                               nthreads = 16
#                              )
#                              
# rt.single.rs <- featureCounts(files      = single.end.bam.file,
#                               annot.ext  = "../annotation/mm10.repbase.gtf",
#                               isGTFAnnotationFile = TRUE,
#                               GTF.featureType = 'exon',
#                               GTF.attrType = 'gene_id',
#                               useMetaFeatures=TRUE,
#                               allowMultiOverlap=TRUE,
#                               countMultiMappingReads=TRUE,
#                               isPairedEnd=FALSE,
#                               reportReads = FALSE,
#                               nthreads = 16
#                              )
#                              
# rt.paired.tag.cnt.matrix    <- rt.paired.rs [['counts']]
# rt.single.tag.cnt.matrix    <- rt.single.rs [['counts']]
# rt.paired.tag.cnt.matrix    <- rt.paired.tag.cnt.matrix[rownames(rt.single.tag.cnt.matrix),]
# rt.tag.cnt.matrix           <- cbind(rt.paired.tag.cnt.matrix,rt.single.tag.cnt.matrix)
# colnames(rt.tag.cnt.matrix) <- gsub(colnames(rt.tag.cnt.matrix),pattern="X",replacement="")
# rt.tag.cnt.matrix           <- rt.tag.cnt.matrix[,sort(colnames(rt.tag.cnt.matrix))]
# 
# if( sum( colnames(rt.tag.cnt.matrix) != colnames(tag.cnt.matrix) ) !=0    ){
#     quit()
# }
# tag.cnt.matrix <- rbind(tag.cnt.matrix,rt.tag.cnt.matrix)


###############
single.end.exp  <- c('4cell.input','ESC.H3K4me3.1','ESC.H3K4me3.2','ESC.input')
suffix.vec      <- c('.all.rmdup.bam','.unique.rmdup.bam')
exp.meta.df     <- foreach(suffix = suffix.vec,.combine='rbind') %do% {
    output.list     <- sprintf(" ls | grep '%s' ",suffix) %>% system(intern = TRUE)
    df1 <- foreach(bam.file = output.list,.combine='rbind') %dopar%{
        cmd.str        <-  sprintf("samtools flagstat %s | head -n 1",bam.file)
        output         <-  system(cmd.str,intern=TRUE)
        tmp            <-  strsplit(x = output,split = "\\s+",perl = TRUE) %>% unlist
        total.tag.cnt  <-  tmp[1] %>% as.integer
        exp.name       <-  gsub(x = bam.file, pattern = suffix,replacement="")
        if((exp.name %in% single.end.exp) == FALSE){
            total.tag.cnt <- total.tag.cnt / 2  # number of fragment
        }
        data.frame(bam.file=bam.file,total.tag.cnt=total.tag.cnt)
    }
    df1
}
rownames(exp.meta.df) <- exp.meta.df$bam.file

###################

save(file='../RData/alignment.RData',list=c('exp.meta.df','tag.cnt.matrix'))










