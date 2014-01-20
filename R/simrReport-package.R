#' simrReport
#'
#' @name simrReport
#' @docType package
NULL

maPlot.simr <- function(fit,coef,p.value=0.05,log2FC=log2(2),reportDir="./reports"){
    exp <- paste(colnames(fit$design)[1],"vs",colnames(fit$design)[coef],sep="_")
    
    edgeRLRT <- glmLRT(fit,coef=coef)
    isAff <- p.adjust(edgeRLRT$table$PV,"BH") < p.value
    isUp <- edgeRLRT$table$logFC >  log2FC
    isDw <- edgeRLRT$table$logFC <  -log2FC

    toPlot <- data.frame(M=edgeRLRT$table$logFC,
                         A=edgeRLRT$A)
    toPlot$side[isAff & isUp] <- 'Up'
    toPlot$side[isAff & isDw] <- 'Down'
    toPlot$side[!((isAff & isDw)|(isAff & isUp))] <- 'NAff'
    
    xlim <- range(toPlot$A[is.finite(toPlot$A)])
    ylim <- range(toPlot$M[is.finite(toPlot$M)])
    ymax <- max(abs(ylim))*1.1
    
    toPlot <- toPlot[!is.na(toPlot$M),]
    toPlot$A[toPlot$A == -Inf] <- min(xlim)-0.5
    toPlot$M[is.infinite(toPlot$M)] <- sign(toPlot$M[is.infinite(toPlot$M)])*ymax

    
    p <- ggplot(toPlot,aes(A,M,color=side))+stat_binhex(bins=125)+
        scale_fill_gradient(low='royalblue3',high='red')+
            labs(title=exp)
    
    img.d <- file.path(reportDir,'images')
    dir.create(img.d,FALSE,TRUE)

    pdf.f <- file.path('images',paste0(exp,'_MAplot.pdf'))
    png.f <- file.path('images',paste0(exp,'_MAplot.png'))

    pdf(file.path(reportDir,pdf.f))
    print(p)
    dev.off()
    png(file.path(reportDir,png.f))
    print(p)
    dev.off()
    return(c(pdf=pdf.f,png=png.f))
}

plotClustering <- function(SE,method=c('spearman','kendall','pearson','euclidian')){
    method <- match.arg(method)
    ## Given a SummarizedExperiment object,
    ## use edgeR to compute a normalized RPKM value
    ## for genes with at least 2 sample with cpm > 5
    dge <- DGEList(assay(SE), group=SE$treatment)
    dge <- calcNormFactors(dge)
    cpm.d <- cpm(dge)
    dge <- dge[ rowSums(cpm.d > 5) >=2, ]
    ## Use the log10 rpkm value
    d <- log10(rpkm(dge, sum(width(SE[rowSums(cpm.d > 5) >=2,]))))
    ## Clustering the data
    if(method == 'euclidian'){
        dd <- d
        r <- range(dd[is.finite(dd)])
        dd[!is.finite(dd) & dd>0] <- r[2]+1
        dd[!is.finite(dd) & dd<0] <- r[1]-1
        x.cluster = hclust(dist(t(dd)))
    } else {
        x.cluster <- hclust(as.dist(1-cor(d,meth=method)))
    }
    ##################################################
    ### Let draw a figure
    ##################################################
    ## Get the color palette from 0 to the max value
    ##zlim <- max(abs(d[is.finite(d)]))
    zlim <- range(d[is.finite(d)])
    cols <- colorRampPalette(c("black","yellow"))(256)
    ## break the region to drawing reagion
    layout(matrix(c(1:3),ncol=1,byrow=TRUE)
           ,heights=c(15,80,10))
    ## Draw the sample dendrogram
    op <- par(mar = c(0,2,4,2))
    plot(as.dendrogram(x.cluster),
         axes=TRUE,
         yaxt='s',
         yaxs='i',
         xaxt='n',
         xaxs='i',
         horiz=FALSE,
         leaflab='none',
         main=sub("(\\w)","\\U\\1",paste(method,"corelation of log10(RPKM)"),perl=TRUE)
         )
    ## Draw the RPKM heatmap
    ## Columns sorted based on the sample clustering
    ## Gene sorted on the overal gene expression level
    par(mar = c(4.5,2,0.2,2))
    image(t(d[order(rowMeans(d)),x.cluster$order]),
          col=cols,
          zlim=zlim,
          xaxt='n',
          yaxt='n'
          )
    ## Add the labels
    labs <- colnames(dge$counts)[x.cluster$order]
    at <- seq(0,1,1/(length(labs)-1))
    axis(side=1,at=at, labels = FALSE)
    text(x=at,y=0-0.025,labels=labs,xpd=TRUE,srt=45,adj=1)
    ## Draw the heatmap scale
    par(mar = c(3,20,0.2,20))
    image(matrix(seq(zlim[1],zlim[2],length=100)),
          col=cols,
          xaxt='n',
          yaxt='n')
    axis(side=1,at=c(0,1),labels=format(zlim,digits=2))
    par(op)
    return(d[order(rowMeans(d)),x.cluster$order])
}
   
##################################################
### Some helper function for publishing tables
##################################################

## This Function adds links to Ensembl
addExternalLink <- function(object,linkOut=NULL,...){
    if('ensembl_gene_id' %in% colnames(object)){
        ids <- c('ensembl_gene_id','external_gene_id')
    }else{
        ids <- c('ID','ID')
    }
    if(!(is.null(linkOut) | is.na(linkOut))){
        object$ExtLink <- hwrite(as.character(object[,ids[2]]), 
                                 link = paste0(linkOut,
                                     as.character(object[,ids[1]])), table = FALSE)
    }
    return(object)
}

addIGVLink <- function(object,gnModel,...){
    if(class(gnModel) == 'GRangesList'){
        gene.ids <- names(gnModel)
        gnModel <- unlist(gnModel)
        names(gnModel) <- gene.ids
    }
    ids <- ifelse('ensembl_gene_id' %in% colnames(object),'ensembl_gene_id','ID')
    gnModel <- gnModel[object[,ids]]
    IGVlinks <- sprintf('http://localhost:60151/goto?locus=%s:%s-%s',
                        seqnames(gnModel),
                        start(gnModel),
                        end(gnModel)
                        )
    object$IGV <- hwrite("IGV",link = IGVlinks, table = FALSE)
    return(object)
}

##This function reorganizes the df
cleanUpDf <- function(object, ...){
    clean.names <- c("ExtLink"="External link",
                     "ID"="ID",
                     "description"="Description",
                     "gene_biotype"="BioType",
                     "logFC"="log2(Fold Change)",
                     "Image"="Expression",
                     "IGV"="To IGV")
    clean.names <- clean.names[names(clean.names) %in% colnames(object)]
    object <- object[,names(clean.names)]
    colnames(object) <- clean.names
    return(object)
}

writeIGVsession <- function(fit,coef,expDesign,IGVgenome,bwDir="./bigwig",reportsRoot="./reports",extraTracks=NULL) {
    
    treat <- colnames(fit$design)[coef]
    
    cntl.bam <- expDesign$files[match(rownames(fit$design)[rowSums(fit$design) == 1],rownames(expDesign))]
    exp.bam <- expDesign$files[match(rownames(fit$design)[fit$design[,treat]==1],rownames(expDesign))]

    cntl <- sub("\\.bam$","",basename(cntl.bam))
    exp <- sub("\\.bam$","",basename(exp.bam))

    ## Recover the bigwig files
    bigwigs <- list.files(bwDir,"\\.bw$",full.names=TRUE)
    bigwigs.loc <- as.vector(sapply(c(cntl,exp),function(x){bw <- grep(x,bigwigs,value=TRUE)
                                                            bw[order(bw,decreasing=TRUE)]
                                                        }))
    
    extraTrk.path = ifelse(!(is.null(extraTracks)|is.na(extraTracks)),
        paste(list.files(extraTracks,pattern='\\.bed$',full=TRUE),collapse=" "),
        "")
    
    IGV.session <- paste("igv_session",
                         "--relative",reportsRoot,
                         "--global",paste0("genome=",IGVgenome),
                         paste(bigwigs.loc,collapse=" "),
                         extraTrk.path
                         )
    IGV.session.file <- paste0(paste(treat,collapse="_"),"_DGE_sesssion.xml")
    write(system(IGV.session,intern=TRUE),file=file.path(reportsRoot,IGV.session.file))
    return(IGV.session.file)
}

### This is one heck of an ugly hack...
makeIGVSessionLink <- function(fit,
                               coef,
                               expDesign,
                               IGVgenome,
                               bwDir="./bigwig",
                               reportsRoot="./reportsDir",
                               extraTracks=NULL,
                               serverRoot=""){
    
    ## Normalize the diffrent paths to their absolute location
    bwDir <- normalizePath(bwDir)
    reportsRoot <- normalizePath(reportsRoot)
    if(!(is.null(extraTracks) | is.na(extraTracks)))
        extraTracks <- normalizePath(extraTracks)
    
    treat <- colnames(fit$design)[coef]

    IGV.session.file <- writeIGVsession(fit=fit,
                                        coef=coef,
                                        expDesign=expDesign,
                                        IGVgenome=IGVgenome,
                                        bwDir=bwDir,
                                        extraTracks=extraTracks,
                                        reportsRoot=reportsRoot)
    
    wd.root <- gsub(paste0(Sys.getenv()['HOME'],'/'),'',getwd())
    
    URI.base <- sub("(.+://).+","\\1",serverRoot)
    serverBaseDir <- sub(".+://(.+)","\\1",serverRoot)
    
    link2server <- paste(serverBaseDir,wd.root,basename(reportsRoot),IGV.session.file,sep='/')
    link2server <- gsub("//","/",link2server)

    paste0("http://localhost:60151/load?file=",URI.base,link2server)
}

gene2name <- function(biomart,martDataset,...){
    mart <- useMart(biomart,martDataset)
    atts <- c('ensembl_gene_id',
              'ensembl_transcript_id',
              'external_gene_id',
              'description',
              'gene_biotype')
    gene2name <-  getBM(attributes=atts,
                        filters = '',
                        values ='',
                        mart=mart)
    gene2name$description <- sub("\\s\\[Source.+","",gene2name$description)
    gene2name
}

getEdgeRdata <- function(fit,gene2name=NULL,coef=2,FC=2,p.val=0.01){
    edgeRLRT <- glmLRT(fit,coef=coef)
    df <- cbind(edgeRLRT$table[,c('logFC','PValue')],adj.p=p.adjust(edgeRLRT$table$PValue,method="BH"))
    if(!is.null(gene2name)){
        df <- cbind(gene2name[match(rownames(df),gene2name$ensembl_gene),],df)
    } else {
        df$ID <- rownames(fit$counts)
        rownames(df) <- rownames(fit$count)
    }
    df$isAff <- df$adj.p < p.val & abs(df$logFC)>=log2(FC)
    return(df)
}

writeCSVfile <- function(treat,df,reportsRoot="./reports"){
    csv.file <- paste0(treat,"_DGE.csv")
    write.csv(df,file=file.path(reportsRoot,csv.file))
    return(csv.file)
}

addCPMPlots <- function(edgeRdata,fit,coef=2,reportsRoot="./reports",
                        nCores=ifelse(!is.na(detectCores()),detectCores()/2,2L)){
    
    if(sum(edgeRdata$isAff)==0){return(NULL)}
    edgeRdata <- data.frame(edgeRdata[edgeRdata$isAff,])
    
    control <- which(rowSums(fit$design)==1)
    treatments <- which(fit$design[,coef]==1)
    countData <- cpm(fit$counts[,c(control,treatments)])

    dir.create(file.path(reportsRoot,"images",colnames(fit$design)[coef]),FALSE,TRUE)

    genes <- edgeRdata[,ifelse('ensembl_gene_id' %in% colnames(edgeRdata),'ensembl_gene_id','ID')]

    ## Lattice seems to have issues with open printing drivers...
    if(!is.null(dev.list())) sapply(dev.list(),dev.off)

    images <- as.data.frame(do.call(rbind,mclapply(genes,function(gene){
        treat <- factor(rep(colnames(fit$design)[c(1,coef)],sapply(list(control,treatments),length)))
        treat <- relevel(treat,colnames(fit$design)[1])
        data <- data.frame(treat=treat,cpm=countData[gene,])
        gene.name <- ifelse('ensembl_gene_id' %in% colnames(edgeRdata),
                            edgeRdata[edgeRdata$ensembl_gene_id == gene,'external_gene_id'],
                            gene)
        p <- stripplot(cpm~treat,
                       data,
                       main=gene.name,
                       jitter.data=TRUE,
                       pch=20,
                       cex=1.25,
                       bg='blue')
        base <- file.path("images",colnames(fit$design)[coef],paste0(gene,"_stripplot"))
        png <- paste0(base,".png")
        pdf <- paste0(base,".pdf")
        png(file.path(reportsRoot,png),150,150)
        print(p)
        dev.off()
        pdf(file.path(reportsRoot,pdf))
        print(p)
        dev.off()
        return(c(png=png,pdf=pdf))
    },mc.cores=nCores)))
    
    edgeRdata$Image <- hwriteImage(images$png,
                                   link=images$pdf, 
                                   table=FALSE,
                                   width=150)
    return(edgeRdata)
}

publishFit <- function(fit,
                       expDesign,
                       gnModel,
                       htmlRep,
                       biomart='ensembl',
                       martDataset='dmelanogaster_gene_ensembl',
                       FC=2,
                       p.val=0.01,
                       linkOut=NULL,
                       bwDir='./reports/bigwig',
                       reportsRoot='./reports',
                       IGVgenome="dmel_5.73",
                       serverRoot="",
                       extraTracks=NULL,
                       nCores=4L){
    ## Opening remarks
    publish(hwrite(paste("A gene is considered significantly change",
                             "if it as a p value lower then",
                             p.val,
                             "and is up or down at least",
                             FC,"fold")),htmlRep)
    
    if(!is.null(biomart) & !is.null(martDataset)){
        gene2name <- gene2name(biomart,martDataset)
    } else {
        gene2name <- NULL
    }
    
    ## Rolling over the different contrast coeficient
    sapply(2:ncol(fit$design),function(coef){
        ## Get the DGE for the coeficient under scrutiny
        treat <- colnames(fit$design)[coef]
        
        edgeRdata <- getEdgeRdata(fit,
                                  coef,
                                  gene2name=gene2name,
                                  FC=FC,
                                  p.val=p.val
                                  )
        
        csv.file <- writeCSVfile(treat,edgeRdata,reportsRoot)
        
        edgeRdata <- addCPMPlots(edgeRdata,
                                 fit,
                                 coef,
                                 reportsRoot=reportsRoot,
                                 nCores=nCores
                                 )
        
        IGVlink <- makeIGVSessionLink(fit=fit,
                                      coef=coef,
                                      expDesign=expDesign,
                                      bwDir=bwDir,
                                      IGVgenome=IGVgenome,
                                      extraTracks=extraTracks,
                                      reportsRoot=reportsRoot,
                                      serverRoot=serverRoot
                                      )

   ################################################## 
   ### Publish to ReportingTools the different object created
   ##################################################
        publish(hwrite(paste('Number of differentially regulated genes in',treat), heading=3),htmlRep)

        if(!is.null(edgeRdata)){
            publish(data.frame("Up regulated"  =sum(edgeRdata$logFC>0 & edgeRdata$isAff),
                               "Down Regulated"=sum(edgeRdata$logFC<0 & edgeRdata$isAff)),htmlRep)
        }else{
            publish(data.frame("Up regulated"  =0,
                               "Down Regulated"=0),htmlRep)
        }

        publish(hwrite(paste('List of differentially regulated genes in',treat), heading=3),htmlRep)
        
        IGVhtml <- hwrite(paste("Click here to open an IGV session with the",
                                colnames(fit$design)[1],"and",paste(treat,collapse=" "),
                                "relative coverages (relative to aligned reads)"),link=IGVlink)
        
        publish(IGVhtml,heading=3,htmlRep)
        
        if(!is.null(edgeRdata)){
            publish(edgeRdata, htmlRep,
                    gnModel=gnModel,
                    linkOut=linkOut,
                    .modifyDF = list(addExternalLink,addIGVLink,cleanUpDf))
        } else {
            publish(hwrite("No genes are differentially regulated",heading=3),htmlRep)
        }
        publish(Link(paste("Click here to download the all of the DGE results for",treat),csv.file),htmlRep)

        return(NULL)
  })

}
