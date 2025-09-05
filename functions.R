renameScores <- function(st, rmLoc=FALSE, rmRaw=TRUE, rmSig=FALSE){
  if(is.data.frame(st)){
    if(!is.null(st$score))
      st$score <- renameScores(st$score, rmLoc=rmLoc, rmRaw=rmRaw, rmSig=rmSig)
    if(!is.null(st$method))
      st$method <- renameScores(st$method, rmLoc=rmLoc, rmRaw=rmRaw, rmSig=rmSig)
    return(st)
  }
  score <- factor(st)
  levels(score) <- gsub("^padj\\.","",levels(score))
  levels(score) <- gsub("Global", "glb", levels(score), ignore.case = TRUE)
  levels(score) <- gsub("Local", "loc", levels(score), ignore.case = TRUE)
  if(rmLoc) levels(score) <- gsub("\\.loc|\\.glb","",levels(score))
  if(rmRaw) levels(score) <- gsub("\\.raw","",levels(score))
  if(rmSig) levels(score) <- gsub("sig\\.","",levels(score))
  score
}

  
.bbhwAll <- function(pbDEA, bulkDEA, pb, verbose=FALSE){
  g <- expand.grid(bin.method=c("PAS","combined","asNA","sig"),
                   correction.method=c("gBH.LSL","binwise","IHW","gBH.TST"),
                   local=c(TRUE,FALSE))
  gi <- intersect(unique(pbDEA$gene), row.names(pb))
  pbDEA <- pbDEA[which(pbDEA$gene %in% gi),]
  pbDEA$ID <- paste0("H",seq_len(nrow(pbDEA)))
  for(i in seq_len(nrow(g))){
    loc <- as.logical(g[i,3])
    name <- paste("padj", g[i,1], gsub("gBH\\.","",g[i,2]),
                  ifelse(loc,"loc","glb"), sep=".")
    if(verbose) print(name)
    res <- tryCatch(
      bbhw(pbDEA, bulkDEA, pb=pb, bin.method=as.character(g[i,1]),
                verbose=verbose, local=loc,
                correction.method=as.character(g[i,2])),
      error=function(e){
        message(paste(name, " failed"))
        return(NULL)
      })
    if(!is.null(res)){
      row.names(res) <- res$ID
      pbDEA[[name]] <- res[pbDEA$ID, "padj"]
    }
  }
  if(TRUE){
    g <- g[which(g[,1]=="sig"),]
    for(i in seq_len(nrow(g))){
      loc <- as.logical(g[i,3])
      name <- paste("padj", g[i,1], gsub("gBH\\.","",g[i,2]),
                    ifelse(loc,"loc","glb"), "raw", sep=".")
      if(verbose) print(name)
      res <- bbhw(pbDEA, bulkDEA, bin.method=as.character(g[i,1]),
                  verbose=verbose, local=loc, useSign = FALSE,
                  correction.method=as.character(g[i,2]))
      row.names(res) <- res$ID
      pbDEA[[name]] <- res[pbDEA$ID, "padj"]
    }
  }
  pbDEA
}


getStats <- function(sl, truth, celltype=rep(1L, length(truth)), roundNd=NULL, noRankAt1=TRUE){
  dplyr::bind_rows(lapply(split(seq_along(truth), celltype), FUN=\(i){
    sl <- sl[i,,drop=FALSE]
    truth <- truth[i]
    if(sum(truth)<3) return(data.frame(nominal=numeric(0), recall=numeric(0), fdr=numeric(0)))
    dplyr::bind_rows(lapply(as.data.frame(sl), \(x){
      o <- order(x)
      if(noRankAt1) o <- intersect(order(x),which(x<1))
      d <- data.frame(nominal=c(x[o],1), label=c(truth[o],NA),
                      recall=c(cumsum(truth[o])/sum(truth),1),
                      fdr=c(cumsum(!truth[o])/seq_along(x[o]), 1-sum(truth)/length(truth)))
      if(!is.null(roundNd)){
        w <- which(d$nominal>0.25)
        d$recall <- round(d$recall, roundNd)
        d$fdr <- round(d$fdr, roundNd)
        d <- d[rev(seq_len(nrow(d))),]
        d <- d[!duplicated(d[,3:4]),]
      }
      d
    }), .id="score")
  }), .id="celltype")
}

plotPR <- function(st, ths=c(0.05,0.1,0.25), facet_scores=TRUE, sqrty=FALSE, noLine=FALSE, sqrtx=TRUE, ...){
  if(!is.null(ths)){
    thsd <- dplyr::bind_rows(lapply(split(st, st[,c("celltype","score")]), FUN=function(x){
      x[unlist(lapply(ths, FUN=\(th){
        w <- which(x$nominal<=th)
        if(length(w)==0) return(c())
        max(w)
      })),]
    }))
    thsd$th <- factor(sapply(thsd$nominal, FUN=function(x) ths[which.min(abs(x-ths))]))
  }
  if(!facet_scores){
    p <- ggplot(st, aes(recall, 1-fdr, colour=score, linetype=score)) + 
      facet_wrap(~celltype, ...) + labs(shape="Nominal\nFDR\nthreshold")
    thaes <- aes(shape=th)
  }else{
    p <- ggplot(st, aes(recall, 1-fdr)) + facet_grid(celltype~score, ...) + 
      labs(shape="Nominal\nFDR\nthreshold", colour="Nominal\nFDR\nthreshold")
    thaes <- aes(colour=th, shape=th)
  }
  if(noLine){
    p <- p + #geom_path(size=0.5, alpha=0.1) + 
      ggrastr::geom_point_rast(size=0.03)
  }else{
    p <- p + geom_path()
  }
  if(!is.null(ths)) p <- p + geom_point(thaes, data=thsd, size=2.5)
  if(sqrty) p <- p + scale_y_sqrt(breaks=scales::pretty_breaks(4))
  if(sqrtx) p <- p + scale_x_sqrt(breaks=c(0,0.2,0.5,1), labels=c("0",".2",".5","1"))
  p + theme_bw()
}


tprPlot <- function(st, ths=c(0.05,0.1,0.25), sqrty=TRUE, sqrtx=FALSE, mergeSeeds=TRUE,
                    legend.inside=TRUE, nrow=2, leg.spacing=0.35, ...){
  st <- st[!grepl("Excitatory",st$celltype),]
  if(!is.null(st$Dataset))
    st$celltype <- factor(st$celltype, unique(st$celltype[order(st$Dataset)]))
  if(is.null(st$seed) || mergeSeeds){
    st$group <- st$score
  }else{
    st$group <- paste(st$score, st$seed)
  }
  if(!is.null(ths)){
    thsd <- dplyr::bind_rows(lapply(split(st, st[,c("celltype","group")]), FUN=function(x){
      x[unlist(lapply(ths, FUN=\(th){
        w <- which(x$nominal<=th)
        if(length(w)==0) return(c())
        max(w)
      })),]
    }))
    thsd$th <- factor(sapply(thsd$nominal, FUN=function(x) ths[which.min(abs(x-ths))]))
  }
  p <- ggplot(st, aes(fdr, recall, colour=score, linetype=score, group=group)) + 
    facet_wrap(~celltype, nrow=nrow, ...) + theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_vline(xintercept=ths, linetype="dashed", colour="darkgrey") +
    ggrastr::geom_point_rast(size=0.03)
  if(!is.null(ths)) p <- p + geom_point(aes(shape=th), data=thsd, size=3)
  if(sqrty) p <- p + scale_y_sqrt(breaks=scales::pretty_breaks(4))
  if(sqrtx){
    p <- p + scale_x_sqrt(breaks=c(0,0.2,0.5,1), labels=c("0",".2",".5","1"))
  }else{
    p <- p + scale_x_continuous(breaks=c(0,0.25,.5,.75,1),
                                labels=c("0",".25",".5",".75","1"))
  }
  if(legend.inside)
    p <- p + theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
      guides(colour=guide_legend(keyheight=leg.spacing, default.unit="cm"),
             shape=guide_legend(keyheight=leg.spacing, default.unit="cm")) +
    theme(legend.spacing.y = unit(0, 'cm'), legend.box.spacing = unit(0,"pt"), 
          legend.box.margin=margin(0,0,0,0), legend.title=element_text(size=10))
  p + labs(y="True Positive Rate (TPR, i.e. recall)", colour="Method",
       x="False Discovery Rate (FDR)", shape="Nominal\nFDR threshold")
}

mergeSts <- function(sim, ms){
  st <- dplyr::bind_rows(list("simulation"=sim, "MS data"=ms), .id="Dataset")
  st <- st[!grepl("Excitatory",st$celltype),]
  st$celltype <- factor(paste(st$Dataset, st$celltype, sep="\n"))
  st
}

plotThStats <- function(x){
  x2 <- x[x$variable=="F1" & !grepl("Unaffected",x$celltype),]
  x2$F1 <- x2$value
  x <- merge(x[x$variable %in% c("precision","recall"),], x2[,c("method","celltype","threshold","F1")])
  ggplot(x, aes(method, value, fill=F1)) + geom_col() + 
    facet_grid(variable~celltype, scales="free_y") + theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.title.x=element_blank()) + 
    scale_y_continuous(breaks=scales::pretty_breaks(2)) +
    scale_fill_viridis_c()
}

plotThStats2 <- function(x, type=c("mean","median")){
  type <- match.arg(type)
  x <- as.data.frame(x[!grepl("Unaffected",x$celltype),])
  x <- x[grep(type, x$variable),]
  x$variable <- gsub(paste0(type,"."), "", x$variable, fixed=TRUE)
  x$variable <- factor(x$variable, c("precision","recall", "F1"))
  ggplot(x, aes(method, value, fill=method)) + geom_col() + 
    geom_segment(aes(y=value-var, yend=value+var)) + 
    facet_grid(variable~celltype, scales = "free_y") +
    theme_bw() + scale_y_continuous(breaks=scales::pretty_breaks(2)) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.title.x=element_blank(), legend.position = "none")
}


computeMetrics <- function(ssres, truth=NULL, th=0.1, scores=NULL, relative=FALSE) {
  library(data.table)
  stopifnot(!is.null(ssres[[1]]$isDEG) || !is.null(truth))
  .fun <- \(m, truth, th, scores=NULL, relative=FALSE) {
    if(is.null(scores))
      scores <- grep("^FDR|^padj\\.", colnames(m), value=TRUE)
    if(is.null(m$isDEG)){
      m <- merge(m, truth[,c("celltype","gene","isDEG")], 
                 by=c("celltype","gene"))
    }
    
    m <- m[which(!is.na(m$isDEG)),]
    m <- m[order(m$celltype, m$PValue),]
    ss <- dplyr::bind_rows(lapply(setNames(scores,scores), FUN=function(x){
      dplyr::bind_rows(lapply(split(m, m$celltype), \(m){
        m2 <- m[which(m[[x]]<th),]
        ret <- data.frame(TP=sum(m2$isDEG), FP=sum(!m2$isDEG))
        ret$precision=as.numeric(ret[1]/sum(ret)) 
        ret$recall=as.numeric(ret[1]/sum(m$isDEG))
        ret$F1 <- 2/(1/ret$precision+1/ret$recall)
        ret
      }), .id="celltype")
    }), .id="method")
    if(!relative) return(ss)
    baseline <- ss[which(ss$method=="FDR.local"),]
    ff <- c("precision","recall","F1")
    row.names(baseline) <- baseline$celltype
    ss[,ff] <- ss[,ff]-baseline[as.character(ss$celltype), ff]
    ss
  }
  ss <- lapply(ssres, \(x) .fun(x, truth, th, scores=scores, relative=relative))
  dt <- rbindlist(ss, idcol = "seed")
  dt[is.na(F1), F1 := 0]
  td <- dt[, .(
    mean.precision = mean(precision, na.rm=TRUE),
    sem.precision = sd(precision, na.rm=TRUE)/sqrt(length(precision)),
    median.precision = median(precision, na.rm=TRUE),
    mad.precision = median(abs(precision-median(precision, na.rm=TRUE))),
    mean.recall = mean(recall),
    sem.recall = sd(recall)/sqrt(length(recall)),
    median.recall = median(recall),
    mad.recall = median(abs(recall-median(recall, na.rm=TRUE))),
    mean.F1 = mean(F1),
    sem.F1 = sd(F1)/sqrt(length(F1)),
    median.F1 = median(F1),
    mad.F1 = median(abs(F1-median(F1, na.rm=TRUE))),
    TP    = as.numeric(median(TP)),
    mean.FP    = as.numeric(mean(FP))
  ), by = .(method, celltype)]
  
  td1 <- melt(
    td,
    measure.vars = paste(rep(c("mean","median"), each=3), c("precision", "recall", "F1"), sep="."),  
    id.vars = c("method", "celltype")
  )
  td2 <- melt(
    td,
    measure.vars = paste(rep(c("sem","mad"), each=3), c("precision", "recall", "F1"), sep="."),  
    id.vars = c("method", "celltype"),
    value.name = "var"
  )
  td1$var <- td2$var
  td1 <- td1[!grepl("Unaffected", td1$celltype),]
  data.table(td1, threshold=th)
}
