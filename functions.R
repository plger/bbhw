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



getThStats <- function(ssres, truth, th=0.1, scores=NULL){
  m <- dplyr::bind_rows(ssres, .id="seed")
  if(is.null(scores))
    scores <- grep("^FDR|^padj\\.", colnames(m), value=TRUE)
  if(is.null(truth$celltype)) truth$celltype <- truth$cluster_id
  if(is.null(truth$isDEG)) truth$isDEG <- !is.na(truth$logFC)
  m <- merge(m, truth[,c("celltype","gene","isDEG")], by=c("celltype","gene"), all.x=TRUE)
  m <- m[which(!is.na(m$isDEG)),]
  m <- m[order(m$celltype, m$PValue),]
  scores <- intersect(scores, colnames(ssres[[1]]))
  ss <- dplyr::bind_rows(lapply(setNames(scores,scores), FUN=function(x){
    dplyr::bind_rows(lapply(split(m, m$celltype), \(m){
      m2 <- m[which(m[[x]]<th),]
      ret <- data.frame(TP=sum(m2$isDEG), FP=sum(!m2$isDEG))
      ret$precision=as.numeric(ret[1]/sum(ret))
      ret$recall=as.numeric(ret[1]/sum(m$isDEG))
      ret$calledSig <- (ret$TP + ret$FP)/length(ssres)
      ret$F1 <- 2/(1/ret$precision+1/ret$recall)
      ret
    }), .id="celltype")
  }), .id="method")
  ss <- reshape2::melt(ss, id.vars=c("method", "celltype", "TP", "FP","calledSig","F1"),
                       measure.vars=c("precision", "recall", "F1"))
  ss2 <- ss[ss$variable=="F1",c("method","value")]
  ss2 <- aggregate(ss2[,2], ss2[,1,drop=FALSE], FUN=mean)
  ss2$x[is.na(ss2$x)] <- 0
  ss$method <- factor(ss$method, ss2$method[order(ss2$x)])
  ss
}

plotThStats <- function(x){
  x <- x[x$variable %in% c("precision","recall"),]
  ggplot(x, aes(method, value, fill=F1)) + geom_col() + 
    facet_grid(variable~celltype, scales="free_y") + theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.title.x=element_blank()) + 
    scale_y_continuous(breaks=scales::pretty_breaks(2)) +
    scale_fill_viridis_c()
}
