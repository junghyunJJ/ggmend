


geteffect <- function(G1,S1){
  beta <- solve(crossprod(cbind(1,S1),cbind(1,S1)))%*% crossprod(cbind(1,S1),G1)
  return(beta[2,1])
}


BggEffectEst <- function(S1,S2,G1,G2){
  Bs1g1 <- geteffect(G1,S1)
  Bs2g2 <- geteffect(G2,S2)
  G2residual <- G2-S2*Bs2g2
  Bs1g2RESID <- geteffect(G2residual,S1)
  Bg1g2RESID_est_ <- Bs1g2RESID/Bs1g1
  return(abs(Bg1g2RESID_est_))
}


Bggpvalue <- function(S1, S2, G1, G2, Bg1g2, maxperms, threads, seed){
  set.seed(seed)

  if(nchar(maxperms) < 3 ) stop("Please reset maxperms")

  for (i in 2:(nchar(maxperms)-1)) {
    z <- as.numeric(sub("([0-9]).*","\\1",maxperms))
    p <- 10^i*z
    limit <- 5/p*z
    Bg1g2NULL <- c()
    if(.Platform$OS.type == "unix" ){
      set.seed(seed)
      Bg1g2NULL <- unlist(parallel::mclapply(1:p, function(i) {
        permorder <- sample(nrow(S1))
        BggEffectEst(S1[permorder],S2,G1[permorder],G2)
      }, mc.cores = threads))
    }else if(.Platform$OS.type == "windows" && threads > 1){
      stop("Please reset threads")
    }else{
      Bg1g2NULL <- unlist(lapply(1:p, function(i) {
        permorder <- sample(nrow(S1))
        BggEffectEst(S1[permorder],S2,G1[permorder],G2)
      }))
    }
    pval <- 1-length(Bg1g2NULL[Bg1g2NULL<Bg1g2])/length(Bg1g2NULL)
    if (pval > limit) break
  }
  return(pval)
}


Discoverbgg <- function(S1, S2, G1, G2, maxperms, threads, seed){

  Bg1g2 <- BggEffectEst(S1, S2, G1, G2)
  Bg2g1 <- BggEffectEst(S2, S1, G2, G1)

  raw_bgg <- c(Bg1g2,Bg2g1)
  Bgg <- raw_bgg[which.max(raw_bgg)]

  if(which.max(raw_bgg) == 1){
    pval <- Bggpvalue(S1, S2, G1, G2, Bg1g2, maxperms, threads, seed)
  }else{
    pval <- Bggpvalue(S2, S1, G2, G1, Bg2g1, maxperms, threads, seed)
  }

  Bs1g1 <- geteffect(G1,S1)
  Bs2g2 <- geteffect(G2,S2)

  res <- data.frame(Bs1g1 = Bs1g1,
                    Bs2g2 = Bs2g2,
                    GGcor = cor(G1,G2),
                    Bg1g2 = Bg1g2,
                    Bg2g1 = Bg2g1,
                    Bgg = Bgg,
                    pval = pval)
  return(res)
}
#' GGmend (Gene-on-Gene effect estimator using Mendelian randomization)
#'
#' ggmend accurately identifies gene-on-gene regulatory effects utilizing a quartet of two genes and their two \emph{cis}-eQTLs}
#' @param ys A \emph{m} by \emph{n} matrix, where \emph{m} is the number of egenes, and \emph{n} is the number of individuals
#' @param xs A \emph{m} by \emph{n} matrix, where \emph{m} is the number of \emph{cis}-acting snps, and \emph{n} is the number of individuals
#' @param effectth effect size filtering threshold. default = 0
#' @param corth correlation filtering threshold. default = 0
#' @param standardize logical value. If false ys and xs transform the genotype data and the expression data for each gene to have mean 0 and variance 1
#' @param maxperms Number of maximum permutation. default = 1000
#' @param threads The number of cores to use. default = 1
#' @param seed set.seed
#' @keywords ggmend
#' @export
#' @examples
#' data(data_ggmend)
#' res <- ggmend(data_ggmend$egenes, data_ggmend$cissnps,
#'               effectth = 0, corth = 0, standardize = TRUE,
#'               maxperms = 1000, threads = 1,seed = 1)
#' netdata <- convertnet(res)
#' plotnet(netdata$nodes, netdata$edges)

ggmend <- function(ys, xs, effectth = 0, corth = 0, standardize = TRUE, maxperms = 1000, threads = 1, seed = 0){

  if(!require(pbapply)){
    warning("instll pbapply r package")
    install.packages("pbapply")
  }

  stopifnot(nrow(ys) == nrow(xs))
  stopifnot(ncol(ys) == ncol(xs))

  message(paste0("n = ", ncol(xs), " / m = ",nrow(xs)))
  message(paste0("Data filtering: |effectsize| > ",effectth," and ", "|correlation| > ",corth))

  if(standardize){
    message("standardizing")
    egenes <- t(scale(t(ys)))
    snps <- t(scale(t(xs)))
  }else{
    egenes <- ys
    snps <- xs
  }


  c <- cor(t(egenes))
  c2 <- matrix(0,nrow(c),ncol(c))
  c2[which(abs(c) > corth, arr.ind=T)] <- 1
  c2[lower.tri(c2,diag=TRUE)] <- 0

  BsgeGenes <- sapply(1:nrow(c), function(i){geteffect(egenes[i,], snps[i,])})

  c3 <- matrix(0,nrow(c),ncol(c))
  c3[which(abs(BsgeGenes) > effectth),]  <- 1
  c3[,which(abs(BsgeGenes) > effectth)] <- c3[,which(abs(BsgeGenes) > effectth)] + 1
  c3[lower.tri(c3,diag=TRUE)] <- 0
  indexesToRun <- data.frame(which(c2 == 1 & c3 == 2, arr.ind = T))

  message(paste0("Number of snp/gene pair: ",nrow(indexesToRun)))


  raw_res <- pbapply::pblapply(1:nrow(indexesToRun), function(indtestnum){
    ind1 <- indexesToRun$row[indtestnum]
    ind2 <- indexesToRun$col[indtestnum]

    S1 <- as.matrix(snps[ind1,])
    S2 <- as.matrix(snps[ind2,])
    G1 <- as.matrix(egenes[ind1,])
    G2 <- as.matrix(egenes[ind2,])

    Discoverbgg(S1, S2, G1, G2, maxperms, threads, seed)
  })

  res <- do.call(rbind, raw_res)
  res$FDR <- p.adjust(res$pval, method = "BH")

  save_res <- data.frame(g1 = rownames(egenes[indexesToRun$row,]),
                         g2 = rownames(egenes[indexesToRun$col,]),
                         indexesToRun,
                         res)
  return(save_res)
}


#' Export network to Cytoscape
#'
#'This function exports a network in edge and node list files in a format suitable for importing to \href{https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html}{visNetwork}
#'
#'
#' @param res output ggmend
#' @param fdrth Significant level of FDR. default = 0.05
#' @examples
#' data(data_ggmend)
#' res <- ggmend(data_ggmend$egenes, data_ggmend$cissnps,
#'               effectth = 0, corth = 0, standardize = TRUE,
#'               maxperms = 1000, threads = 1,seed = 1)
#' netdata <- convertnet(res)
#' plotnet(netdata$nodes, netdata$edges)
convertnet <- function(res, fdrth = 0.05){
  sig_res <- res[res$FDR < fdrth,]
  message("filtering < " ,fdrth,"...")
  sig_res$idx <- apply(sig_res[,c("Bg1g2","Bg2g1")],1,which.max)
  edges <- rbind(sig_res[sig_res$idx == 1,c("g1","g2","Bgg")],
                 sig_res[sig_res$idx == 2,c("g2","g1","Bgg")])

  colnames(edges) <- c("from","to","Bgg")

  edges[,c(1,2)] <- apply(edges[,c(1,2)],2,as.character)

  nodes <- data.frame(id = sort(unique(c(edges$from,edges$to))),
                      label = sort(unique(c(edges$from,edges$to))))
  message("nodes: ",nrow(nodes)," / edges: ",nrow(edges))
  return(list(edges = edges, nodes=nodes))
}


#' Network visualization for GGMend
#'Network visualization using \href{https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html}{visNetwork}
#' @param nodes See online documentation \href{https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html}{visNetwork}
#' @param edges See online documentation \href{https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html}{visNetwork}
#' @examples
#' data(data_ggmend)
#' res <- ggmend(data_ggmend$egenes, data_ggmend$cissnps,
#'               effectth = 0, corth = 0, standardize = TRUE,
#'               maxperms = 1000, threads = 1,seed = 1)
#' netdata <- convertnet(res)
#' plotnet(netdata$nodes, netdata$edges)
plotnet <- function(nodes, edges){

  if(!require(visNetwork)){
    warning("instll visNetwork r package")
    install.packages("visNetwork")
  }

  require(visNetwork, quietly = TRUE)


  visNetwork(nodes, edges) %>%
    visNodes(size = 10,color = "gold", shape = "circle") %>%
    visEdges(color = "grey", arrows ="to") %>%
    visOptions(manipulation = TRUE)%>%
    visInteraction(navigationButtons = TRUE) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
}
