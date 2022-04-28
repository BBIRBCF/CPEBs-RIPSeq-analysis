
weighted.canberra <- function(rank1, rank2, K){
        m1 <- pmin(rank1,K+1)
        m2 <- pmin(rank2,K+1)

        return(sum(abs(m1 - m2)/(m1 + m2)))
}

weighted.canberra.mat <- function(rankmat, K){
    ms <- apply(rankmat, 2,function(rank1) pmin(rank1,K+1))
    COR <- sapply(1:dim(ms)[2], function(k) sapply(1:dim(ms)[2], function(j)  return(sum(abs(ms[,k] - ms[,j])/(ms[,k] + ms[,j])))))
    rownames(COR) <- colnames(COR) <-  colnames(rankmat)
    return(COR)
}
