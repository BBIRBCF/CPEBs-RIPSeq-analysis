## Index gene sets in the expression matrix used
indexFunction <- function(y, gsetsel ="Hallmarks",minsize=10, maxsize=300,
                          gspath = "/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/",
                          geneRandom = FALSE, forPRgsea = FALSE){
 gsets <- paste0(gspath, list.files(gspath))
 gsets <- gsets[regexpr("\\.gmt", gsets) > 0]
 gsets <- gsets[regexpr("v3",gsets)>0]
 gsets <- gsets[regexpr(gsetsel,gsets)>0]

 auxg <- strsplit(readLines(gsets),"\t")
 GSET <- sapply(auxg,function(x) x[3:length(x)])
 names(GSET) <- sapply(auxg,function(x) x[1])
# GSET <- lapply(GSET,tolower)
 if(forPRgsea)  index <- GSET
 else{
#     rownames(y) <- tolower(rownames(y))
     index <- sapply(GSET, function(x) which(rownames(y)%in%x))
     if(geneRandom)
         index <- sapply(index, function(x) sample(seq_len(dim(y)[1]), length(x)))

     len <- sapply(index,length)
     index <- index[len > minsize & len < maxsize]
 }
 return(index)
}

