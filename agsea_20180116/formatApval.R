

formatApval <- function(x){
    st = "ns"
    st <- ifelse(x < 0.25, "+", st)
    st <- ifelse(x < 0.10, "*", st)
    st <- ifelse(x < 0.05, "**", st)
    st <- ifelse(x < 0.01, "***", st)
    if(any(is.na(st))) st[is.na(st)] <- "nc"
    return(st)
}
