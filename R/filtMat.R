#'Filtering species
#'@description Filtering low abundance/occurence species
#'
#'@param data : data is a matrix. Row is sample, column is species.
#'@param th : th is an abundance by given an integer. 
#'@param nsample ; nsample is total sample.
#'@param minth : max.itr is minimum rate of occurence in samples. 
#'
#' @export

filtMat <- function(x=data, th=10, nsample=nrow(x),
                 minth=0.05, maxth=0.95){
    freq <- colSums(x>th)/nsample
    cat( sprintf('ASVs frequency range : minimum is %s, maximum is %s\n', 
                 round(min(freq), 2), round(max(freq), 2)))
    cat(sprintf('Discard ASVs in which appeared less than %s samples\n', round(nsample*minth)))
    y <- x[, freq>minth & freq<maxth]
    z <- y[rowSums(y)>0, ]
    cat( sprintf('Remained sample/asvs is %s/%s\n', nrow(z), ncol(z)))
    return(z) }