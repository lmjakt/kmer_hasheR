dyn.load( paste(dirname(sys.frame(1)$ofile), "src/kmer_hash.so", sep="/") )

## the kmer positions should be sorted naturally as they are added
## by traversal of the sequence. 
make.kmer.hash <- function(seq, k, do.sort=FALSE){
    .Call("make_kmer_h_index", as.character(seq),
          as.integer(k), as.integer(do.sort))
}

kmer.pos <- function(ex.ptr, opt.flag){
    tmp <- .Call("kmer_positions", ex.ptr, as.integer(opt.flag))
    names(tmp) <- c('kmer', 'pos', 'pair.pos')
    if(!is.null(tmp$pos)){
        tmp$pos <- t(tmp$pos)
        colnames(tmp$pos) <- c("i", "pos")
    }
    if(!is.null(tmp$pair.pos)){
        tmp$pair.pos <- t(tmp$pair.pos)
        colnames(tmp$pair.pos) <- c("i", "x", "y")
    }
    tmp
}

