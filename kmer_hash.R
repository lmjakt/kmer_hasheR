dyn.load( paste(dirname(sys.frame(1)$ofile), "src/kmer_hash.so", sep="/") )

## the kmer positions should be sorted naturally as they are added
## by traversal of the sequence. 
make.kmer.hash <- function(seq, k, do.sort=FALSE){
    .Call("make_kmer_h_index", as.character(seq),
          as.integer(k), as.integer(do.sort))
}

## hash.ptr: a suitable pointer to a hash object..
## params: k, source, source.no
## seq: nucleotide sequences
## returns the external pointer.
count.kmers <- function(seq, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers", hash.ptr, params, seq);
}


## params: k, report_n, prefix_bits, max memory, min quality, max_read_n
count.kmers.fq <- function(fq.file, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers_fastq", hash.ptr, params, fq.file);
}

## params: k, report_n, prefix_bits, max memory, min quality, max_read_n
count.kmers.fq.sh <- function(fq.file, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers_fastq_sh", hash.ptr, params, fq.file);
}

## params: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
count.kmers.fq.sh.mt <- function(fq.file, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers_fastq_sh_mt", hash.ptr, params, fq.file);
}

count.kmers.fq.sh.rp <- function(fq.file, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers_fastq_sh_rp", hash.ptr, params, fq.file);
}


kmer.spec.kt <- function(ptr, max.count){
    .Call("kmer_spectrum_ktree", ptr, as.integer(max.count))
}

kmer.spec.sh <- function(ptr, max.count){
    .Call("kmer_spectrum_suffix_hash", ptr, as.integer(max.count))
}


kmer.pos <- function(ex.ptr, opt.flag){
    tmp <- .Call("kmer_positions", ex.ptr, as.integer(opt.flag))
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

kmer.pairs <- function(ptr.a, ptr.b){
    tmp <- t(.Call("kmer_pair_pos", ptr.a, ptr.b))
    colnames(tmp) <- c("a", "b")
    tmp
}
