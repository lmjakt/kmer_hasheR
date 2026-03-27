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
## count.kmers.fq.sh.mt <- function(fq.file, params, hash.ptr=NULL){
##     params <- as.integer(params)
##     .Call("count_kmers_fastq_sh_mt", hash.ptr, params, fq.file);
## }

## params:
## k, prefix_bits, min_q, thread_n, max_reads, max_mem, source_n, source
## NOTE: the prefix_bits does not seem to matter very much to timing, but
##       larger prefix_bits seem to take more memory.
##       Note that it seems to me that a khash_t can only have up to
##       2^32 buckets; it's thus possible that suffix_bits > 32
##       could cause crashes. This may also explain why an initial experiment
##       using a straight khash_t( uint64_t ) failed to complete (seem to hang).
count.kmers.fq.sh.rp <- function(fq.file, params, hash.ptr=NULL){
    params <- as.integer(params)
    .Call("count_kmers_fastq_sh_rp", hash.ptr, params, fq.file);
}

seq.kmer.depth.sh <- function(hash.ptr, seq, k){
    .Call("seq_kmer_depth_sh", hash.ptr, as.character(seq), 
          as.integer(k))
}

kmer.spec.kt <- function(ptr, max.count){
    .Call("kmer_spectrum_ktree", ptr, as.integer(max.count))
}

kmer.spec.sh <- function(ptr, max.count){
    .Call("kmer_spectrum_suffix_hash", ptr, as.integer(max.count))
}

kmer.spec.sh.n <- function(ptr, max.count, comb, comb.inner, source.min){
    .Call("kmer_spectrum_suffix_hash_n", ptr, as.integer(max.count),
          as.integer(comb), as.integer(comb.inner), as.integer(source.min))
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
