require("Biostrings")

seq <- readDNAStringSet("test.fa")
seq.2 <- readDNAStringSet("fLopPis1.1.hap2.fa")

source("kmer_hash.R")
ptr.1 <- make.kmer.hash(seq, 10)
pos.1 <- kmer.pos(ptr.1, opt.flag=7)

ptr.2 <- make.kmer.hash(seq, 32)
pos.2 <- kmer.pos(ptr.2, opt.flag=7)

## this does take rather a lot of time. That's not that surprising
## it has to be said...
system.time(
    ptr.3 <- make.kmer.hash(seq.2[["SUPER_6"]], 32)
)
##  user  system elapsed 
## 7.721   1.212   8.934
## previously:
##  user  system elapsed 
## 9.574   0.856  10.429 

## this doesn't work due to ofverflow in R;
## as it attempts to set up a table 9.7e9 columns
## and 3 rows. That is simply too much to deal with.
## we would need to use some form of minimiser to
## deal with this; But for now we can simply subset the sequence to something more reasonable.
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=7)
)
##   user  system elapsed 
## 15.922   0.036  15.959

## consider to add a range option to the kmer.pos function...
## but for now lets just try:
pos.4 <- make.kmer.hash(subseq(seq.2[["SUPER_6"]], 3e6, 9e6), 32)

## I can generate x and y positions for a dotplot by:
pts <- with(pos.3, tapply( pos[,'pos'], pos[,'i'], function(x){
    l <- length(x)
    if(l == 1)
        return(matrix(nrow=0, ncol=2))
    m <- matrix(nrow=l * (l-1)/2, ncol=2)
    r <- 0
    for(i in 2:l-1){
        for(j in (i+1):l)
            m[(r <- r + 1), ] <- x[c(i,j)]
    }
}))

pos$pos[,2] <- pos$pos[,2] - 10
pos.1 <- with(pos, pos[ pos[,'i'] == 1, 2])
sseq.1 <- as.character( subseq(seq[rep(1,length(pos.1))], pos.1, pos.1 + 9))
## sseq.1 are all the same. But
start.1 <- start( matchPattern(sseq.1[1], seq[[1]]) )
## but I have 15 starts?
all( as.character( subseq( seq[rep(1, length(start.1))], start.1, start.1 + 9)) == sseq.1[1] ) ## TRUE?

sum( start.1 %in% pos$pos[,2])
which(pos$pos[,2] %in% start.1 )

seq.ch <- as.character(seq[[1]])

i <- with(pos, which(pos[,2] %in% start.1))
pos$kmer[ pos$pos[i,'i'] ]
##  [1] "GTTAAAAAAA" "GTTAAAAAAA" "GTTAAAAAAA" "GTTAAAAAAA" "GTTAAAAAAA"
##  [6] "GTTAAAAAAA" "GTTAAAAAAA" "GTTAAAAAAA" "GTTAAAAAAA" "TAGTAACAAA"
## [11] "TTGTAGAAAA" "TTGTAGAAAA" "TTGTAGAAAA" "TTGTAGAAAA" "TTGTAGAAAA"

with(pos, as.character(subseq(seq[rep(1, length(i))], pos[i,2], pos[i,2]+9)))
##      SUPER_1      SUPER_1      SUPER_1      SUPER_1      SUPER_1      SUPER_1 
## "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" 
##      SUPER_1      SUPER_1      SUPER_1      SUPER_1      SUPER_1      SUPER_1 
## "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" 
##      SUPER_1      SUPER_1      SUPER_1 
## "AGAAGCCCCC" "AGAAGCCCCC" "AGAAGCCCCC" 
## it's a bit strange; I'm not getting all the sequences. I probably need to check my
## the offsets calculated. But these are unfortunately 64 bit.. 
