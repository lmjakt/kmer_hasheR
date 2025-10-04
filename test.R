require("Biostrings")

seq <- readDNAStringSet("test.fa")
sseq <- as.character(seq[1])
seq.2 <- readDNAStringSet("fLopPis1.1.hap2.fa")
seq.3 <- as.character(seq.2[["SUPER_6"]])
rm(seq.2)

source("kmer_hash.R")
ptr.1 <- make.kmer.hash(seq, 10, do.sort=FALSE)
pos.1 <- kmer.pos(ptr.1, opt.flag=15)

ptr.2 <- make.kmer.hash(sseq, 17)
pos.2 <- kmer.pos(ptr.2, opt.flag=7)

## this does take rather a lot of time. That's not that surprising
## it has to be said...
source("kmer_hash.R")
system.time(
    ptr.3 <- make.kmer.hash(seq.3, 32)
)
##    user  system elapsed
##  9.683   1.652  11.335 
## 10.421   2.248  12.670
## 10.598   2.269  12.867 
##  7.721   1.212   8.934
## previously:
##  user  system elapsed 
## 9.574   0.856  10.429 

## kmer sequences only:
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=1)
)
## I seem to get very variable numbers for this
##   user  system elapsed 
## 81.621   2.075  83.699
## 77.609   1.776  79.387
## first time slow, second time faster?
## 12.401   0.136  12.537 

ptr.4 <- make.kmer.hash(seq.3, 32)

system.time(
    pos.4 <- kmer.pos(ptr.4, opt.flag=1)
)
## first time:
##    user  system elapsed 
##  11.240   0.296  11.536  



## pos only
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=2)
)
##  user  system elapsed 
## 3.697   0.496   4.193 
##  3.427   0.348   3.775 

## counts only
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=8)
)
##  user  system elapsed 
## 0.969   0.136   1.105 
## 0.761   0.000   0.762 

## kmer, pos, count:
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=1 + 2 + 8)
)
##   user  system elapsed 
## 51.900   0.808  52.710 
## 17.209   0.560  17.771 

## this doesn't work due to ofverflow in R;
## as it attempts to set up a table 9.7e9 columns
## and 3 rows. That is simply too much to deal with.
## we would need to use some form of minimiser to
## deal with this; But for now we can simply subset the sequence to something more reasonable.
system.time(
    pos.3 <- kmer.pos(ptr.3, opt.flag=15-4)
)
##   user  system elapsed
##   user  system elapsed
## 78.505   2.759  81.265 
## 51.900   0.808  52.710 
## why so much more than before? 
## 15.922   0.036  15.959


## consider to add a range option to the kmer.pos function...
## but for now lets just try:
ptr.4 <- make.kmer.hash(subseq(seq.2[["SUPER_6"]], 3e6, 5e6), 32)

## this may still fail due to overflow... 
pos.4 <- kmer.pos(ptr.4, opt.flag=6)

## but this is not going to make an easy drawing...
nrow(pos.4$pair.pos) / 1e6
## [1] 625.307

plot.new()
plot.window(xlim=range(pos.4$pos[,'pos']), ylim=range(pos.4$pos[,'pos']))
with(pos.4, segments( pair.pos[,'x'], pair.pos[,'y'], pair.pos[,'x'] + 32, pair.pos[,'y'] + 32))

plot.pairs <- function(pos, xlim, ylim, k=32){
    b1 <- with(pos, pair.pos[,'x'] >= xlim[1] & pair.pos[,'x'] <= xlim[2])
    b2 <- with(pos, pair.pos[,'y'] >= ylim[1] & pair.pos[,'y'] <= ylim[2])
    pairs <- as.data.frame(pos$pair.pos[ b1 & b2, ])
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    axis(1)
    axis(2)
    with(pairs, segments(x, y, x+k, y+k))
    invisible(pairs)
}

xlim <- c(3e5, 4e5)
plot.pairs(pos.4, xlim=xlim, ylim=xlim, k=32)

xlim <- c(3.27e5, 3.4e5)
tmp <- plot.pairs(pos.4, xlim=xlim, ylim=xlim, k=32)



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
