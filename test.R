source("kmer_hash.R")

## no Biostrings package. and on the plane!
read.fasta <- function(fn){
    lines <- readLines(fn)
    id.i <- grep("^>", lines)
    beg.i <- id.i + 1
    end.i <- c(id.i - 1, length(lines))[-1]
    seq <- mapply(function(b, e){
        paste(lines[b:e], collapse="")
    }, beg.i, end.i)
    names(seq) <- sub("^>([ ]+).*$", "\\1", lines[id.i])
    seq
}

## A -> T 0 -> 2
## C -> G 1 -> 3
## G -> C 3 -> 1
## T -> A 2 -> 0
## (v + 2) % 4

rev.comp <- function(seq){
    nucs <- c("A", "C", "T", "G")
    seq.2b <- lapply(seq, function(x){ bitwAnd(3, bitwShiftR(utf8ToInt(x), 1)) })
    seq.2b <- lapply(seq.2b, function(x){ rev((x + 2) %% 4) })
    sapply(seq.2b, function(x){ paste(nucs[1+x], collapse="") })
}

seq <- read.fasta("test.fa")
seq.rc <- rev.comp(seq)

nchar(seq) ## 59940

## The following sequence is not provided by the repository
## consider replacing with a different one.
require("Biostrings")
seq.2 <- readDNAStringSet("fLopPis1.1.hap2.fa")
seq.3 <- as.character(seq.2[["SUPER_6"]])
rm(seq.2)
nchar(seq.3) / 1e6 ## 40.00614

ptr.1 <- make.kmer.hash(seq, 10, do.sort=FALSE)
ptr.2 <- make.kmer.hash(seq.rc, 10, do.sort=FALSE)

## This crashes; don't do unless you are fixing the issue
k.pairs <- kmer.pairs(ptr.1, ptr.2)

pos.1 <- kmer.pos(ptr.1, opt.flag=15)

ptr.2 <- make.kmer.hash(sseq, 17)
pos.2 <- kmer.pos(ptr.2, opt.flag=7)

## count only
source("kmer_hash.R")
counts.ptr <- count.kmers(seq, c(21, 0, 2), NULL)
counts.ptr <- count.kmers(seq, c(21, 1, 2), counts.ptr)

counts.k <- kmer.pos(counts.ptr, opt.flag=1 + 2 + 8)
counts <-  do.call(rbind, with(counts.k, tapply(pos[,'pos'], pos[,'i'], eval)))
counts <- data.frame(kmer=counts.k$kmer, n=counts[,1])
o <- order(counts$n, decreasing=TRUE)
head( counts[o,], n=30 ) ## the top ones are telomore sequences. That's nice.
tail( counts[o,], n=30 )
plot( counts$n[o], type='l' )

system.time(
    counts.ptr.1 <- count.kmers.fq("test.fastq.gz", c(21, 0, 2), NULL)
)
##  user  system elapsed 
## 0.049   0.000   0.050 

counts.k <- kmer.pos(counts.ptr.1, opt.flag=1 + 2 + 8)
counts <-  do.call(rbind, with(counts.k, tapply(pos[,'pos'], pos[,'i'], eval)))
counts <- data.frame(kmer=counts.k$kmer, n=counts[,1])
o <- order(counts$n, decreasing=TRUE)
head(counts[o,], n=30)

system.time(
    counts.ptr.1 <- count.kmers.fq("test.fastq.gz", c(21, 1, 2), counts.ptr.1)
)
##  user  system elapsed 
## 0.038   0.001   0.038 

counts.k <- kmer.pos(counts.ptr.1, opt.flag=1 + 2 + 8)


## we can use both fasta files and fq files..
## repeat.fq contains reads that are repeats of:
## 10 repeats of 5 repeats of ACTGG
## so 50 repeats of ACTGG,
## 49 repeats of CTGGA
## the other kmers.
source("kmer_hash.R");
ptr.0 <- count.kmers.fq("repeat.fq", c(5, 1, 0, 100, 30, 100))

ptr.0 <- count.kmers.fq("repeat.fq", c(5, 1, 0, 100, 30, 100))
spc.0 <- kmer.spec.kt(ptr.0, 10000) 
## these are actually the correct numbers:
## one kmer is counted 5000 times (ACTGG)
## 4 kmers are counted 4900 times (the alternate frames)
## and there are no other counts.
table(spc.0)
## spc.0
##    0    1    4 1019 
## 9998    1    1    1 

ptr.1 <- count.kmers.fq("repeat.fq", c(5, 1, 1, 100, 30, 100))
spc.1 <- kmer.spec.kt(ptr.0, 10000)

## try the combination of array and hash.. 
source("kmer_hash.R");
ptr.0 <- count.kmers.fq.sh("repeat.fq", c(20, 1, 20, 100, 30, 100))

spc.1 <- kmer.spec.sh(ptr.0, 10000)

### and try with a large file to see how the memory grows...
### first with the assembly
ptr.0 <- count.kmers.fq.sh("fLopPis1.1.hap1.fa", c(20, 1, 20, 300, 30, -1), NULL)

system.time(
    ptr.0 <- count.kmers.fq.sh("fLopPis1.1.hap1.fa", c(20, 1, 20, 300, 30, -1), ptr.0)
)
##    user  system elapsed 
## 273.290   5.904 279.206 

## and that uses about 10GB of memory.
## for a total of 6.58e08 k-mers

## this now uses a bit more memory, around 11.5 GB, and it seemed slower?
system.time(
    ptr.1 <- count.kmers.fq.sh("fLopPis1.1.hap1.fa", c(20, 1, 24, 300, 30, -1), NULL)
)
##    user  system elapsed 
## 334.044   8.078 342.145 

system.time(
    ptr.0 <- count.kmers.fq.sh("fLopPis1.1.hap1.fa", c(20, 1, 20, 300, 30, -1), NULL)
)
##    user  system elapsed 
## 353.992   9.592 363.600 
## Difference is not big enough to be sure if it is meaningful. Memory use would also seem
## to be similar.

time.2 <- system.time(
    ptr.3 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(20, 1e5, 20, 300, 30, -1), NULL)
)

## this seems to settle in at 4 seconds / 100000 reads even though the number of new k-mers
## does decrease substantially. The memory use seems to be in the region of 10GB, so, not
## so bad. For four seconds, the 4.4e8 reads would take 4 to 5 hours.

## the khashes here should be 16 times smaller; one would think that would improve the performance somewhat..
## but didn't seem to make much of a difference when I tried for the genome itself
## and most hashes will not have that many entries anyway.. 
time.3 <- system.time(
    ptr.3 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(20, 1e5, 24, 300, 30, 50e6), NULL)
)
##     user   system  elapsed 
## 1938.152   43.177 1981.429 

time.4 <- system.time(
    ptr.4 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(20, 1e5, 26, 300, 30, 50e6), NULL)
)
## strangely the process seems to go from fast to slower? I don't quite understand that.
##     user   system  elapsed 
## 2084.360   49.133 2133.628 

time.5 <- system.time(
    ptr.5 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(20, 1e5, 18, 300, 30, 50e6), NULL)
)
##     user   system  elapsed 
## 1378.953   40.206 1419.227 
##
## This is unexpected; it seems that the hash is really vey efficient and that I might be better off
## having more bits in the suffix.
## but my experience of just using the hash was not very good, so whether this actually generalises I don't
## know..

## run overnight.. 
time.6 <- system.time(
    ptr.6 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(20, 1e5, 18, 300, 30, -1), NULL)
)
##     user   system  elapsed 
## 13184.63   482.86 13668.27 
## 3.7 hours, for a total of 434,000,000
## reads: 31,752.37 / second
## if we can use 10 threads, maybe 22 minutes.
## memory use about 11 GB.
## total number of k-mers: 7.73e+08
## From a total of 37,793,452,976
## ~38e9 words; average coverage 48
## most frequent word: CACACACACACACACACACA
## counts:             15,576,025 times

source("kmer_hash.R");

time.7 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.7 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 18, 30, 30, 1e6, 2^14), NULL)
)

## Unfortunately, I'm not getting much more than 3 to 400% CPU usage.
## and the time taken gets worse with more threads?

## increase the total number of hashes; thus dividing up the prefixes more evenly?
time.8 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.8 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 22, 30, 30, 1e6, 2^14), NULL)
)


time.9 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.8 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 16, 30, 30, 1e6, 2^14), NULL)
)


## count.kmers.fq.sh.mt is broken with a deadlock at the moment. Needs fixing..

time.10 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.8 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 16, 5, 30, 1e6, 2^14), NULL)
)

time.12 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.8 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 16, 1, 30, 1e6, 2^14), NULL)
)

## lets have a massive number of threads
time.13 <- system.time(
    ## params are: k, report_n, prefix_bits, thread_n, min_quality, max_read_n, queue_buffer_size
    ptr.8 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 16, 1, 200, 1e6, 2^14), NULL)
)


## naive, no threading at all; much faster? why?
time.11 <- system.time(
    ptr.9 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 18, 300, 30, 1e6), NULL)
)

source("kmer_hash.R"); ptr.9 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 18, 300, 30, 1e6), NULL)

source("kmer_hash.R"); ptr.7 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 18, 30, 2, 1e5, 100), NULL)

source("kmer_hash.R"); ptr.7 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 18, 30, 3, 1e6, 100), NULL)

ptr.9 <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 18, 300, 30, 1e5), NULL)

rp.time <- sapply(1:10, function(t){
    system.time(count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 18, 30, t, 1e6, 100), NULL))
})

## this time include a larger number of reads.. 
rp.time.2 <- sapply(2^(0:6), function(t){
    system.time(count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 18, 30, t, 1e7, 100), NULL))
})


plot(2^(0:6), rp.time.2['elapsed',])
plot(2^(0:6), 1e7/rp.time.2['elapsed',])
points(1:10, 1e6/rp.time["elapsed",], col='red')

## That unfortunately suggests that things get worse with larger numbers of reads.
## lets consider the effect of increasing the prefix size for a given number of threads (32)

p.bits <- c(18, 20, 22, 24)
t <- 32
rp.time.3 <- sapply(p.bits, function(p){
    system.time(count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, p, 30, t, 1e7, 100), NULL))
})
## There's no real pattern here:
## user.self  819.722 844.490 885.525 948.337
## sys.self    36.245  38.732  40.064  34.051
## elapsed     44.043  41.178  43.753  44.821
## user.child   0.000   0.000   0.000   0.000
## sys.child    0.000   0.000   0.000   0.000


ptr.1 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 30, t, 1e7, 100), NULL)
## this takes more time than I would like. It should be possible to multithread this as well.
spc <- kmer.spec.sh(ptr.1, max.count=10000)

ptr.st <- count.kmers.fq.sh("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e6, 22, 300, 30, 1e7), NULL)
spc.st <- kmer.spec.sh(ptr.st, max.count=10000)

plot(spc.st, type='l')
plot(spc.st[2:10], type='l')
plot(spc[2:10], type='l')

plot(spc, spc.st)

plot(log10(spc[2:200]), type='l')

## how does the number of different k-mers affect the situation?
## change the min quality to 20
t <- 32
system.time(
    ptr.1 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, t, 1e7, 100), NULL)
)
##     user   system  elapsed 
## 1070.797   38.641   54.131
## not super bad.

## but quitting takes a long time to recover all the broken memory
source("kmer_hash.R")

t <- 32
system.time(
    ptr.q20 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, t, -1, 100), NULL)
)
##      user    system   elapsed 
## 40115.858  1230.218  1897.518 
1897.518 / 60 ## 31.6253

t <- 33
system.time(
    ptr.q20 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, t, -1, 100), NULL)
)
##      user    system   elapsed 
## 42062.939  1288.534  1356.804 
1897.518 / 1356.804 ## 1.39
## 32 threads is almost 40% slower..
## 22.6 minutes instead of 32 ;-)

## extend a hash..
system.time(
    ptr.1 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, 33, 1e7, 100), NULL)
)
##     user   system  elapsed 
## 1051.646   45.360   35.275

spc.1 <- kmer.spec.sh(ptr.1, 1e4)

## if this works,, then it should exactly double
## the spc.1 counts
## 
system.time(
    ptr.2 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, 33, 1e7, 100), ptr.1)
)
##    user  system elapsed 
## 979.580  31.207  33.984 
1e7 / 33.984 ## 294256.1

spc.2 <- kmer.spec.sh(ptr.2, 1e4)
spc.3 <- kmer.spec.sh(ptr.1, 1e4)
## These now look complete..
## the problem was that we have to always use kh_exist when iterating
## over entries (I guess we can just get the bucket otherwise).

## anyway, try
t <- 33
system.time(
    ptr.q20 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 22, 20, t, -1, 100), NULL)
)
##      user    system   elapsed 
## 39740.561  1181.841  1259.443 
## 
system.time(
    ptr.q20 <- count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R2.fastq", c(21, 22, 20, t, -1, 100), ptr.q20)
)
##      user    system   elapsed 
## 35299.155  1185.282  1141.009 

(1259.443 + 1141.009) / 60 ## 40 minutes

system.time(
    spc.q20 <- kmer.spec.sh(ptr.q20, 1e4)
)
##   user  system elapsed 
## 10.794   1.992  12.786 

get.peaks <- function(x){
    k <- 3:length(x)
    j <- k-1
    i <- j-1
    1 + which(x[j] > x[i] & x[j] > x[k])
}

sum(spc.q20) / 1e6 ## 1065.572
sum(spc.q20[ -(1:10) ]) / 1e6 ## 693.9576

plot(spc.q20, type='l')

plot(spc.q20[2:200], type='l', ylim=c(0, 2e7))
spc.q20.p <- get.peaks(spc.q20)


k <- 3:length(spc.q20)
j <- k-1
i <- j-1
peaks.i <- which( spc.q20[j] > spc.q20[i] & spc.q20[j] > spc.q20[k] )

plot(5:200, spc.q20[5:200], type='l')
abline(v=peaks.i+1)
abline(v=peaks.i[1]/2)

## it's quite possible that a prime number (or at least even)
## number of threads will give a better performance..
rp.time.2 <- sapply(30:33, function(t){
    system.time(count.kmers.fq.sh.rp("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 18, 30, t, 1e7, 100), NULL))
})
## user.self  849.476 861.308 858.477 904.650
## sys.self    35.974  37.074  27.031  37.150
## elapsed     35.153  31.987  48.260  30.666

1e7 / rp.time.2['elapsed',]
## [1] 284470.7 312627.0 207210.9 326094.0
## 31 threads better than 32; 33 way better
326094.0 / 207210.9 ## 1.57373

## and this while we have quite some activity from aruna as well... 

## one consumer to keep things simple.. 
## source("kmer_hash.R"); ptr.7 <- count.kmers.fq.sh.mt("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(21, 1e5, 18, 1, 30, 100, 100), NULL)

time.7 <- system.time(
    spc <- kmer.spec.sh(ptr.6, max.count=10000)
)

plot(spc[2:200], type='l')

system.time(
    counts.ptr.ass <- count.kmers.fq("fLopPis1.1.hap1.fa", c(17, 1, 0, 300, 30, -1), NULL)
)

system.time(
    counts.ptr.ass <- count.kmers.fq("fLopPis1.1.hap1.fa", c(17, 1, 0, 300, 30, -1), NULL)
)



time.1b <- system.time(
    counts.ptr.ass <- count.kmers.fq("fLopPis1.1.hap2.fa", c(17, 1, 0, 300, 30, -1), counts.ptr.ass)
)

##
## Count k-mers from all reads from a 147 Gb fastq (not compressed)
## this will most likely take several hours. It's not very well optimised
time.2 <- system.time(
    counts.ptr.1 <- count.kmers.fq("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(17, 1e6, 0, 300, 30, -1), NULL)
)
##
##     user   system  elapsed 
## 3774.407  387.426 4162.010 
##
time.3 <- system.time(
    counts.ptr.1 <- count.kmers.fq("breiflabb_hiseq_qtrimmed_nophix_R1.fastq", c(17, 1e6, 0, 300, 30, -1), counts.ptr.1)
)

system.time(
    counts.g <- kmer.spec.kt(counts.ptr.ass, 10000)
)
##   user  system elapsed 
## 38.688   1.332  40.022 

plot(log10(counts.g[2:100]), type='l')

system.time(
    counts.r <- kmer.spec.kt(counts.ptr.1, 10000)
)
##   user  system elapsed 
## 38.414   1.295  39.710 

plot(counts.r[4:400], type='l')



## this for:
## 434000000 reads -> 7.12e+08 kmers (1.18e+05 new): 9.391e+06 clicks 9 sec / 1000000 reads:  9.391e+00 clicks / read (IO: 4.390e+05)
## total allocated: 1 (7.12e+08)  new: 0
## prefix / suffix bits: 0 / 36
## Estimated memory use: 2.75e+11
## Total number of words counted 39122080676 max count: 19206319
## Most common k-mer: CACACACACACACACACA

## This is way slower than it should be ?
## it obtains counts from 2^36 positions
## that is about 68e9 different positions. Even
## at 1e9 per second it's going to take more than a minute.
system.time(
    counts.n <- kmer.spec.kt(counts.ptr.2, 1000)
)
##    user  system elapsed 
## 157.796   7.776 165.577 

plot(log10(counts.n))

plot(counts.n[3:100])



### This seems to go into some form of infinite loop. I'm not sure what happens, but it got
### to use about 160GB of data and 100% of CPU for several days. Somewhere the process is hanging.
system.time(
    ptr.3 <- make.kmer.hash(seq.3, 32)
)
##  user  system elapsed 
## 7.980   1.703   9.684 

## That is for 40 Mbp. For the complete Lophius genome
## we might then expect something like 20 * 10 -> 200
## seconds. That's not too bad.

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
