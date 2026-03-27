#!/bin/bash

## stop on first error.
set -e

## all optimisation turned off: -O2 --> ... 

## Note that you need to set the value of RINC to the location where
## your appropriate R header files can be found. This can be found
## by looking at what is printed when you run the usual R CMD SHLIB
## command.

##RINC="/usr/lib64/R/include"
RINC="/home/lmj/R/R-4.4.1/include"



gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic -Og  -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_hash.c -o kmer_hash.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic  -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_pos.c -o kmer_pos.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_reader.c -o kmer_reader.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_tree.c -o kmer_tree.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c suffix_hash.c -o suffix_hash.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c thread_queue.c -o thread_queue.o
gcc -I${RINC} -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_util.c -o kmer_util.o
gcc -shared -L/usr/lib64/R/lib -flto=auto -o kmer_hash.so kmer_hash.o kmer_pos.o kmer_reader.o kmer_tree.o suffix_hash.o thread_queue.o kmer_util.o -L/usr/lib64/R/lib -lR
