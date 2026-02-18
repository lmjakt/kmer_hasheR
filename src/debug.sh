#!/bin/bash

## all optimisation turned off: -O2 --> ... 

gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic -Og  -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_hash.c -o kmer_hash.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic  -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_pos.c -o kmer_pos.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_reader.c -o kmer_reader.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_tree.c -o kmer_tree.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c suffix_hash.c -o suffix_hash.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c thread_queue.c -o thread_queue.o
gcc -I"/usr/lib64/R/include" -DNDEBUG   -I/usr/local/include    -fpic   -Og -Wall -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -fstack-protector-strong -funwind-tables -fasynchronous-unwind-tables -fstack-clash-protection -Werror=return-type -flto=auto -g  -c kmer_util.c -o kmer_util.o
gcc -shared -L/usr/lib64/R/lib -flto=auto -o kmer_hash.so kmer_hash.o kmer_pos.o kmer_reader.o kmer_tree.o suffix_hash.o thread_queue.o kmer_util.o -L/usr/lib64/R/lib -lR
