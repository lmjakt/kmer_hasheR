# kmer_hash

An extension to `R` that uses `kvec` and `khash` from `klib`
[https://github.com/attractivechaos/klib](https://github.com/attractivechaos/klib)
to create a k-mer index of a given sequence.

This index is stored as an external pointer from which k-mer positions
and off diagonal coordinates can be obtained. This makes it possible
to create classical dot plots over large regions.

The functions are currently limited to values of k up to 32 due to the
hashing of 64 bit integers representing two bit encoded sequences.
Regions containing N's are excluded; regions containing ambiguity codes
will be incorrectly handled the hashing function only checks for
"N" and assumes any non-N characters are from the set of A, C, T, G
(lower or upper case).

# Usage and functions

To use the functions you need to be able to compile the source code. 
I have only tested this on Linux, but I expect that it will also work
in MacOS. You should be able to compile it in Windows, but I'm not sure
whether the tools needed are freely available or not.

1. Clone the repository to a suitable location. You should obtain
   `README.md` (this file), `kmer_hash.R`, `test.R`, `test.fa`,
   `src/kmer_hash.c` in addition to several header files from
   `klib`. The `test.R` file contains code I used to test the
   functions; it may provide some example usage, but it is *not* up to
   date and includes expressions written for earlier versions of the
   functions and which may not work now. The code also includes
   references to files not present in the repository.
2. Compile the source (`src/kmer_hash.c`) by invoking `R CMD SHLIB kmer_hash.c` in the
   source directory. This should create a shared object (`kmer_hash.so`) which
   needs to loaded using `dyn.load()` (see next point).
3. From an active `R` session, source the `kmer_hash.R` file
   (`source("<dir to file>/kmer_hash.R")`). This will call `dyn.load()` on the shared
   object and create the functions described below.
   
Note that the functions should not be considered as stable; they may change in the
future without this document being updated.

## `make.kmer.hash(seq, k, do.sort=FALSE)`

Returns an external pointer containing the hash and some other information
(eg. the value of `k`).

### Arguments

1. `seq`: a single sequence; this may be in any format that can be converted to
   `R`'s standard character vector format. The hash will be constructed from
   the first element of the resulting character vector. The sequence should contain
   only valid nucleotides or N's but may be upper or lower case. Non-N ambiguity
   codes are not valid. The sequence must be longer than or equal to `k`.
2. `k`: a numeric value giving the length of the k-mer; `k` must be within the range
   of 1 to 32 (inclusive).
3. `do.sort`: if `TRUE` then positions for individual `k-mers` will be sorted from
   small to large. This should in fact be pointless as the current method for defining
   the positions should guarantee that they are already sorted. I implemented the option
   without considering this; I'm keeping the code as it is as it might be useful in the
   future.

   
## `kmer.pos(ex.ptr, opt.flag)`

Returns a named `R` `list` object containing the following elements:

1. `kmer`: a character vector giving the sequences of all k-mers obtained.
2. `pos`: an integer matrix with two columns ('i', 'pos') giving the kmer
   index (`i`, position in the vector k-kmer) and positions in the sequence (`pos`).
   The number of rows of this matrix should be $1 + l - k$, where $l$ is the
   length of the sequence and $k$ is `k` (the length of the k-mer).
3. `pair.pos`: off diagonal pairs of positions in a matrix with three
   columns (`i`, `x`, `y`) that indicate the kmer index (`i`), and two
   positions where the kmer was found where $x < y$ (i.e. it gives the
   upper triangle of off diagonal alignments).
4. `count`: the number of occurences of the k-mer in the sequence. One value for
   each kmer as in `kmer`.
   
   Pair positions are only given for k-mers found more than once in the sequence.
   However, common k-mers will give rise to large numbers of pairs ($n * (n-1)/2$)
   and this can exceed the size of matrices that `R` can handle. For example
   using a 32-mer index of a 40 million base pair sequence from *L. piscatorius*
   resulted more than $9\times10^9$ pairs[^mem]. Attempting to obtain such tables
   will result in an error, but does not seem to crash the `R` session. It seems
   likely that such errors will also result in a memory leak as the error
   will be raised when attempting to allocate memory `R` for the `R` matrix
   and will not provide any opportunity for freeing memory allocated for the
   positions.
      
Each of the elements may be `NULL` depending on the value of `opt.flag`.

[^mem]: This is from memory; the actual numbers may have been larger by some factor
	of 10.

### Arguments

1. `ex.ptr`: an external pointer to a k-mer hash as created by `make.kmer.hash()`.
2. `opt.flag`: a bitwise flag indicating which elements of the list should be returned
	interpreted as follows:
	- 1: Return a list of kmer sequences.
	- 2: Return k-mer positions.
	- 4: Return pair positions.
	- 8: Return counts.

# Performance, limitations and future

## Performance

Creating a 32-mer hash for ca. 40 million base pairs took
approximately 10 seconds on a single core of an "Intel(R) Xeon(R) Gold
6248R CPU @ 3.00GHz"; obtaining coordinates for these (excluding pair
positions) took around 15 seconds. It appears that the first time that
kmer-sequences are retrieved is much more time consuming (taking up to
80 seconds). This is likely because `R` stores strings as a hash;
the first time it comes across a string it will need to create a new
entry in the hash. In general you should not need to obtain the k-mer
sequences if positions are retrieved as you can then get the k-mers
from the original sequence.

I have not checked memory requirements; but these will be in excess of
$4 * l$ ($l$ is the length of the sequence) for the hash itself. Obtaining
the positions (`kmer.pos()`) will triple this requirement as this will
create a matrix containing the `kmer` index and the position for every 
position in the sequence. In addition, there will be a temporary requirement
that exceeds this as `kmer.pos()` collects data using `kvec` data structures
that are copied to `R` data structures afterwards. This might seem wasteful
but to avoid this requires creating a custom memory allocator for `R`.
It may be possible to do this using `allocVector3`[^mem2]; I believe this
allows you to transfer already allocated memory blocks to `R`'s management
(i.e. become garbage collected). However, I've not found much documentation
about this function.

[^mem2]: The name of the allocator function is similar to this but may be
	different.


## Limitations and future

### Limitations

- `valgrind` did not reveal any memory leaks, but I expect that memory
  leaks *will* result from attempting to create excessively large pair
  pos tables (or any other operation leading to excessive sizes).
- The `kmer.pos()` function should be modified to check the size of data before using
  `R`'s memory allocation functions.
- `k` is limited to 32; this is somewhat small, but sufficient for my current
  needs. This limitation arises from encoding `k-mers` using 64 bit unsigned
  integers; these can be obtained and `hashed` much faster than `const char`
  objects. I suspect that it should be possible to make use of arrays of 
  integers, but I haven't checked how to do this with `khash`.

### Future

I have written this in order to look at some repetitive
sequences. Whether or not I extend the functions depends on how useful
they turn out to be compared to other tools. I'm fairly sure that
similar functionality exists elsewhere, and maybe even within existing
`R` packages, but for this use-case it was easier to re-implment than
to go looking for these. 

If I have more complex requirements that I cannot find implemented
elsewhere then I may provide additional functionality. There are
a number of obvious improvements which can be made, but I will only
implement these if I have a need. If you have some such need, then,
feel free to discuss.
