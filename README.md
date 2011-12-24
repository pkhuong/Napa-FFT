Napa-FFT
========

Introduction
------------

The Napa-FFT library implements a cache-oblivious Fast Fourier
Transform (and its inverse) for complex valued signals, in mostly
portable Common Lisp.

The specialized routine generator and the transform functions were
written by Paul Khuong.  Andy Hefner provided the signal windowing
code; he also wrote the documentation for Bordeaux-FFT, from which
this document is heavily inspired.

The interface is not identical to that of Bordeaux-FFT, but is meant
to be as easy to use in the common case, while still offering some
lower-level, lower-overhead, facilities.

The library's name reflects its American origins, and was suggested by
Zach Beane.

Installation
------------

Napa-FFT can be loaded by using ASDF; the system is defined in
`napa-fft.asd`.  The system is typically registered with ADSF by
adding a symbolic link to that file in a directory under
`asdf:*central-registry*`.  If Quicklisp is installed, [FILL ME].

Once registered with ASDF, Napa-FFT can be loaded by executing
`(asdf:oos 'asdf:load-op "napa-fft")`, or, if Quicklisp is installed,
`(ql:quickload "napa-fft")`.

Napa-FFT generates a fairly large function at compile-time; it may
take a long while to compile `interface.lisp` (the `(fill-vectors 8)`
form), or even exhaust the heap on small machines.  Lowering the value
passed to `fill-vector` to 4 should yield a sufficiently simpler
function.

Normal FFT functions
--------------------

Napa-FFT only supports a single type of input and output vectors:
arrays of double complexes values.  Scratch space and configuration
data can be cached through instances of `fft-instance`.

_It is important to note that Napa-FFT currently only supports input
sizes that are powers of two; I have vague plans to improve on the
matter, but I wouldn't bet on it happening anytime soon._

> Type `complex-sample`: `(complex double-float)`

> Type `complex-sample-array`: `(simple-array complex-sample 1)`

> Type `fft-instance`: `(defstruct fft-instance ...)`

Most Napa-FFT users will only use a single function, `FFT` or `SFFT`
(simple FFT), along with `make-fft-instance`.  

> Function `sfft (src &optional direction dst)`

Transforms the input `src`.  `src` is coerced to a
`complex-sample-array` as needed, and `direction` defaults to
`forward`.  When `dst` is provided, the results are written in that
`complex-sample-array`.

> Function `fft (instance src &optional dst)`

Applies the the transformation specified in `instance` to the
`complex-sample-array` `src`, returning the result as a
`complex-sample-array` of the same length as `src`.  When `dst` is
non-`nil`, the result is written in it.

> Function `make-fft-instance (size direction &optional scale)`

Returns a new `fft-instance` to perform a transformation on
`complex-sample-array`s of length `size`, `:forward` or `:backward`
depending on `direction`.  If `scale` is not provided (or `nil`),
forward transformations do not scale the results, while backward ones
normalize by `(/ 1d0 size)`; otherwise any scaling factor can be
provided (as a double float value)

> Function `reset-fft-instance (instance &key size direction scale)`

Mutates the `fft-instance` to specify the `size`, `direction` or
`scale` of the transformation.  When `scale` is provided but `nil`, it
is reset according to `direction`; this is usually the desired
behavior when changing `direction`.

> Variable `*fft-instance*`

This variable is bound to the `fft-instance` that `sfft` reuses across
invocations.  In a threaded setting, it should be bound to `nil` when
threads are spawned; `sfft` will initialize it on demand.

Windowing
---------

See the Bordeaux-FFT manual for now.

Low-level interface
-------------------

> Function `%fft (size dst dst-offset src src-offset direction tmp1
> tmp2 tmp-offset)` 

This function applies a transformation on a vector `size` samples in
`src`, starting from `src-offset`, and writes the result in `dst`,
starting from `dst-offset`.  The direction is 1 for forward
transformations, and -1 for backward.  `tmp1` and `tmp2` are two
`complex-sample-array` used as scratch space of length `size`,
starting from `tmp-offset`.

Note that no scaling is performed, even for backward transformations.

> Function `%fft-scale (size dst dst-offset src src-offset direction
>  scale tmp1 tmp2 tmp-offset) `

This function is the same as the previous one, except `scale` is a
double float value by which the result is multiplied.

This is useful to normalize the result of backward transformations.

> Function `find-fft-function (size direction &key scale)`

Returns two values: an `FFT-function` and a
`complex-sample-array` of twiddle factors.  `size` is the size of the
transformation to perform and `direction` 1 or -1 for a forward or
backward transform.

When `scale` is false, the function must be called with the following
arguments:

1. twiddle factor array
2. size
3. destination vector
4. destination offset
5. source vector
6. source vector offset
7. temporary vector 1
8. temporary vector 2
9. temporary vectors offset

When `scale` is true, the function has an additional argument:

1. twiddle factor array
2. size
3. destination vector
4. destination offset
5. source vector
6. source vector offset
7. scale factor (`double-float`)
8. temporary vector 1
9. temporary vector 2
10. temporary vectors offset

Performance
-----------

Napa-FFT mostly implements the 6-step Fast Fourier Transform
algorithm.  It is a variant of Cooley-Tukey's algorithm, and, more
specifially, of Gentleman and Sande's 4-step algorithm [1] that is
better suited to computers with caches (or with slow tape storage).  I
find the algorithm is best described in an article by David Bailey [2]
on FFTs in external memory.

The idea is to decompose the transformation as two sets of smaller
transforms.  If we consider the input as a matrix (stored in row-major
order), we can first transform each column separately, transpose the
result, apply an element-wise multiplication with a matrix of twiddle
factors, and transform each column again.

The access patterns of that algorithm are fairly bad on all but the
smallest inputs.

Rather than transforming by columns, the 6-step version first
transposes the input.  In a row-major layout, values in the same row
are adjacent in memory, and each sub-transformation is thus performed
on an small range of addresses.  Of course, the transposition must be
reversed at the end to obtain a correct result vector.

This algorithm, when transpositions are implemented by recursing on
submatrices, like the transformation itself, ensures asymptotically
cache-optimal access patterns.

The constant factors associated with this deeply recursive structure
can be daunting.  Napa-FFT also includes a simple code generator for
fixed-size transpositions and transformations.  In the default
configuration, transformations and transpositions up to size 256 are
compiled ahead of time (and exploited by operations on larger
inputs).

There are hand-written generators for transformations of size 2 to 16.
The two larger ones (8 and 16) simply implement the four-step
algorithm and aren't particularly well-tuned.

On my 2.8 GHz X5660, Napa-FFT is faster than Bordeaux-FFT on all input
sizes.  For at most 256 elements, the speed-up factor varies between 2
and 4.  The speed-up then slowly decreases, hitting bottom at
approximately 1.2 for 2048 elements.  For larger inputs (up to 64M
elements), the speed-up factor varies between 1.3 and 1.6.  Scaled
(e.g. backward) transformations should be also slightly faster with
Napa-FFT, since it performs the scaling as part of the last set of
recursive fourier transformations.

On small or medium-sized inputs, it is usually around 3 times as slow
as FFTW.  Larger inputs (e.g. 1 million elements) favour Napa-FFT's
layout slightly, with a slowdown on the of 2, even hitting 1.6 on
humongous inputs (64M elements).

Several things could be improved, performance-wise.  In particular:

 * the hand-written FFTs aren't particularly good;
 * it only exploits SIMD operations indirectly, via CL's native
   complex arithmetic;
 * windowing could be applied as part of the initial transpose;
 * scaling could be applied in the final transpose, which is
   bandwidth-bound and arithmetic-free.

Overall, however, obtaining performance so close to the state of the
art in pretty much straight CL and relatively little code is
surprising, even for a bandwidth-bound task.

Paul Khuong

References:

[1] Gentleman W. M., and G. Sande, "Fast Fourier transforms—for fun
and profit," Proc. AFIPS 29, 563–578 (1966),
http://www.computer.org/portal/web/csdl/doi/10.1109/AFIPS.1966.83

[2] Bailey, David H., "FFTs in external or hierarchical memory,"
J. Supercomputing 4 (1), 23–35 (1990), http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/fftq.pdf

[3] Bordeaux-FFT, http://vintage-digital.com/hefner/software/bordeaux-fft/manual.html

[4] Fastest Fourier Transform in the West (FFTW), http://www.fftw.org/
