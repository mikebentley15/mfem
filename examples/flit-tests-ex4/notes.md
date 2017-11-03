# Investigation of Example 04

I have created a file called `narrow_down.py` that semi-automates the process
of running all of these combinations.  Here I document the experiments I have
performed.

## Problem Description

For example 04, all Intel compilations disagreed with the basline compilation
output, namely `gcc -O0`.  The reproducibility problems from the Intel
compilations resulted in simply the Intel linker.  Using the GCC linker does
not result in reproducibility problems.

One possability for this difference is that the Intel linker could be linking
against different implementations of standard functionality, such as the math
library.  This was prompty checked, and it was found that the shared libraries
used are all exactly the same, including the same file path on the file system.
There was one exception, and that is `libdl.so` that Intel includes and is not
included from `gcc`.  The `libdl.so` library provides functionality to
dynamically load symbols from shared libraries at runtime, bypassing the
dynamic linker.  However, it has been confirmed that none of the symbols of
`libdl.so` are used in the test executable.  That was a dead end.

Example 04 was not the only test that exibited this kind of behavior.  However,
Example 04 is the only one to have been confirmed to be caused by the linker.
I would not be surprised to find out the other examples are the same problem,
and perhaps even optimization of the same functionality within MFEM.  All tests
affected this way are:

1. Example 04
2. Example 05
3. Example 09
4. Example 10

In past manual experimentation, I showed that the file called
`../../mesh/nurbs.cpp` was the only one necessary to keep in the executable in
order to get the same bad result.  This was actually narrowed down to a single
function within the `nurbs.cpp` file.

The assumption was that the problem was limited to this single file.  To test
that assumption, the selection was inverted to have only `nurbs.cpp` compiled
and linked into the shared library and everything else compiled and linked into
the executable.  According to the assumption, this should have produced the
good baseline result.  However, it still produced the exact same bad result at
before.  So the problem is not limited to this single identified function or
identified file, but the problem is manifested with this single file alone.

From this point, a new experiment needs to be performed in order to obtain a
new hypothesis of how the Intel linker does these unsafe optimizations.  That
is the scope of this work here.

## Methodology

In this set of experiments, I take one set of files and compile them into a
shared library using the `gcc` linker.  I then take the remainder of the files
and link them together to create an executable with the Intel linker.

If some of the files are stored in a shared library linked with a different
linker than the Intel linker, then when linked together with the test
executable, the symbols in the shared library are off-limits from the
optimizations.

One thing to note is that if certain symbols exist both in the executable and
in the shared library, the default behavior of the dynamic linker is to use the
symbols from the executable instead of those from the shared library.  This
functionality can be overridden, but it is a bit complicated to do.

There is a Makefile template called `Makefile.in`.  This template has places
for the set of files to be linked into the shared library with the baseline
linker and the set of files to be linked into the executable with the
problematic linker.  The python script `narrow_down.py` takes this template and
creates a new Makefile by inserting these two sets as well as the filename of
the Makefile (in order to ensure the executable gets compiled and not skipped).
Unique filenaming techniques are used so that this procedure can occur in
parallel.

So far, there have been three approaches to determine the effects

1. Compile only one file into the shared library and all others into the
   executable.  Called `only_one_gt()`.
2. Compile only one file into the executable and all others into the shared
   library.  Called `only_one_dev()`.
3. Compile two files into the executable and all others into the shared
   library.  Called `all_two_combos_dev()`.
4. Compile all combinations of "interesting files" into the shared library, and
   all others into the executable.  Called `all_interesting_files()`.  The
   interesting files are those found from the previous experiments.

## Results

Here I describe the results of each of these approaches:

### Results for `only_one_gt()`

The hope was that removing only one source file from the executable would cause
it to not generate the bad result.  No such luck.  All of these produced the
bad result, so the problem is not so simple that it can be isolated to a single
object file.

### Results for `only_one_dev()`

This is a different approach.  The hope was that I could identify a single file
that when linked by itself, the problem would still be there.  All other files
are compiled into the shared library, thus limiting the linker's capabilities
only to this one object file.  The results are a bit surprising.  The following
files still demonstrated the problem when compiled alone and linked against the
large shared library.

1. `tests/Example04.cpp`
2. `../../mesh/nurbs.cpp` (this is not surprising since it was found by binary
   search manually)
3. `../../fem/intrules.cpp`

This was interesting since all three files behaved the same in terms of the
linker.  You need to only choose _any one of these files_ to link into the
executable to get the problem.  I still do not know what this means.

This is also surprising that only a single file is necessary to cause the
problem, since when any one of these three files are the only ones left out of
the linking into the executable, the problem still persists.  Surely the
problem is more complicated than that.

### Results for `all_two_combos_dev()`

Here, I compile all combinations of two files into the test executable while
keeping all others in the shared library.  The thought here is that there may
be some files not identified in `only_one_dev()` that could be shown to be
problematic or interesting.

Since there are 49 files, the amount of combinations are
  $\frac{49 \cdot 48}{2} = 1176$,
which is equivalent to `choose(49, 2)` in combinatorics.

Unsurprisingly, all files identified from `only_one_dev()` triggered the
problem to persist, regardless of the second file included as well.  It is nice
that at least that was consistent!

But are there some interesting pairings that cuase problems that do not include
the three identified files from `only_one_dev()`?  The answer is actually yes!
The following pairings of files in the executable caused the problem to
persist.

1. `../../fem/fe.cpp` with `../../linalg/densemat.cpp`
2. `../../fem/fe.cpp` with `../../mesh/mesh.cpp`

Why does `../../fem/fe.cpp` not trigger the problem by itself, but does when
paired with either `../../linalg/densemat.cpp` or `../../mesh/mesh.cpp`?  Why
does the pairing of `../../linalg/densemat.cpp` and `../../mesh/mesh.cpp` not
cause the problem to show itself?

The plot thickens...

Note, I do not search over all combinations of three.  With 49 files, this
amounts to `choose(49, 3) = 18424`.  Doing all combinations of two took a long
time, but doing combinations of three would take more than 15 times as long.
Maybe I'll try it later, but not now.  Granted, I could reduce the search space
by removing the already found interesting files.  Hopefully, I can narrow down
on the problem without needing to do this search.

### Results fro `all_interesting_files()`

From the previous experiments, the following files have been identified as
interesting:

1. `tests/Example04.cpp`
2. `../../mesh/nurbs.cpp`
3. `../../fem/intrules.cpp`
4. `../../fem/fe.cpp`
5. `../../linalg/densemat.cpp`
6. `../../mesh/mesh.cpp`

Keep in mind that we have a total of 49 source files involved here.

The goal with this experiment is to see what combinations of the above
"interesting files" are required to be put into the shared library to make the
problem go away.  So, all combinations of the above files are tried while
putting all other files into the executable.  The following combinations were
found to make the problem disappear:

Combo #1:

1. `../../fem/fe.cpp`
1. `../../fem/intrules.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

Combo #2:

1. `../../fem/fe.cpp`
1. `../../fem/intrules.cpp`
1. `../../linalg/densemat.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

Combo #3:

1. `../../fem/fe.cpp`
1. `../../fem/intrules.cpp`
1. `../../mesh/mesh.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

Combo #4:

1. `../../fem/intrules.cpp`
1. `../../linalg/densemat.cpp`
1. `../../mesh/mesh.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

Combo #5:

1. `../../fem/fe.cpp`
1. `../../fem/intrules.cpp`
1. `../../linalg/densemat.cpp`
1. `../../mesh/mesh.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

Only these 5 combinations out of the $2^6 = 64$ possible combinations caused
the problem to disappear.

There are a few observations that can be made from these combinations.  The
first observation is that Combo #1 trumps Combos #2, #3, and #5, meaning that
all of the files in Combo #1 are found in #2, #3, and #5.  That means that
really the only distinct combos presented here are Combos #1 and #4.

Now, these two distinct combos have the following files in common:

1. `../../fem/intrules.cpp`
1. `../../mesh/nurbs.cpp`
1. `tests/Example04.cpp`

The files that are different are

- Combo #1: `../../fem/fe.cpp`
- Combo #4: `../../linalg/densemat.cpp` and `../../mesh/mesh.cpp`

If you think about it, given the previous results, this is entirely expected.
We found that if one of the files found from `only_one_dev()` are in the
executable, the problem persists.  This means all three need to be in the
shared library to prevent problems.  Furthermore, if one of the two pairings of
files found from `all_two_combos_dev()` are in the executable, then the problem
persists as well.  This means that either `../../fem/fe.cpp` needs to not be in
the shared library or both `../../linalg/densemat.cpp` and
`../../mesh/mesh.cpp` need to be in the shared library.  This is exactly what
this experiment shows.  In fact, it proves that there are no more files of
interest, and therefore there is no need to do an experiment of all
combinations of three files.

## Different versions of Intel compiler

Moving to a different version of the intel compiler (namely from 17.0 to
16.0.3), the problem seems to have different characteristics.  I don't quite
know how to interpret this result.

Now trying on the CAB cluster at Livermore since the intel license ran out on
Fractus.  The versions of compilers in this setting (which could change if I
looked into how to use different compilers there simultaneously) are:

- Intel Compiler: 16.0.1 20151021
- GCC: 4.9.2 20150212
- Clang: 4.0.1

One of the difficulties about the module system on the Livermore computers is
that when you load a newer version of GCC, the newer version of Clang gets
unloaded and visa-versa.  Same iwth loading a different version of the intel
compiler, it makes the GCC and Clang compilers go back to the system installed
defaults, which is not desirable.

These versions are almost the same to the poster, except for Clang is at 4.0.1
instead of 4.0.0.  I hope that does not make much of a difference.  Through
running the experiments again, it looks like it doesn't.  The results from
Fractus match exactly those from CAB, which is good news.

## Difference in Symbols

I took one of the compiled problematic ones from the `only_one_gt()` experiment.  Specifically, I chose to have all objects into the shared library and only `../../fem/intrules.cpp` compiled into the executable.  This I did with linking using Intel and another time with GCC.  I then wanted to see what the difference was in symbols between the two.

This experiment is brought on because when looking that the verbose output for
the link command generated by the intel compiler, I noticed a number of math
libraries from Intel that were statically linked against the executable.
Specifically, libraries that I never mentioned in my original command to the
intel compiler.

There were 12 different types of symbols found in each resultant executable:

1. **B**: In the uninitialized data section
1. **D**: In the initialized data section
1. **R**: In the read-only data section
1. **T**: In the text (code) section
1. **U**: _Undefined symbol_
1. **V**: Weak object
1. **W**: Weak symbol -- not a weak object
1. **b**: In the uninitialized data section
1. **d**: In the initialized data section
1. **r**: In the read-only data section
1. **t**: In the text (code) section
1. **w**: Weak symbol -- not a weak object, default value initialized

One interesting thing to do is to see what symbols are in the gcc compiled version and not in the intel compiled version.  Those symbols are these:

1. `cos@@GLIBC_2.2.5`
1. `pow@@GLIBC_2.2.5`
1. `sin@@GLIBC_2.2.5`

All three of these symbols are of type **U** which means they are undefined and
would need to be resolved by the dynamic linker at runtime.  The fact that
these symbols (under these exact names) are not in the intel compiled binary
hints that maybe different versions of these functions are being used.  In
fact, the associated symbols in the Intel compiled binary are:

1. `cos`
1. `pow`
1. `sin`

respectively, all of type **T** which means they are fully defined and in the
text (code) section of the binary.  These three symbols were statically linked
in instead of dynamically linked to the GLIBC implementation.

It looks like these functions are coming from the following static library on CAB:

`/opt/intel-16.0/linux/compiler/lib/intel64_lin/libimf.a`

Knowing that, I can probably test out the two different versions of each of
these three function.  I could also try to do this same analysis on
compilations with the other files of interest to see if there are any more
surprises.

## Things To Do And Try

1. Demonstrate this problem happens even when the FLiT testing framework is not
   used.
    - Verified for example 4
    - **Could not verify for example 5!!!**  Need to revisit
    - **Could not verify for example 9!!!**  Need to revisit
    - Verified for example 10
2. Show that the other three tests are also victims of the Intel linker
   optimizations.
    1. Example 05
    2. Example 09
    3. Example 10
3. Invoke the linker directly and make sure I can cause this to happen
4. Use the linker directly in these experiments
5. Run the linker with the `-v` option to show verbose information.  It should
   describe exactly what the linker is doing, and hopefully will give the
   answer.
6. Diff the dissassembly of one of the minimal bad executables with the same
   minimal good executable.
7. Run the experiment again, and make sure those that do not demonstrate the
   same problem are exactly identical to the answer from the baseline
   compilation.  These experiments were only checking bitwise equivalence with
   the problem answer and never checking against the baseline answer.

