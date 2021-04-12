# Single Quartet Check

A small program that determines the correctness of quartet trees in regard to some reference tree.

This repository is a fork of [lutteropp/quartet_check](https://github.com/lutteropp/quartet_check).

## How does it work?

Let's say we want to check if the tree `((A,B),(C,D));` is correct. In this case, the program calculates the lowest common ancestor of `(A,B)`, `(A,C)`, `(A,D)`, `(B,C)`, `(B,D)` and `(C,D)` in the reference tree. It is then assumed that the quartet tree is correct if and only if the distance (number of edges between the nodes) between `lcu(A,B)` and `lcu(C,D)` is larger than the distance between `lcu(A,C)` and `lcu(B,D)` as well as the distance between `lcu(A,D)` and `lcu(B,C)`.

## Compilation

You need gcc 7 to compile the code. (To be honest, I'm not quite sure about that. However, gcc 7 worked for me, while gcc 10 did not.)

```
mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make
```
If everything worked correctly, you should find the binary under `build/single_quartet_check`.

## Usage

```
./single_quartet_check <FASTA file> <quartet tree file> <reference tree file>
```

- FASTA file: Input sequences in FASTA format
- quartet tree file: Quartet trees in Newick notation, one tree per line
- reference tree file: Reference tree in Newick notation

### Output

For each quartet tree in the input file, the program outputs either `1` (if the tree is correct) or `0` (if it is incorrect) to stdout, one result per line.

## License

Single Quartet Check is provided under GPL-3.0 license. It is based on [Quartet Check](https://github.com/lutteropp/quartet_check), which was originally created by [Sarah Lutteropp](https://github.com/lutteropp).

This program depends on [Genesis](https://github.com/lczech/genesis), which is also licensed under GPL-3.0.
