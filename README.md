# LLBounds

A project for computing lower bounds for the Lichnerowicz Laplacian on homogeneous spaces G/H with G simple and Einstein standard metric.

Running the code requires an installation of both Sage and LiE.

## Usage

Open sage in a terminal and type:

    load("Sym2Bounds.sage")

Then enjoy the implemented functionalities, or navigate to the subfolders you are interested in and load the .sage files there.

The Sage function `WeylCharacterRing` is the standard to handle compact semisimple Lie groups. Typing

    G = WeylCharacterRing("A2xB3",style="coroots")
    
initializes `G` as the compact Lie group A2xB3. The optional parameter refers to the basis in which highest weights are given. For example `G(1,0,1,0,0)` will correspond to the tensor product of standard representations of A2 and B3.

Most of the implemented functions require a *branching rule* as input. This can be obtained using the Sage function `branching_rule`.

LiE instead uses a *restriction matrix* for branching. Since LiE works with row vector rather than column vectors, this matrix is the transpose of the one defined in the article.

## Documentation

### Sym2LowerBounds

`Sym2LowerBounds(G,H,b,startscan=0,endscan=-1)` takes as input the WeylCharacterRings `G` and `H` corresponding to a **simple** Lie group *G* and a **semisimple** subgroup *H*, a branching rule `b` (alternatively the restriction matrix as used in LiE, given as a list of lists), and the two optional arguments `startscan` and `endscan` which refer to the range of Casimir eigenvalues of *G* for which Fourier modes should be scanned. If left blank, the appropriate upper bounds is automatically determined using the *crude estimate*.

The output consists of fibrewise estimates for *A\*A* and *q(R)* on the standard homogeneous space *G/H* as well as Fourier-mode-wise lower bounds for the Lichnerowicz Laplacian (if *q(R)>E* does not follow from the first step). Additionally, the program tries to find Killing tensors in each Fourier mode by comparing *Hom(.,Sym^2_0)* with *Hom(.,Sym^2)*.

This is printed to the console and also to a `.txt` file, named by the Cartan types of *G* and *H*.

### Sym2LowerBoundsWithTorus

`Sym2LowerBoundsWithTorus(G,Hss,rm,startscan=0,endscan=-1)` takes as input the WeylCharacterRings `G` and `Hss` corresponding to a **simple** Lie group *G* and the **semisimple** part of a **non-semisimple** subgroup *H*, plus a LiE restriction matrix `rm` (given as a list of lists) for the embedding of *H* in *G*. The optional arguments and the output are as above.

### Sym2LowerBoundsFullFlag

`Sym2LowerBoundsFullFlag(cartantype,startscan=0,endscan=-1)` takes as input just the Cartan type string `cartantype` of a **simple** Lie group *G*. The optional arguments and the output are as above, except that the fibrewise estimates are not printed for each *H*-isotype to avoid cluttering.

## Subfolders

There is a subfolder containing a `.sage` file for each family and for the two classes of exceptions. Make sure to comment/uncomment or alter the `for` loops at the end of each file that starts the computation as needed.

## Troubleshooting

Here are some of the most common errors that may occur.

### Object table overflow

This can be remedied by raising the LiE parameter `maxobjects` (standard 99999). Type, for example:

    lie.eval("maxobjects 999999")
    
Sometimes a higher value is needed. However, this may lead to longer computation times.


### Hash table overflow

There seems to be no way to manipulate the hash table used in the LiE interface. We tried to be parsimonious with the amount of LiE objects in use, but it *will* exceed the system capabilities at some point. We recommend dividing up the Fourier modes by use of the `startscan` and `endscan` parameters.

### Non-virtual decomposition failed

This seems to be a problem within LiE itself. A reproducible example is the branching of the *D124*-module \[2,2,0,0,0,...\] to *E8*.
