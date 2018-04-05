# WordGraph2Index

Generate a probabilistic word index from a given word lattice in *Standard Lattice Format*:
...
S. Young, J. Odell, D. Ollason, V. Valtchev, and P. Woodland. The HTK
Book: Hidden Markov Models Toolkit V2.1. Cambridge Research Laboratory
Ltd, Mar. 1997.
...


## Compilation

WordGraph2Index is implemented in C++ and STL, and depends on the
package [zlib >=1.2.11](http://zlib.net/).



make distclean; make CFLAGS=-m32  # compile for 32 bits
make distclean; make CFLAGS=-m64  # compile for 64 bits


## Usage

wordGraphTools -i d06-030-02-01.wg -d 0.01 -z w > ll

## Acknowledgments

If you find useful this toolkit in your research, please cite:

```
@article{wgidx16,
 title = "{HMM} word graph based keyword spotting in handwritten
                  document images ",
 journal = "Information Sciences ",
 volume = "370-371",
 pages = "497 - 518",
 year = "2016",
 note = "",
 issn = "0020-0255",
 doi = "http://dx.doi.org/10.1016/j.ins.2016.07.063",
 url = "http://www.sciencedirect.com/science/article/pii/S0020025516305461",
 author = "Alejandro H. Toselli and Enrique Vidal and Ver{\'o}nica
                  Romero and Volkmar Frinken"
}
```
