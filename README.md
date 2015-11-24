# BH_diagonalize

A Julia exact diagonalization code for the [Bose-Hubbard model](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model) with a focus on entanglement entropy.

The basis is enumerated using the method of [Szabados et al., 2011](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf).
See also [Zhang et al., 2011](http://arxiv.org/pdf/1102.4006v1.pdf).


## Requirements

* [ArgParse](https://github.com/carlobaldassi/ArgParse.jl) (`Pkg.add("ArgParse")`)
* [JeszenszkiBasis](https://github.com/0/JeszenszkiBasis.jl) (`Pkg.clone("https://github.com/0/JeszenszkiBasis.jl.git")`)


## Examples

* `julia BH_main.jl --help`
* `julia BH_main.jl --out output.dat --ee-all 2 4 4`
* `julia BH_main.jl --out output.dat --site-max 4 --ee-all 3 6 6`
