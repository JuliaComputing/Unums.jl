module UnumX

import Base: print, show, showcompact, convert, zero, one
import Base: +, -, *, /, abs, ==, <, <=

export Unum, Unum22, Ubound22, Bnum, Bbound, isexact, esmax, fsmax, emax, fmax


include("interval.jl")
include("bnum.jl")
include("unum.jl")


end # module
