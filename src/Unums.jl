module Unums

using Compat

import Base: print, show, showcompact, convert, zero, one
import Base: +, -, *, /, ^, abs, ==, <, <=, sqrt, signbit, in, isfinite, isnan

export Unum, Unum22, Ubound22, Unum34, Ubound34, Bnum, Bbound, isexact, esmax, fsmax, emax, fmax, limits, print_bits


include("bnum.jl")
include("unum.jl")


end # module
