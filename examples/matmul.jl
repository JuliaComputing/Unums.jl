using UnumX

import UnumX: Ubound34

X = randn(4,4)
U = map(Ubound34,X)

X*X
U*U
