using UnumX

import UnumX: esmax, fsmax

function all_exact_pos(U)
    E = Array(U,0)

    for es in 1:esmax(U)
        for fs in 1:fsmax(U)
            for e = 0:(1<<es-1)
                for f = 0:(1<<fs-1)
                    push!(E,U(false,e,f,false,es,fs))
                end
            end
        end
    end
    E
end

E = all_exact_pos(Unum22)

F = [convert(Float64,x) for x in E]
Fu = unique(F)
sort!(Fu)

R = [convert(Rational{Int},x) for x in E]
Ru = unique(R)
sort!(Ru)
