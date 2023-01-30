using Graphs
using SparseArrays

include("ModularEvolutionScript.jl")

"""Takes an interaction matrix, converts it to a list of edges
Loses activation/inhibition information"""
function interaction2edgelist(tNet::AbstractMatrix)
    xi, yi, vi = findnz(sparse(tNet))    
    return [(xi[i], yi[i]) for i ∈ eachindex(xi)], vi
end

function nfeedbackloops(J::AbstractMatrix)
    E, F = interaction2edgelist(J)
    E = Edge.(E)
    G = SimpleDiGraph(E)
    
    C = simplecycles_limited_length(G, 10)
    pfbl = 0
    nfbl = 0
    for loop in C
        v = 1
        for i in loop
            v *= F[i]
        end
        if v == 1
            pfbl += 1
        else 
            nfbl += 1
        end
    end
    return pfbl, nfbl
end

npositivefeedbackloops(J::AbstractMatrix) = nfeedbackloops(J)[1]
nnegativefeedbackloops(J::AbstractMatrix) = nfeedbackloops(J)[2]

