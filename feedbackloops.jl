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
    
    C = simplecycles_limited_length(G, 10, 10^8)
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
ntotalfeedbackloops(J::AbstractMatrix) = sum(nfeedbackloops(J))

# function that maps the difference between 
# number of positive and negative feedback loops 
# between 0 and 1 
function ratiomap(p, n)
    t = p+n
    1 - ((p-n)/(t))^2 
end


function PfbAndNfbScore(J::AbstractMatrix) 
    n, p = nfeedbackloops(J)
    t = n + p
    r = ratiomap(n, p)
    return t * r
end

"""=
Takes in a network, and returns score for selecting for
positive and negative feedback.
Rationale: 
Q := ((p-n)/(p+n))^2 transforms nPfb / nNfb to the interval [0, 1] 
1 (if the ratio is unbalanced) and 0 (if the ratio is 1)

1 - Q transforms Q to
1 (if Q is 0) and 0 (if Q is 1)

(1-Q)*(p+n) gives a higher score for: 
- increasing the number of feedback loops (p+n)
- getting the ratio nPfb / nNfb closer to 1 
 
But (1-Q) * (p+n) simplifies to 
4*n*p/(p+n)
"""
PfbAndNfbFitness(p, n) = 4*n*p/(p+n)
PfbAndNfbFitness(J::AbstractMatrix) = PfbAndNfbFitness(nfeedbackloops(J)...)
