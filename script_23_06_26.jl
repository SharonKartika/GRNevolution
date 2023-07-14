#select 4 node networks for monopositive states based on various cost functions

include("ModularEvolutionScript.jl")
include("analysisutils.jl")
include("CSBscript.jl")
using CSV, LinearAlgebra

nn = 4
ngoodnetworks = 20
costs = ["product", "sum", "nrmsd"]
cost = costs[3]

costfunc(J) = getPscore(calcFreq(solveCSB(J,100)),cost)

# J = rand(-1:1, nn, nn)
for _ in 1:100
    J = ones(Int,nn, nn) - I
    X, Xs = evolveNetBlind(J, costfunc)
    net, netscore = X[end], Xs[end]
    plot(Xs)
        
    topofile = interaction2topo(net)
    CSV.write("BestNetworks/fromFullyConnected$(cost)_$(round(rand(),digits=5)).topo",
        topofile,
        delim='\t')
end
