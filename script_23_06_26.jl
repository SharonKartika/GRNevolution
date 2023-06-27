#select 4 node networks for monopositive states based on various cost functions

include("ModularEvolutionScript.jl")
include("analysisutils.jl")
include("CSBscript.jl")
using CSV

nn = 4
ngoodnetworks = 20
costs = ["product", "sum", "nrmsd"]
cost = costs[3]

costfunc(J) = getPscore(calcFreq(solveCSB(J,100)),cost)

J = rand(-1:1, nn, nn)
X, Xs = evolveNetBlind(J, costfunc)
net, netscore = X[end], Xs[end]

topofile = interaction2topo(net)
CSV.write("BestNetworks/$(cost)_$(round(rand(),digits=5)).topo",
    topofile,
    delim='\t')

