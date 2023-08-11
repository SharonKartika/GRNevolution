include("CSBscript.jl")
include("analysisutils.jl")
using CSV


filedir = "TopoFilesToAnalyze/"
files = readdir(filedir) 
Ds = CSV.File.("$filedir".*files) .|> DataFrame

@time for i in eachindex(Ds)
    D = Ds[i]
    J = topo2interaction(D)

    sol = solveCSB(J, 100000) |> calcFreq
    sol[!, "Frust"] = [calcFrust(J, sol[i, 1]) for i in 1:size(sol)[1]]
    sol
    rename!(sol, ["Sequence"=>"states", "RelFreq"=>"Avg0", "Frust"=>"frust0"])

    CSV.write("TopoFileResults/"*files[i][begin:end-5]*"_CSB_results.csv", sol)
end