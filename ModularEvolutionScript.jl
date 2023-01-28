using StatsBase

"""Scoring function. Takes """
function testFitnessFunction(J::AbstractMatrix)
    return sum(J)
end

function best(M, fitnessFunction::Function)
    scores = fitnessFunction.(M)
    i = argmax(scores)
    return M[i], scores[i]
end

function createMutant(J::AbstractMatrix; nr::Int=0)
    l = length(J)
    (nr == 0) ? (nr = l ÷ 10) : () #10% mutation rate 
    J[sample(1:l, nr, replace=false)] = rand(-1:1, nr)
    # ensure that all the nodes are connected.  
    # see if the network has empty crosses. If yes, add random element in the row.
    # nzeros = zeros(Int, n)
    # for i in 1:n
    #     if (J[i, :] == nzeros) && (J[:, i] == nzeros)
    #         J[i, rand(1:n)] = rand([-1, 1])
    #     end
    # end
    return J
end

function createMutants(J::AbstractMatrix;
                        nmutants::Int=6,
                        nr::Int=0)
    return createMutant.([copy(J) for _ in 1:nmutants], nr=nr)
end

"""
Generates mutants, selects best.
Maintains list of seen nodes: prevents backtracking.
Gets slow as `seen` gets larger (quadratically?).

Returns the sequence of `selected` Networks, 
and their corresponding `selectedScore`.
"""
function evolveNet(J::AbstractMatrix,
                    fitnessFunction::Function;
                    nmutants::Int=10,
                    nr::Int=0,
                    niter::Int=100)
    C = copy(J)
    seen = [C]
    selected = [C]
    selectedscore = [fitnessFunction(C)]
    for i in 1:niter
        M = createMutants(C, nmutants=nmutants, nr=nr)
        M = setdiff(M, seen)
        if length(M) == 0
            continue
        end
        C, Cs = best(M, fitnessFunction)
        push!(seen, M...)
        push!(selected, C)
        push!(selectedscore, Cs)
    end
    return selected, selectedscore
end

"""
Generates mutants, selects best. 
Can get stuck in cycles and local optima.

Returns the sequence of `selected` Networks, 
and their corresponding `selectedScore`.
"""
function evolveNetBlind(J::AbstractMatrix,
                        fitnessFunction::Function;
                        nmutants::Int=10,
                        nr::Int=0,
                        niter::Int=100)
    C = copy(J)
    selected = [C]
    selectedscore = [fitnessFunction(C)]
    for i in 1:niter
        M = createMutants(C, nmutants=nmutants, nr=nr)
        C, Cs = best(M, fitnessFunction)
        push!(selected, C)
        push!(selectedscore, Cs)
    end
    return selected, selectedscore
end

"""
Generates mutants, selects best.
Keeps track of mutants in the last `k` generations.
Prevents short term backtracking.  

Returns the sequence of `selected` Networks, 
and their corresponding `selectedScore`.
"""
function evolveNetTabu(J::AbstractMatrix,
                    fitnessFunction::Function;
                    nmutants::Int=10,
                    nr::Int=0,
                    niter::Int=100,
                    k::Int = 10)
    C = copy(J)
    seen = [C]
    selected = [C]
    selectedscore = [fitnessFunction(C)]
    for i in 1:niter
        M = createMutants(C, nmutants=nmutants, nr=nr)
        M = setdiff(M, seen)
        if length(M) == 0
            continue
        end
        C, Cs = best(M, fitnessFunction)
        push!(seen, M...)
        if (i > k)
            deleteat!(seen, 1:(length(seen) ÷ k))
        end
        print(length(seen), " ")
        push!(selected, C)
        push!(selectedscore, Cs)
    end
    return selected, selectedscore
end
        
