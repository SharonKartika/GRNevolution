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

function convertrange(ri, x1, x2, y1, y2)
    runit = (ri - x1) / (x2 - x1)
    rf = runit * (y2 - y1) + y1
    return rf
end

"index, start, finish, length"
function pgb(i, s, f, l=20)
    if (f-s+1) < l
        l = f-s+1
    end 
    j = (f-s+1) ÷ l
    if (i == s)
        println("|"*("-"^(l+1))*"|")
        print("|=")
    end
    if (i % j) == 0
        print("=")
        if (i == f)
            print("|\n")
        end    
    end
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
        pgb(i, 1, niter)
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
        
"""
Generates mutants, selects if score greater all previous scores.
Monotonous fitness increase. Maintains list of previous scores

Returns the sequence of `selected` Networks, 
and their corresponding `selectedScore`.
"""
function evolveNetMonotone(J::AbstractMatrix,
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
        if all(Cs .> selectedscore) 
            push!(selected, C)
            push!(selectedscore, Cs)
        end
        pgb(i, 1, niter)
    end
    return selected, selectedscore
end
