# take a network as input
# simulate it 
# return output
using Plots, DataFrames


"""constrains x between -1 and 1"""
function constrain(x, a::Float64=-1.0, b::Float64=1.0)
    if x > b
        return b
    elseif x < a
        return a
    else
        return x
    end
end

σ(x) = ((1 / (1 + exp(-10x))) - 0.5) * 2

"""Take a network, evaluates it, returns the trajectory"""
function simulateCSBtraj(J, N=100)
    Δt = 0.1 #time step
    nn = size(J, 2) #number of nodes
    X = rand(-1:1.0:1, nn) #expression levels of each node 
    XN = zeros(nn, N)
    A = zeros(nn)
    for k in 1:N
        A = J'X
        X = X + (Δt .* A) ./ (abs.(A) .+ 1)
        X = constrain.(X)
        XN[:, k] = X
    end
    return XN
end


#TODO: rewrite as a special case of simulateCSBtraj
"""Take a network, evaluates it, returns the trajectory"""
function simulateCSB(J, asString=true)
    N = 1000 #number of steps
    Δt = 0.1 #time step
    nn = size(J, 2) #number of nodes
    X = rand(-1.0:1.0, nn) #expression levels of each node 
    A = zeros(nn)
    for _ in 1:N
        A = J'X
        X = X + (Δt .* A) ./ (abs.(A) .+ 1)
        X = constrain.(X)
    end
    Xint = (x -> constrain(x, 0.0, 1.0)).(X) .|> round .|> Int
    if asString
        Xint = Xint .|> string
        return *(Xint...)
    end
    return Xint
end

solveCSB(J, n=1000; asString=true) = [simulateCSB(J, asString) for _ in 1:n]

function averageOutput(J, n=1000)
    R = solveCSB(J, n, asString=false)
    s = zeros(size(R[1]))
    for r in R
        s += r
    end
    s ./= length(R)
    return s
end

function getMonoScores(J, evalfunc)
    D1 = cFreq(evalfunc(J))
    S = getMonopositiveStrings(size(J, 1))
    f(x) = x in S
    D1 = D1[f.(D1.Sequence), :]
    sort!(D1)
    D1
end

function getAllBinaryStrings(n)
    x = String[]
    function f(s, x, n, i)
        if i < n
            f(s * "0", x, n, i + 1)
            f(s * "1", x, n, i + 1)
        else
            push!(x, s)
        end
    end
    return f("", x, n, 0)
end

J = rand(-1:1, 3, 3)
x = rand(0:0.001:1, size(J)[1])

# can be calculated for a network, state pair 
function calcFrust(J, x::Vector{Float64})
    n, m = size(J)
    Ef = 0 # #frustrated edges
    Et = sum(abs.(J)) # #edges (nonzero)
    for i in 1:n, j in 1:m  
        if J[i, j] * x[i] * x[j] < 0
            Ef += 1
        end
    end
    Ef / Et 
end

exprLevelsStringToVec(x::String) = split(x, "") .|> x->parse(Float64, x)

calcFrusta(J, x::String) = calcFrust(J, exprLevelsStringToVec(x))