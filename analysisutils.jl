using SparseArrays
"""Returns a list of strings containing `1`s and `0`s,
which are monopositive"""
function getMonopositiveStrings(n)
    [("0"^(i - 1)) * "1" * ("0"^(n - i)) for i in 1:n]
end

function interaction2topo(tNet::AbstractMatrix)
    xi, yi, vi = findnz(sparse(tNet))
    lets = 'A':'Z'
    nmap(x) = (3 - x) ÷ 2 # to map 1->1 and -1->2
    df = DataFrame(Source=lets[xi], Target=lets[yi], Type=nmap.(vi))
    df
end

"""Takes in the result of simulation
(a list of strings of outputs), and finds the relative
frequencies of each unique state."""
function calcFreq(dfr)
    D = proportionmap(dfr)
    dfFreq = DataFrame(Sequence=collect(keys(D)), RelFreq=collect(values(D)))
    #=for when a monopositive state is not present 
    in the output of an evaluation.=#
    S = getMonopositiveStrings(length(dfr[1,1]))
    for i in S
        if !(i in dfFreq.Sequence)
            push!(dfFreq, (i, 0.0))
        end
    end 
    dfFreq
end

"find indices of elements of `a` that are present in `b`"
function ainb(a, b)
    [i for i in 1:length(a) if a[i] in b]
end


function getPscore(dfFreq, scoretype::String="sum")
    nn = length(dfFreq[1, 1]) #number of nodes in network
    reqStates = getMonopositiveStrings(nn) # list of required states
    indexOfReq = ainb(dfFreq.Sequence, reqStates)
    reqFreqs = dfFreq[indexOfReq, "RelFreq"]
    if (scoretype == "product")
        return prod(reqFreqs)
    elseif (scoretype == "sum")
        return sum(reqFreqs)
    elseif (scoretype == "productbysum")
        return prod(reqFreqs) / sum(reqFreqs)
    elseif (scoretype == "sumbyproduct")
        return sum(reqFreqs) / prod(reqFreqs)
    elseif (scoretype == "nrmsd")
        function ismonopositive(s::String)
            return s in reqStates
        end
        Dmono = filter(:Sequence => ismonopositive, dfFreq).RelFreq
        DNmono = filter(:Sequence => x -> !ismonopositive(x), dfFreq).RelFreq
        f = 1/nn 
        score = sum((Dmono .- f) .^ 2) + sum(DNmono .^ 2)
        return -sqrt(score)
    else
        error("The string passed does not match any known scoring function")
    end
end
