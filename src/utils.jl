using SymbolicUtils

"""
    split_freqs(freqs::Array, ububs::Int, dbubs::Int)
It is useful to have the contraction frequencies in a single list, 
since it allows more easily to go over all possible permutation. 
However, other functions need the frequencies to be split into up- and down-bubbles.
The function `split_freqs` splits a single frequency array `[μ1, μ2, ..., μl, ν1, ν2, ..., νr]` 
into an array `ububs = [μ1, μ2, ..., μl]` of up-bubbles and an array `dbubs = [ν1, ν2, ..., νr]` of down-bubbles.

Arguments:
    `freqs`  - unified array of frequencies, should be of the form [μ1, μ2, ..., μl, ν1, ν2, ..., νr]
    `l` - number of up-modes
    `r` - number of down-modes

Returns:
    freqs_up - array of only the up-modes
    freqs_dn - array of only the down-modes
"""
function split_freqs(freqs::Array, l::Int, r::Int)
    if length(freqs) < l + r
        error("Number of frequencies does not match number of modes!")
    end
    freqs_up = freqs[1:l]
    freqs_dn = freqs[l+1:l+r]
    return freqs_up, freqs_dn
end

"""
    count_modes(diagram)
Counts the total number of up-modes and down-modes in a given diagram. 
Due to the way the diagrams are constructed now, this will be deprecated soon.
"""
function count_modes(diagram)
    bubs = tuple(map(sum, zip(diagram...)))[1]
    return bubs[1], bubs[2] 
end