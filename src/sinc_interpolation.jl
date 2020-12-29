export sinc_interpolate, sinc_interpolate_sum

"""
    sinc_interpolate(arr, new_size)

Calculates the `sinc` interpolation of an `arr` on a new array size
`new_size`.
This method is based on FFTs and therefore implicitly assumes periodic
boundaries and a finite frequeny support.

# Examples
```julia-repl
julia> sinc_interpolate([1.0, 2.0, 3.0, 4.0], 8)
8-element Array{Float64,1}:
 1.0
 1.085786437626905
 2.0
 2.5
 3.0
 3.914213562373095
 4.0
 2.5

julia> sinc_interpolate([1.0  2.0; 3.0 4.0], (4,4))
4×4 Array{Float64,2}:
 1.0  1.5  2.0  1.5
 2.0  2.5  3.0  2.5
 3.0  3.5  4.0  3.5
 2.0  2.5  3.0  2.5
```
"""
function sinc_interpolate(arr::AbstractArray{T}, new_size) where T<:Complex
    arr_f = fftshift(fft(arr))
    out_f = zeros(eltype(arr_f), new_size)
    center_set!(out_f, arr_f)
    out = ifft(ifftshift(out_f)) ./ length(arr) .* length(out_f)
    return out
end

function sinc_interpolate(arr::AbstractArray{T}, new_size) where T<:Real
    # array of same shape but with Complex element type
    arr_nt = Complex.(arr) 
    return real(sinc_interpolate(arr_nt, new_size))
end


function downsample(arr::AbstractArray{T, N}, new_size::P) where {T<:Real, N, P}
    if N == 1 && ~(P <:AbstractArray)
        new_size = [new_size]
    end
    arr_f = fftshift(fft(arr))
    arr_f_n = center_extract(arr_f, new_size)
    
    arr_n = real(ifft(ifftshift(arr_f_n))) ./ length(arr) .* length(arr_f_n)
    return arr_n
end

"""
    sinc_interpolate_sum(arr, new_length)

Calculates the `sinc` interpolation of an 1D `arr` on a new array size
`new_size`. 
This method is slow, because of an explicit sum evalulation and
not a FFT based evaluation.

# Examples
```julia-repl
julia> sinc_interpolate_sum([1.0, 2.0, 3.0, 4.0], 8)
8-element Array{Float64,1}:
 1.0
 1.7825353626292277
 2.0
 2.1220659078919377
 3.0
 4.1592491794681985
 4.0
 2.073561544282979
```
"""
function sinc_interpolate_sum(arr, new_size)
    out = zeros(eltype(arr), new_size)
    T = 1 / size(arr)[1]
    t_arr = range(0, 1.0, length=size(out)[1]+1)[1:size(out)[1]]
    for (j, t) in enumerate(t_arr)
        v = zero(eltype(arr))
        for n = 1:size(arr)[1] 
            v += arr[n] * sinc((t - (n-1) * T) / T)
        end
        out[j] = v 
    end
    return out
end
