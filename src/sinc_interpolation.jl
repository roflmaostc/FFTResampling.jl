export sinc_interpolate, sinc_interpolate_sum

using Base.Cartesian
using GeneralizedGenerated

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
    @assert new_size ≥ size(arr)

    # go to fourier space
    arr_f = fftshift(fft(arr))
    # create fourier space new array
    out_f = zeros(eltype(arr_f), new_size)
    # set the old array into the new 0-padded array
    center_set!(out_f, arr_f)
    # go back to real space and apply proper value scaling
    out = ifft(ifftshift(out_f)) ./ length(arr) .* length(out_f)
    return out
end

 # for real arrays, take real part due to numerical inaccuracies
function sinc_interpolate(arr::AbstractArray{T}, new_size) where T<:Real
    # array of same shape but with Complex element type
    arr_nt = Complex.(arr) 
    return real(sinc_interpolate(arr_nt, new_size))
end


"""
    make_hermitian(arr)

Takes an array `arr` and appends rows, cols, ... if necessary
so that `arr` is a hermitian array which preserves Parseval's theorem.

# Examples
```julia-repl
julia> FFTInterpolations.make_hermitian([1.0 2.0])
1×3 Array{Float64,2}:
 0.5  2.0  0.5

julia> FFTInterpolations.make_hermitian([1.0 2.0; 3.0 4.0])
3×3 Array{Float64,2}:
 0.5  1.0  0.0
 1.5  4.0  1.5
 0.0  1.0  0.5

julia> FFTInterpolations.make_hermitian([1im 2.0; 3.0im 4.0im; 5.0 6.0im])
3×3 Array{Complex{Float64},2}:
 0.0+0.5im  2.0+0.0im  2.5-0.0im
 0.0+1.5im  0.0+4.0im  0.0-1.5im
 2.5+0.0im  0.0+6.0im  0.0-0.5im
```
"""
@gg function make_hermitian(arr::AbstractArray{T, N}) where {T, N}
    N2 = N - 1
    quote 
        # if the size is odd, we need to add a slice to that dimension
        # at this slice we then fill new values, so that `arr` is hermitian
        mf(x) = x % 2 == 0 ? x +1 : x
        size_new = map(mf, size(arr))
        # new array with new size 
        arr_new = zeros(eltype(arr), size_new)
        
        # fill the new array with the old one 
        same_ind = map(x -> 1:x, size(arr))
        arr_new[same_ind...] = arr 
        
        # copy again to have a out array (might be necessary because of reassigning)
        arr_out = copy(arr_new)
        
        # now modify values at the edges to get hermitian property
        # and Parsveval's theorem correct
        for d = 1:$N
            # @show d
            # special case when singleton dimension 
            if size(arr_new)[d] == 1 || size(arr)[d] % 2 == 1
                continue
            end
            
            # the N2 dimensional loops, loop over this
            iter_nloops = map(x -> 1:x, collect(size_new))
            # delete at dimension d because we have only N2 dimensions to loop over
            deleteat!(iter_nloops, d)
            
            # @show iter_nloops
            @nloops $N2 i (k -> iter_nloops[k]) begin
                # indices to access the left slice
                # left slice is the existing slice which will be copied
                inds_l = collect(@ntuple $N2 i)
                insert!(inds_l, d, 1)
                
                # right slice which will be copied from the left slice
                inds_r = similar(inds_l)
                
                # because of the hermitian property, we reverse the slices
                for j in 1:length(inds_r)
                    # in this case, we want to get the right hand slice
                    if j == d
                        inds_r[j] = size_new[j]
                    else
                        inds_r[j] = size_new[j] - inds_l[j] + 1
                    end
                end
                
                # assign the complex conjugate 
                arr_out[inds_l...] = 0.5 * arr_new[inds_l...]
                arr_out[inds_r...] = 0.5 * conj(arr_new[inds_l...])
            end
    
        end
        return arr_out
    end
end




"""
    downsample(arr, new_size)

Downsample an array `arr` to the new size `new_size`.
This is calculated by cutting a centered frequency window from the frequency
spectrum and going back to real space with an `ifft`.

# Examples
```julia-repl
julia> FFTInterpolations.downsample([1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [3])
3-element Array{Float64,1}:
 0.5
 0.5
 0.5

julia> FFTInterpolations.downsample([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [6])
6-element Array{Float64,1}:
 0.6666666666666666
 0.0
 0.6666666666666666
 0.0
 0.6666666666666666
 0.0
```
"""
function downsample(arr::AbstractArray{T, N}, new_size::P) where {T<:Complex, N, P}
    if N == 1 && ~(P <:AbstractArray)
        new_size = [new_size]
    end
    arr_f = fftshift(fft(arr))
    arr_f_n = center_extract(arr_f, new_size)
    
    arr_n = real(ifft(ifftshift(arr_f_n))) ./ length(arr) .* length(arr_f_n)
    return arr_n
end

function downsample(arr::AbstractArray{T}, new_size) where {T<:Real}
    arr_nt = Complex.(arr) 
    return real(downsample(arr_nt, new_size))
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
