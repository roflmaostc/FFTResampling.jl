export sinc_interpolate, sinc_interpolate_sum

"""
    sinc_interpolate(arr, new_size)

Calculates the `sinc` interpolation of an `arr` on a new array size
`new_size`.
This method is based on FFTs and therefore implicitly assumes periodic
boundaries and a finite frequeny support.

# Examples
```jldoctest
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
function sinc_interpolate(arr::AbstractArray{T}, new_size; real=false) where T<:Complex
    if typeof(new_size) <: Number
        @assert new_size ≥ size(arr)[1] && ndims(arr) == 1
    else
        @assert new_size ≥ size(arr)
    end 
    # go to fourier space
    arr_f = fftshift(fft(arr))
    # create fourier space new array
    out_f = zeros(eltype(arr_f), new_size)
   
    # in this case, change arr_f to be a hermitian array
    if real
        arr_f = make_hermitian(arr_f)
        #it can happen, that arr_f was now padded with an extra row
        # but there is no extra column in out_f
        # therefore, cut it
        inds = []
        for (a, o) in zip(size(arr_f), size(out_f))
            push!(inds, 1:min(a, o))
        end
        arr_f = arr_f[inds...]
    end
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
    return real(sinc_interpolate(arr_nt, new_size, real=true))
end




 # some remarks for real valued downsampling which are not backed up by literature
 # In the case when we want to downsample a odd sized array
 # we need to erase a single slice of frequencies.
 # But iffting yields a complex and non real result.
 # However, if we add the slices from the highest and positive and highest negative
 # frequencies, we get an real result after iffting.

"""
    downsample(arr, new_size)

Downsample an array `arr` to the new size `new_size`.
This is calculated by cutting a centered frequency window from the frequency
spectrum and going back to real space with an `ifft`.

# Examples
```jldoctest
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
        new_size = collect(new_size)
    end
    arr_f = fftshift(fft(arr))
    
    # if the new_size[d] is even, we need to add the highest positive frequency
    # of the initial spectrum
    # to the highest negative one. In that way, we get a purely real result
    for d = 1:N
        if new_size[d] % 2 == 1 || new_size[d] == size(arr)[d] ||
            (size(arr)[d] % 2 == 0 && new_size[d] == (size(arr)[d] + 1))
            continue
        end
        
        # construct the slice containing the highest positive frequency
        inds_extract = []
        inds_assign = []
        for i = 1:N
            a,b = get_indices_around_center(size(arr)[i], new_size[i])
            if i == d
                # b+1 is the highest positive frequency which 
                # will be cut off (under certain conditions)
                push!(inds_extract, min(b+1, size(arr)[i]))
                push!(inds_assign, a)
            else
                push!(inds_extract, a:min(b+1, size(arr)[i]))
                push!(inds_assign, a:min(b+1, size(arr)[i]))
            end
        end
        # add the highest positive frequency slice to the highest negative
        arr_f[inds_assign...] += arr_f[inds_extract...]
    end
   
    # do the cutting in Fourier space
    arr_f_n = center_extract(arr_f, new_size)
    # back to real space and renormalization
    arr_n = ifft(ifftshift(arr_f_n)) ./ length(arr) .* length(arr_f_n)
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
```jldoctest
julia> sinc_interpolate_sum([1.0, 2.0, 3.0, 4.0], 8)
8-element Array{Float64,1}:
 1.0
 1.7825353626292277
 2.0
 2.1220659078919377
 3.0
 4.1592491794681985
 4.0
 2.0735615442829793
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
