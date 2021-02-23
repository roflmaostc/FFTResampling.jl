export sinc_interpolate, sinc_interpolate_sum, downsample
export resample


"""
    resample(arr, new_size [, normalize])

Calculates the `sinc` interpolation of an `arr` on a new array size
`new_size`.
This method is based on FFTs and therefore implicitly assumes periodic
boundaries and a finite frequeny support.
`normalize=true` by default multiplies by an appropriate factor so that 
the average intensity stays the same.
If `size(new_size)[i] > size(arr)[i]`, we apply zero padding in Fourier space.
If `size(new_size)[i] < size(arr)[i]`, we cut out a centered part of the
Fourier spectrum.
We apply some tricks at the boundary to increase accuracy of highest frequencies. 

# Examples
```jldoctest
julia> resample([1.0, 2.0, 3.0, 4.0], 8)
8-element Array{Float64,1}:
 1.0
 1.085786437626905
 2.0
 2.5
 3.0
 3.914213562373095
 4.0
 2.5

julia> resample([1.0  2.0; 3.0 4.0], (4,4))
4×4 Array{Float64,2}:
 1.0  1.5  2.0  1.5
 2.0  2.5  3.0  2.5
 3.0  3.5  4.0  3.5
 2.0  2.5  3.0  2.5

julia> downsample([1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [3])
3-element Array{Float64,1}:
 0.5
 0.5
 0.5

julia> downsample([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [6])
6-element Array{Float64,1}:
  1.0
 -0.3333333333333333
  1.0
 -0.3333333333333333
  1.0
 -0.3333333333333333

julia> resample([1 2 3; 4 5 6], (3, 2))
3×2 Array{Float64,2}:
 1.0   3.0
 3.25  5.25
 3.25  5.25
```
"""
function resample(arr::AbstractArray{T, N}, new_size, normalize=true; take_real=true ) where {T<:Real, N}
    # go to fourier space
    arr_f = fftshift(fft(arr))
    # create fourier space new array
    
    # the idea is the following
    # 1) first we handle all the upsampling dimensions. By doing so, we leave the size
    # or increase it. downsampling is done afterwards. Hence we use max to create 
    # a new array which is a least as large as the initial one or larger
    # 2) we then handle all the downsampling dimensions
    new_size_interp = Tuple(max(x[1], x[2]) for x in zip(size(arr), new_size))
    out_f = zeros(eltype(arr_f), new_size_interp)


    # change arr_f to be a hermitian array because we want a purely real result
    # after iffting
    arr_f = make_hermitian(arr_f)
    # it can happen, that arr_f was now padded with an extra row
    # but there is no extra column in out_f
    # therefore, cut it
    inds = []
    for (a, o) in zip(size(arr_f), size(out_f))
        push!(inds, 1:min(a, o))
    end
    arr_f = arr_f[inds...]
    
    # set the old array into the new 0-padded array
    center_set!(out_f, arr_f)
    arr_f = out_f

    # if the new_size[d] is even, we need to add the highest positive frequency
    # of the initial spectrum
    # to the highest negative one. In that way, we get a purely real result
    arr_out_f = add_high_frequencies(size(arr_f), arr_f, new_size, N)
    # return arr_out_f

    # do the cutting in Fourier space
    arr_f_n = center_extract(arr_out_f, new_size)
    # back to real space 
    arr_out = ifft(ifftshift(arr_f_n))
    

    if normalize
        arr_out .*= length(arr_out) ./ length(arr)
    end
    if take_real
        return real(arr_out)
    else
        return arr_out
    end
end

 # for a complex signal: split into real and imaginary part and solve for each 
 # part individually
function resample(arr::AbstractArray{T}, new_size, normalize=true) where T<:Complex
    # array of same shape but with Complex element type
    arr_r = resample(real(arr), new_size, normalize)
    arr_i = 1im .* resample(imag(arr), new_size, normalize)
    return arr_r .+ arr_i
end


"""
    sinc_interpolate(arr, new_size [, normalize])

Calculates the `sinc` interpolation of an `arr` on a new array size
`new_size`.
This method is based on FFTs and therefore implicitly assumes periodic
boundaries and a finite frequeny support.
`normalize=true` by default multiplies by an appropriate factor so that 
the average intensity stays the same.

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
function sinc_interpolate(arr::AbstractArray{T, N}, new_size, normalize=true; take_real=true) where {T<:Real, N}
    if typeof(new_size) <: Number
        @assert new_size ≥ size(arr)[1] && ndims(arr) == 1
    else
        @assert new_size ≥ size(arr)
    end 
    # go to fourier space
    arr_f = fftshift(fft(arr))
    # create fourier space new array
    out_f = zeros(eltype(arr_f), new_size)
   
    # change arr_f to be a hermitian array because we want a purely real result
    # after iffting
    arr_f = make_hermitian(arr_f)
    # it can happen, that arr_f was now padded with an extra row
    # but there is no extra column in out_f
    # therefore, cut it
    inds = []
    for (a, o) in zip(size(arr_f), size(out_f))
        push!(inds, 1:min(a, o))
    end
    arr_f = arr_f[inds...]
    
    # set the old array into the new 0-padded array
    center_set!(out_f, arr_f)
    # go back to real space and apply proper value scaling
    out = ifft(ifftshift(out_f)) 
    if normalize
        out .*= length(out_f) ./ length(arr)
    end
    if take_real
        return real(out)
    else
        return out
    end
end

 # for a complex signal: split into real and imaginary part and solve for each 
 # part individually
function sinc_interpolate(arr::AbstractArray{T}, new_size, normalize=true) where T<:Complex
    # array of same shape but with Complex element type
    arr_r = sinc_interpolate(real(arr), new_size, normalize)
    arr_i = 1im .* sinc_interpolate(imag(arr), new_size, normalize)
    return arr_r .+ arr_i
end




 # some remarks for real valued downsampling which are not backed up by literature
 # In the case when we want to downsample a odd sized array:
 # If we add the slices from the highest and positive and highest negative
 # frequencies, we get an real result after iffting.

"""
    downsample(arr, new_size [, normalize])

Downsample an array `arr` to the new size `new_size`.
This is calculated by cutting a centered frequency window from the frequency
spectrum and going back to real space with an `ifft`.
`normalize=true` by default multiplies by an appropriate factor so that 
the average intensity stays the same.

# Examples
```jldoctest
julia> downsample([1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [3])
3-element Array{Float64,1}:
 0.5
 0.5
 0.5

julia> downsample([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [6])
6-element Array{Float64,1}:
  1.0
 -0.3333333333333333
  1.0
 -0.3333333333333333
  1.0
 -0.3333333333333333
```
"""
function downsample(arr::AbstractArray{T, N}, new_size::P, normalize=true; take_real=true) where {T<:Real, N, P}
    if N == 1 && ~(P <:AbstractArray)
        new_size = collect(new_size)
    end
    arr_f = fftshift(fft(arr))
    # if the new_size[d] is even, we need to add the highest positive frequency
    # of the initial spectrum
    # to the highest negative one. In that way, we get a purely real result
    arr_out_f = add_high_frequencies(size(arr), arr_f, new_size, N)
    # return arr_out_f

    # do the cutting in Fourier space
    arr_f_n = center_extract(arr_out_f, new_size)
    # back to real space 
    arr_out = ifft(ifftshift(arr_f_n))
    if normalize
        arr_out .*= length(arr_out) ./ length(arr)
    end

    if take_real
        return real(arr_out)
    else
        return arr_out
    end
end
 
 # for a complex signal: split into real and imaginary part and solve for each 
 # part individually
function downsample(arr::AbstractArray{T}, new_size, normalize=true) where {T<:Complex}
    arr_r = downsample(real(arr), new_size, normalize)
    arr_i = 1im .* downsample(imag(arr), new_size, normalize)
    return arr_r .+ arr_i
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
