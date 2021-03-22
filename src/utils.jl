function make_hermitian!(arr, old_size)

    fix_corner = -1 
    # we loop over all dimensions which need to be fixed with hermitian property
    for dim = 1:length(old_size)
        # if initial array is odd, we don't need to fix hermitian property
        # if we downsample, don't fix 
        if old_size[dim] % 2 == 1 || size(arr, dim) ≤ old_size[dim]
            continue
        else
            # increase corner fix counter
            fix_corner += 1 
            ind_l, ind_r = get_indices_around_center(size(arr, dim), old_size[dim])
            ind_r += 1
            arr_left_slice = slice(arr, dim, ind_l)
            arr_left_slice ./= 2
            s_inds = slice_indices(arr, dim, ind_r)
            
            region_reverse_indices = region_revert_inds(arr, old_size)
            region_assign_indices = copy(region_reverse_indices) 
            region_reverse_indices[dim] = 1:1
            region_assign_indices[dim] = ind_r:ind_r
            arr[region_assign_indices...] = conj.(reverse_all(arr_left_slice[region_reverse_indices...]))
            # arr[s_inds...] = conj.(arr_left_slice[reverse_inds...])
        end
    end
 ## corner was divided too often, fix it
    if fix_corner ≥ 0
        @show "lol"
        corner = get_indices_around_center.(size(arr), old_size)
        corner_l = map(first, corner)
        corner_r = map(x-> 1 + x[2], corner)
        arr[corner_l...] *= 2 ^fix_corner
        arr[corner_r...] = conj(arr[corner_l...])
    end
    return arr
end

"""
    make_hermitian(arr)

Takes an array `arr` and appends rows, cols, ... if necessary
so that `arr` is a hermitian array which preserves Parseval's theorem.

# Examples
```jldoctest
julia> FFTResampling.make_hermitian([1.0 2.0])
1×3 Array{Float64,2}:
 0.5  2.0  0.5

julia> FFTResampling.make_hermitian([1.0 2.0; 3.0 4.0])
3×3 Array{Float64,2}:
 0.5  1.0  0.0
 1.5  4.0  1.5
 0.0  1.0  0.5

julia> FFTResampling.make_hermitian([1im 2.0; 3.0im 4.0im; 5.0 6.0im])
3×3 Array{Complex{Float64},2}:
 0.0+0.5im  2.0+0.0im  2.5-0.0im
 0.0+1.5im  0.0+4.0im  0.0-1.5im
 2.5+0.0im  0.0+6.0im  0.0-0.5im
```
"""
function make_hermitian(arr::AbstractArray{T, N}) where {T, N}
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
    @inbounds for d = 1:N
        # special case when singleton dimension 
        if size(arr_new)[d] == 1 || size(arr)[d] % 2 == 1
            continue
        end
    
        # extract indices represent the left slice we need to copy
        extract_indices = slice_indices(arr_new, d, 1)
        y = arr_new[extract_indices...] ./ 2
        y_rev = reverse_all(y)
        
        # assign indices represent the right sliuce we need to assign the
        # conjugated and halved one
        assign_indices = slice_indices(arr_new, d, size(arr_out)[d])

        # assign the halved
        arr_out[extract_indices...] = y
        # assign the conjugated, halved one
        arr_out[assign_indices...] = conj.(y_rev)
    end

    return arr_out
end


"""
    get_indices_around_center(i_in, i_out)
A function which provides two output indices `i1` and `i2`
where `i2 - i1 = i_out`
The indices are chosen in a way that the set `i1:i2`
cuts the interval `1:i_in` in a way that the center frequency
stays at the center position.
Works for both odd and even indices
"""
function get_indices_around_center(i_in, i_out)
    if (mod(i_in, 2) == 0 && mod(i_out, 2) == 0 
     || mod(i_in, 2) == 1 && mod(i_out, 2) == 1) 
        x = (i_in - i_out) ÷ 2
        return 1 + x, i_in - x
    elseif mod(i_in, 2) == 1 && mod(i_out, 2) == 0
        x = (i_in - 1 - i_out) ÷ 2
        return 1 + x, i_in - x - 1 
    elseif mod(i_in, 2) == 0 && mod(i_out, 2) == 1
        x = (i_in - (i_out - 1)) ÷ 2
        return 1 + x, i_in - (x - 1)
    end
end


"""
    center_extract(arr, new_size_array; view=false)
Extracts a center of an array. 
`new_size_array` must be list of sizes indicating the output
size of each dimension. Centered means that a center frequency
stays at the center position. Works for even and uneven.
If `length(new_size_array) < length(ndims(arr))` the remaining dimensions
are untouched and copied.
`view=true` returns a view of `arr`. 
# Examples
```jldoctest
julia> FFTResampling.center_extract([1 2; 3 4], [1]) 
1×2 Array{Int64,2}:
 3  4

julia> FFTResampling.center_extract([1 2; 3 4], [1, 1])
1×1 Array{Int64,2}:
 4

julia> FFTResampling.center_extract([1 2 3; 3 4 5; 6 7 8], [2 2])
2×2 Array{Int64,2}:
 1  2
 3  4
```
"""
function center_extract(arr::AbstractArray, new_size_array; view=false)
    new_size_array = collect(new_size_array)

    # we construct two lists
    # the reason is, that we don't change higher dimensions which are not 
    # specified in new_size_array
    out_indices1 = [get_indices_around_center(size(arr)[x], new_size_array[x]) 
                    for x = 1:length(new_size_array)]
    
    out_indices1 = [x[1]:x[2] for x = out_indices1]
    
    # out_indices2 contains just ranges covering the full size of each dimension
    out_indices2 = [1:size(arr)[i] for i = (1 + length(new_size_array)):ndims(arr)]
    if view
        return @view arr[out_indices1..., out_indices2...]
    end
    return arr[out_indices1..., out_indices2...]
end


"""
    center_set!(arr_large, arr_small)
Puts the `arr_small` central into `arr_large`.
The convention, where the center is, is the same as the definition
as for FFT based centered.
Function works both for even and uneven arrays.
# Examples
```jldoctest
julia> FFTResampling.center_set!([1, 1, 1, 1, 1, 1], [5, 5, 5])
6-element Array{Int64,1}:
 1
 1
 5
 5
 5
 1
```
"""
function center_set!(arr_large, arr_small)
    out_is = []
    for i = 1:ndims(arr_large)
        a, b = get_indices_around_center(size(arr_large)[i], size(arr_small)[i])
        push!(out_is, a:b)
    end

    #rest = ones(Int, ndims(arr_large) - 3)
    arr_large[out_is...] = arr_small
    
    return arr_large
end


"""
    center_pos(x)
Calculate the position of the center frequency.
Size of the array is `x`
# Examples
```jldoctest
julia> FFTResampling.center_pos(3)
2
julia> FFTResampling.center_pos(4)
3
```
"""
function center_pos(x::Integer)
    # integer division
    return div(x, 2) + 1
end



"""
    slice(arr, dim, index)

Return a `N` dimensional slice (where one dimensions has size 1) of the N-dimensional `arr` at the index position
`index` in the `dim` dimension of the array.
It holds `size(out)[dim] == 1`.

# Examples
```jldoctest
julia> x = [1 2 3; 4 5 6; 7 8 9]
3×3 Array{Int64,2}:
 1  2  3
 4  5  6
 7  8  9

julia> FFTResampling.slice(x, 1, 1)
1×3 view(::Array{Int64,2}, 1:1, :) with eltype Int64:
 1  2  3

julia> FFTResampling.slice(x, 2, 3)
3×1 view(::Array{Int64,2}, :, 3:3) with eltype Int64:
 3
 6
 9

```

"""
function slice(arr::AbstractArray{T, N}, dim::Integer, index::Integer) where {T, N}
    inds = slice_indices(arr, dim, index)
    return @view arr[inds...]
end

function slice_indices(arr::AbstractArray{T, N}, dim::Integer, index::Integer) where {T, N}
    f(x) = x[1] == dim ? (index:index) : (1:x[2])
    inds = map(f, enumerate(size(arr)))
    return inds
end




"""
    reverse_all(arr; view=true)

Reverse an array `arr` over all dimensions.
`view=true` returns a view instead of a new array.

# Examples
```jldoctest
julia> FFTResampling.reverse_all([1 2 3 4 5])
1×5 Array{Int64,2}:
 5  4  3  2  1

julia> FFTResampling.reverse_all([1 2 3 4 5; 7 8 9 10 11])
2×5 Array{Int64,2}:
 11  10  9  8  7
  5   4  3  2  1

julia> FFTResampling.reverse_all([1; 2; 3; 4; 5])
5-element Array{Int64,1}:
 5
 4
 3
 2
 1
```
"""
function reverse_all(arr::AbstractArray; view=true)
    return @view arr[reverse_all_indices(arr)...]
end

function reverse_all_indices(arr)
    out = []
    for i in size(arr)
        push!(out, i:-1:1)
    end
    return out 
end



function dft_1D(arr)
    N = length(arr)
    X = zeros(ComplexF64, N)
    
    for k = 0:N-1
        s = zero(X[1])
        for n = 0:N-1
            s += arr[n + 1]  * exp(-1im*2*π*k *n/ N)
        end
        X[k+1] = s
    end
    return X
end


function region_revert_inds(arr, old_size)
    inds = [1:1 for i = 1:ndims(arr)] 
    for d = 1:ndims(arr)
        ind_l, ind_r = get_indices_around_center(size(arr, d), old_size[d])
        inds[d] = (ind_l+1):(ind_r+1)
    end
    return inds
end




"""
This function ensures the constraints a spectrum must have so that after
iffting the result is purely real.

For example, inspect `ffshift(fft(randn((6,6))))` and you see that the highest negative
frequencies obey a certain symmetry. This can be generalized to N dimensions.
If we cut the array so that it'll be even in the that d dimension (new_size[d] % 2 ==0)
we need to manually restore the symmetry to guarantee a purely real result.
"""
function add_high_frequencies(size_arr, arr_out_f, new_size, N)
    fix_mul = false
    for d = 1:N
        if new_size[d] % 2 == 1 || new_size[d] == size_arr[d] ||
            (size_arr[d] % 2 == 0 && new_size[d] == (size_arr[d] + 1))
            continue
        end
        
        # construct the slice containing the highest positive frequency
        inds_extract = []
        inds_assign = []

        for i = 1:N
            a,b = get_indices_around_center(size_arr[i], new_size[i])
            if i == d
                # b+1 is the highest positive frequency which 
                # will be cut off (under certain conditions)
                iend = min(b+1, size_arr[i])
                push!(inds_extract, iend:iend)
                push!(inds_assign, a:a)
            else
                push!(inds_extract, a:min(b+1, size_arr[i]))
                push!(inds_assign, a:min(b+1, size_arr[i]))
            end
        end

        # add the highest positive frequency slice to the highest negative
        # but use the same arr_out_f
        # in that way we create a matrix which has the demanded symmetry of FFT to produce a real result
        arr_out_f[inds_assign...] = (arr_out_f[inds_assign...] .+ arr_out_f[inds_extract...])
    end
    
    return arr_out_f
end

