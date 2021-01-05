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
    center_extract(arr, new_size)
Extracts a center of an array. 
`new_size` must be list of sizes indicating the output
size of each dimension. Centered means that a center frequency
stays at the center position. Works for even and uneven.
If `length(new_size) < length(size(arr))` the remaining dimensions
are untouched and copied.
# Examples
```julia-repl
julia> center_extract([[1,2] [3, 4]], [1])
1×2 Array{Int64,2}:
 2  4
julia> center_extract([[1,2] [3, 4]], [1, 1])
1×1 Array{Int64,2}:
4
```
"""
function center_extract(arr, index_arrays)
    index_arrays = collect(index_arrays)

    # we construct two lists
    # the reason is, that we don't change higher dimensions which are not 
    # specified in index_arrays
    out_indices1 = [get_indices_around_center(size(arr)[x], index_arrays[x]) 
                    for x = 1:length(index_arrays)]
    
    out_indices1 = [x[1]:x[2] for x = out_indices1]

    out_indices2 = [1:size(arr)[length(out_indices1) + i] for i = (1 + size(index_arrays)[1]):ndims(arr)]
    return arr[out_indices1..., out_indices2...]
end


"""
    center_set!(arr_large, arr_small)
Puts the `arr_small` central into `arr_large`.
The convention, where the center is, is the same as the definition
as for FFT based centered.
Function works both for even and uneven arrays.
# Examples
```julia-repl
julia> center_set!([1, 1, 1, 1, 1, 1], [5, 5, 5])
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
```julia-repl
julia> center_pos(3)
2
julia> center_pos(4)
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
```julia-repl
julia> x = randn((3, 3))
3×3 Array{Float64,2}:
 -0.925331   1.79456     0.465846
 -0.18492    0.705636    0.328199
  1.36222   -0.0132336  -0.586589

julia> FFTInterpolations.slice(x, 2, 2)
3×1 Array{Float64,2}:
  1.794557861500336
  0.7056355732334497
 -0.013233577444161712

julia> FFTInterpolations.slice(x, 1, 1)
1×3 Array{Float64,2}:
 -0.925331  1.79456  0.465846
```

"""
function slice(arr::AbstractArray{T, N}, dim::Integer, index::Integer) where {T, N}
    @assert 1 ≤ dim ≤ N

    inds = fill(1:1, N)

    for (i, v) in enumerate(size(arr))
        inds[i] = (dim == i ? (index:index) : (1:v))
    end

    return arr[inds...]
end
