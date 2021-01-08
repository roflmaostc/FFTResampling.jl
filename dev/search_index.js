var documenterSearchIndex = {"docs":
[{"location":"#FFTInterpolations.jl","page":"FFTInterpolations.jl","title":"FFTInterpolations.jl","text":"","category":"section"},{"location":"","page":"FFTInterpolations.jl","title":"FFTInterpolations.jl","text":"Here you can find the docstrings of all functions.","category":"page"},{"location":"","page":"FFTInterpolations.jl","title":"FFTInterpolations.jl","text":"sinc_interpolate\nsinc_interpolate_sum\nFFTInterpolations.downsample","category":"page"},{"location":"#FFTInterpolations.sinc_interpolate","page":"FFTInterpolations.jl","title":"FFTInterpolations.sinc_interpolate","text":"sinc_interpolate(arr, new_size [, normalize])\n\nCalculates the sinc interpolation of an arr on a new array size new_size. This method is based on FFTs and therefore implicitly assumes periodic boundaries and a finite frequeny support. normalize=true by default multiplies by an appropriate factor so that  the average intensity stays the same.\n\nExamples\n\njulia> sinc_interpolate([1.0, 2.0, 3.0, 4.0], 8)\n8-element Array{Float64,1}:\n 1.0\n 1.085786437626905\n 2.0\n 2.5\n 3.0\n 3.914213562373095\n 4.0\n 2.5\n\njulia> sinc_interpolate([1.0  2.0; 3.0 4.0], (4,4))\n4×4 Array{Float64,2}:\n 1.0  1.5  2.0  1.5\n 2.0  2.5  3.0  2.5\n 3.0  3.5  4.0  3.5\n 2.0  2.5  3.0  2.5\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.sinc_interpolate_sum","page":"FFTInterpolations.jl","title":"FFTInterpolations.sinc_interpolate_sum","text":"sinc_interpolate_sum(arr, new_length)\n\nCalculates the sinc interpolation of an 1D arr on a new array size new_size.  This method is slow, because of an explicit sum evalulation and not a FFT based evaluation.\n\nExamples\n\njulia> sinc_interpolate_sum([1.0, 2.0, 3.0, 4.0], 8)\n8-element Array{Float64,1}:\n 1.0\n 1.7825353626292277\n 2.0\n 2.1220659078919377\n 3.0\n 4.1592491794681985\n 4.0\n 2.0735615442829793\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.downsample","page":"FFTInterpolations.jl","title":"FFTInterpolations.downsample","text":"downsample(arr, new_size [, normalize])\n\nDownsample an array arr to the new size new_size. This is calculated by cutting a centered frequency window from the frequency spectrum and going back to real space with an ifft. normalize=true by default multiplies by an appropriate factor so that  the average intensity stays the same.\n\nExamples\n\njulia> downsample([1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [3])\n3-element Array{Float64,1}:\n 0.5\n 0.5\n 0.5\n\njulia> downsample([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [6])\n6-element Array{Float64,1}:\n  1.0\n -0.3333333333333333\n  1.0\n -0.3333333333333333\n  1.0\n -0.3333333333333333\n\n\n\n\n\n","category":"function"},{"location":"","page":"FFTInterpolations.jl","title":"FFTInterpolations.jl","text":"Additionally we have some utility functions:","category":"page"},{"location":"","page":"FFTInterpolations.jl","title":"FFTInterpolations.jl","text":"FFTInterpolations.make_hermitian\nFFTInterpolations.reverse_all\nFFTInterpolations.slice\nFFTInterpolations.center_extract\nFFTInterpolations.center_set!\nFFTInterpolations.center_pos\nFFTInterpolations.get_indices_around_center","category":"page"},{"location":"#FFTInterpolations.make_hermitian","page":"FFTInterpolations.jl","title":"FFTInterpolations.make_hermitian","text":"make_hermitian(arr)\n\nTakes an array arr and appends rows, cols, ... if necessary so that arr is a hermitian array which preserves Parseval's theorem.\n\nExamples\n\njulia> FFTInterpolations.make_hermitian([1.0 2.0])\n1×3 Array{Float64,2}:\n 0.5  2.0  0.5\n\njulia> FFTInterpolations.make_hermitian([1.0 2.0; 3.0 4.0])\n3×3 Array{Float64,2}:\n 0.5  1.0  0.0\n 1.5  4.0  1.5\n 0.0  1.0  0.5\n\njulia> FFTInterpolations.make_hermitian([1im 2.0; 3.0im 4.0im; 5.0 6.0im])\n3×3 Array{Complex{Float64},2}:\n 0.0+0.5im  2.0+0.0im  2.5-0.0im\n 0.0+1.5im  0.0+4.0im  0.0-1.5im\n 2.5+0.0im  0.0+6.0im  0.0-0.5im\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.reverse_all","page":"FFTInterpolations.jl","title":"FFTInterpolations.reverse_all","text":"reverse_all(arr)\n\nReverse an array arr over all dimensions.\n\nExamples\n\njulia> FFTInterpolations.reverse_all([1 2 3 4 5])\n1×5 Array{Int64,2}:\n 5  4  3  2  1\n\njulia> FFTInterpolations.reverse_all([1 2 3 4 5; 7 8 9 10 11])\n2×5 Array{Int64,2}:\n 11  10  9  8  7\n  5   4  3  2  1\n\njulia> FFTInterpolations.reverse_all([1; 2; 3; 4; 5])\n5-element Array{Int64,1}:\n 5\n 4\n 3\n 2\n 1\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.slice","page":"FFTInterpolations.jl","title":"FFTInterpolations.slice","text":"slice(arr, dim, index)\n\nReturn a N dimensional slice (where one dimensions has size 1) of the N-dimensional arr at the index position index in the dim dimension of the array. It holds size(out)[dim] == 1.\n\nExamples\n\njulia> x = [1 2 3; 4 5 6; 7 8 9]\n3×3 Array{Int64,2}:\n 1  2  3\n 4  5  6\n 7  8  9\n\njulia> FFTInterpolations.slice(x, 1, 1)\n1×3 view(::Array{Int64,2}, 1:1, :) with eltype Int64:\n 1  2  3\n\njulia> FFTInterpolations.slice(x, 2, 3)\n3×1 view(::Array{Int64,2}, :, 3:3) with eltype Int64:\n 3\n 6\n 9\n\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.center_extract","page":"FFTInterpolations.jl","title":"FFTInterpolations.center_extract","text":"center_extract(arr, new_size_array)\n\nExtracts a center of an array.  new_size_array must be list of sizes indicating the output size of each dimension. Centered means that a center frequency stays at the center position. Works for even and uneven. If length(new_size_array) < length(ndims(arr)) the remaining dimensions are untouched and copied.\n\nExamples\n\njulia> FFTInterpolations.center_extract([1 2; 3 4], [1]) \n1×2 Array{Int64,2}:\n 3  4\n\njulia> FFTInterpolations.center_extract([1 2; 3 4], [1, 1])\n1×1 Array{Int64,2}:\n 4\n\njulia> FFTInterpolations.center_extract([1 2 3; 3 4 5; 6 7 8], [2 2])\n2×2 Array{Int64,2}:\n 1  2\n 3  4\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.center_set!","page":"FFTInterpolations.jl","title":"FFTInterpolations.center_set!","text":"center_set!(arr_large, arr_small)\n\nPuts the arr_small central into arr_large. The convention, where the center is, is the same as the definition as for FFT based centered. Function works both for even and uneven arrays.\n\nExamples\n\njulia> FFTInterpolations.center_set!([1, 1, 1, 1, 1, 1], [5, 5, 5])\n6-element Array{Int64,1}:\n 1\n 1\n 5\n 5\n 5\n 1\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.center_pos","page":"FFTInterpolations.jl","title":"FFTInterpolations.center_pos","text":"center_pos(x)\n\nCalculate the position of the center frequency. Size of the array is x\n\nExamples\n\njulia> FFTInterpolations.center_pos(3)\n2\njulia> FFTInterpolations.center_pos(4)\n3\n\n\n\n\n\n","category":"function"},{"location":"#FFTInterpolations.get_indices_around_center","page":"FFTInterpolations.jl","title":"FFTInterpolations.get_indices_around_center","text":"get_indices_around_center(i_in, i_out)\n\nA function which provides two output indices i1 and i2 where i2 - i1 = i_out The indices are chosen in a way that the set i1:i2 cuts the interval 1:i_in in a way that the center frequency stays at the center position. Works for both odd and even indices\n\n\n\n\n\n","category":"function"}]
}
