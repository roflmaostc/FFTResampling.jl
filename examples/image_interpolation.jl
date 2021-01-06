### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 64af5a38-490a-11eb-34bc-37017a602974
using Revise, FFTInterpolations, Images, TestImages, Interpolations, PlutoUI

# ╔═╡ aec408dc-490a-11eb-0d5f-b167a9599994
md"#### Image Interpolation with FFT based sinc interpolation"

# ╔═╡ 12f4c2d0-5028-11eb-37b4-b19ad9a0603d
begin
	x = range(-10.0, 10.0, length=129)[1:end-1]
	y = x'
	arr = log.(abs.(sinc.(sqrt.(x .^2 .+ y .^2))) .+1)
	arr_interp = sinc_interpolate(arr[1:end, 1:end], (200, 200));
	arr_interp2 = sinc_interpolate(arr[1:end, 1:end], (1024, 1024));
	arr_interp3 = sinc_interpolate(arr[1:end, 1:end], (4096, 4096));
end

# ╔═╡ 12d67b86-5028-11eb-1f49-1b1d3b86b4cd
colorview(Gray, arr)

# ╔═╡ dfbfa426-5049-11eb-140a-e9d9cf94ccbe
colorview(Gray, arr_interp3)

# ╔═╡ 12bff398-5028-11eb-04e7-5708064794a6
md"#### Napari.jl is a wrapper for a Python Image viewer called Napari"

# ╔═╡ 12a83f3c-5028-11eb-0505-2f152e052077
#begin
#napari.view_image(arr)
#	napari.view_image(arr_interp)
#	napari.view_image(arr_interp2)
#	napari.view_image(arr_interp3)
#end

# ╔═╡ e02c2626-504b-11eb-3dde-4d133fd0388c
md"### Interpolation of a normal image"

# ╔═╡ 7a1e1dc8-490a-11eb-0b27-2100ada566c0
begin
	img = testimage("fabio_gray_256");
	img_a = Float64.(img);
end

# ╔═╡ edb98f40-504b-11eb-2dca-9dcbc42ff24e
md"### Sinc Interpolated image"

# ╔═╡ 82c3e9bc-490a-11eb-3214-e198e103a989
img_s_i = colorview(Gray, sinc_interpolate(img_a, (4096, 4096)))

# ╔═╡ d22ca4da-490a-11eb-1334-cf3ffdc24123
md"#### Cube interpolatated image"

# ╔═╡ 966033f4-490a-11eb-0b33-c7fff58b19e2
begin
	k = BSpline(Constant())
	k = BSpline(Cubic(Line(OnGrid())))
	itp = interpolate(img_a, (k, k))
	arr_c_i = itp(LinRange(1, size(img_a)[1], 4096), LinRange(1, size(img_a)[2], 4096))
	img_c_i = colorview(Gray, arr_c_i)
end

# ╔═╡ fd0a8c24-504b-11eb-3e44-75fb8475f9da
md"### Slider to switch between images"

# ╔═╡ 9245616e-490b-11eb-1648-c56a0155a111
img_all = cat(img_s_i, img_c_i, dims=3);

# ╔═╡ 4c3e9134-490b-11eb-03a4-4d58d4d942da
md"
$(@bind index Slider(1:2))
"

# ╔═╡ 6a0cbf10-490b-11eb-10bd-359a3425da42
begin
	img_all[800:2500, 800:2500, index]
end

# ╔═╡ Cell order:
# ╠═64af5a38-490a-11eb-34bc-37017a602974
# ╠═aec408dc-490a-11eb-0d5f-b167a9599994
# ╠═12f4c2d0-5028-11eb-37b4-b19ad9a0603d
# ╠═12d67b86-5028-11eb-1f49-1b1d3b86b4cd
# ╠═dfbfa426-5049-11eb-140a-e9d9cf94ccbe
# ╟─12bff398-5028-11eb-04e7-5708064794a6
# ╠═12a83f3c-5028-11eb-0505-2f152e052077
# ╠═e02c2626-504b-11eb-3dde-4d133fd0388c
# ╠═7a1e1dc8-490a-11eb-0b27-2100ada566c0
# ╟─edb98f40-504b-11eb-2dca-9dcbc42ff24e
# ╠═82c3e9bc-490a-11eb-3214-e198e103a989
# ╟─d22ca4da-490a-11eb-1334-cf3ffdc24123
# ╠═966033f4-490a-11eb-0b33-c7fff58b19e2
# ╟─fd0a8c24-504b-11eb-3e44-75fb8475f9da
# ╟─9245616e-490b-11eb-1648-c56a0155a111
# ╟─4c3e9134-490b-11eb-03a4-4d58d4d942da
# ╟─6a0cbf10-490b-11eb-10bd-359a3425da42
