### A Pluto.jl notebook ###
# v0.12.18

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
using Revise, FFTInterpolations, Images, TestImages, Interpolations, PlutoUI, FFTW

# ╔═╡ aec408dc-490a-11eb-0d5f-b167a9599994
md"#### Image Interpolation with FFT based sinc interpolation"

# ╔═╡ 18d15d82-50ee-11eb-1a86-abae54ad0050
begin
	x = range(-10.0, 10.0, length=129)[1:end-1]
	x_exact = range(-10.0, 10.0, length=2049)[1:end-1]
	y = x'
	y_exact = x_exact'
	arr = sinc.(sqrt.(x .^2 .+ y .^2))
	arr_exact = sinc.(sqrt.(x_exact .^2 .+ y_exact .^2))
	arr_interp = sinc_interpolate(arr[1:end, 1:end], (131, 131));
	arr_interp2 = sinc_interpolate(arr[1:end, 1:end], (512, 512));
	arr_interp3 = sinc_interpolate(arr[1:end, 1:end], (1024, 1024));
	arr_ds = downsample(arr_interp, (128, 128))
end

# ╔═╡ 12d67b86-5028-11eb-1f49-1b1d3b86b4cd
colorview(Gray, arr_ds)

# ╔═╡ 8999c1ec-5124-11eb-233b-93c53f928a04
≈(arr_ds, arr, rtol=1e-15)

# ╔═╡ dfbfa426-5049-11eb-140a-e9d9cf94ccbe
colorview(Gray, arr_interp3)

# ╔═╡ 12bff398-5028-11eb-04e7-5708064794a6
md"#### Napari.jl is a wrapper for a Python Image viewer called Napari"

# ╔═╡ 12a83f3c-5028-11eb-0505-2f152e052077
begin
#using Napari
#napari.view_image(arr)
	#napari.view_image(arr_interp)
	#napari.view_image(arr_interp2)
	#napari.view_image(arr_interp3)
	#napari.view_image(permutedims(cat(arr, arr_ds, arr .- arr_ds, dims=3), [3, 2, 1]))
end

# ╔═╡ d22ca4da-490a-11eb-1334-cf3ffdc24123
md"#### Cube interpolatated image"

# ╔═╡ 966033f4-490a-11eb-0b33-c7fff58b19e2
begin
	#k = BSpline(Constant())
	k2 = BSpline(Cubic(Line(OnGrid())))
	itp = interpolate(arr_interp, (k2, k2))
	arr_c_i = itp(LinRange(1, size(arr)[1], 512), LinRange(1, size(arr)[2], 512))
	img_c_i = colorview(Gray, arr_c_i)
end

# ╔═╡ fd0a8c24-504b-11eb-3e44-75fb8475f9da
md"### Slider to switch between images
The movement of the image originates from the reason that also in images the coordinate system changes. See the 1D plot example as reference.
"

# ╔═╡ 9245616e-490b-11eb-1648-c56a0155a111
img_all = cat(img_c_i, arr_interp2, dims=3);

# ╔═╡ 4c3e9134-490b-11eb-03a4-4d58d4d942da
md"
$(@bind index Slider(1:2))
"

# ╔═╡ 6a0cbf10-490b-11eb-10bd-359a3425da42
begin
	colorview(Gray, img_all[200:300, 200:300, index])
end

# ╔═╡ Cell order:
# ╠═64af5a38-490a-11eb-34bc-37017a602974
# ╟─aec408dc-490a-11eb-0d5f-b167a9599994
# ╠═18d15d82-50ee-11eb-1a86-abae54ad0050
# ╠═12d67b86-5028-11eb-1f49-1b1d3b86b4cd
# ╠═8999c1ec-5124-11eb-233b-93c53f928a04
# ╠═dfbfa426-5049-11eb-140a-e9d9cf94ccbe
# ╟─12bff398-5028-11eb-04e7-5708064794a6
# ╠═12a83f3c-5028-11eb-0505-2f152e052077
# ╟─d22ca4da-490a-11eb-1334-cf3ffdc24123
# ╠═966033f4-490a-11eb-0b33-c7fff58b19e2
# ╟─fd0a8c24-504b-11eb-3e44-75fb8475f9da
# ╠═9245616e-490b-11eb-1648-c56a0155a111
# ╟─4c3e9134-490b-11eb-03a4-4d58d4d942da
# ╠═6a0cbf10-490b-11eb-10bd-359a3425da42
