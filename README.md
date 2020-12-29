# FFTInterpolations.jl

| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |


This package provides a simple sinc interpolation package written in Julia.
The FFTW based method `sinc_interpolate` requires a periodic, bandwidth limited (and properly Nyquist sampled) signal.


## Installation
`FFTInterpolations.jl` is available for all version equal or above Julia 1.0. It is mainly tested under Linux but should also work on Windows.
It can be installed with the following command

```julia
julia> ] add https://github.com/roflmaostc/FFTInterpolations.jl
```


[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg 
[docs-dev-url]: https://roflmaostc.github.io/FFTInterpolations.jl/dev/ 

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg 
[docs-stable-url]: https://roflmaostc.github.io/FFTInterpolations.jl/stable/

[CI-img]: https://github.com/roflmaostc/FFTInterpolations.jl/workflows/CI/badge.svg
[CI-url]: https://github.com/roflmaostc/FFTInterpolations.jl/actions?query=workflow%3ACI 

[codecov-img]: https://codecov.io/gh/roflmaostc/FFTInterpolations.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/roflmaostc/FFTInterpolations.jl
