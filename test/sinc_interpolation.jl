@testset "sinc interpolation fft" begin

    function test_interpolation_sum_fft(N_low, N)
	    x_min = 0.0
	    x_max = 16π
	    
	    xs_low = range(x_min, x_max, length=N_low+1)[1:N_low]
	    xs_high = range(x_min, x_max, length=N)[1:end-1]
	    f(x) = sin(0.5*x) + cos(x) + cos(2 * x) + sin(0.25*x)
	    arr_low = f.(xs_low)
	    arr_high = f.(xs_high)

	    xs_interp = range(x_min, x_max, length=N+1)[1:N]
	    arr_interp = sinc_interpolate(arr_low, N)

	    xs_interp_s = range(x_min, x_max, length=N+1)[1:N]
	    arr_interp_s = FFTInterpolations.sinc_interpolate_sum(arr_low, N)

        @test ≈(arr_interp[2*N ÷10: N*8÷10], arr_high[2* N ÷10: N*8÷10], rtol=0.05)
        @test ≈(arr_high[2*N ÷10: N*8÷10], arr_interp_s[2*N ÷10: N*8÷10], rtol=0.05)
    end

    test_interpolation_sum_fft(128, 1000)
    test_interpolation_sum_fft(129, 1000)
    test_interpolation_sum_fft(120, 1531)
    test_interpolation_sum_fft(121, 1211)
end


@testset "downsampe" begin
    function test_downsample(N_low, N)
	    x_min = 0.0
	    x_max = 16π
	    
	    xs_low = range(x_min, x_max, length=N_low+1)[1:N_low]
	    f(x) = sin(0.5*x) + cos(x) + cos(2 * x) + sin(0.25*x)
	    arr_low = f.(xs_low)

	    xs_interp = range(x_min, x_max, length=N+1)[1:N]
	    arr_interp = sinc_interpolate(arr_low, N)

	    xs_interp_s = range(x_min, x_max, length=N+1)[1:N]
	    arr_interp_s = FFTInterpolations.sinc_interpolate_sum(arr_low, N)

        arr_ds = FFTInterpolations.downsample(arr_interp, N_low)
        @test ≈(arr_ds, arr_low)
    end

    test_downsample(128, 1000)
    test_downsample(128, 1232)
    test_downsample(128, 255)
    test_downsample(253, 254)
    test_downsample(253, 1001)
    test_downsample(99, 100101)
end


