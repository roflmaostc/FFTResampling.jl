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


@testset "sinc downsample" begin
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

    
    function test_real(s, s_n)
        x = randn(s)
        x = x .+ 0.0im
        y = FFTInterpolations.downsample(x, s_n)
        @test all(imag.(y) .< 1e-14)
    end
    
    test_real((12, 13, 14, 15), (12, 12, 12, 12))
    test_real((12, 13, 14, 15), (12, 13, 13, 13))
    test_real((12, 13, 14, 15), (12, 12, 13, 14))

    test_real((7,7), (6,6))
    test_real(100, 99)
    test_real(101, 89)
    test_real(100, 98)
    test_real(101, 90)


end



@testset "FFT sinc downsample and upsample together in 2D" begin


    function test_2D(in_s, out_s)
        x = range(-10.0, 10.0, length=in_s[1] + 1)[1:end-1]
        y = range(-10.0, 10.0, length=in_s[2] + 1)[1:end-1]'
	    arr = sinc.(sqrt.(x .^2 .+ y .^2))
	    arr_interp = sinc_interpolate(arr[1:end, 1:end], out_s);
	    arr_ds = downsample(arr_interp, in_s)
        @test arr_ds ≈ arr
    end

    test_2D((128, 128), (150, 150))
    test_2D((128, 128), (151, 151))
    test_2D((129, 129), (150, 150))
    test_2D((129, 129), (151, 151))
    
    test_2D((150, 128), (151, 150))
    test_2D((128, 128), (151, 153))
    test_2D((129, 128), (150, 153))
    test_2D((129, 128), (129, 153))


    x = range(-10.0, 10.0, length=129)[1:end-1]
    x2 = range(-10.0, 10.0, length=130)[1:end-1]
    x_exact = range(-10.0, 10.0, length=2049)[1:end-1]
    y = x'
    y2 = x2'
    y_exact = x_exact'
    arr = sinc.(sqrt.(x .^2 .+ y .^2))
    arr2 = sinc.(sqrt.(x .^2 .+ y .^2))
    arr_exact = sinc.(sqrt.(x_exact .^2 .+ y_exact .^2))
    arr_interp = sinc_interpolate(arr[1:end, 1:end], (131, 131));
    arr_interp2 = sinc_interpolate(arr[1:end, 1:end], (512, 512));
    arr_interp3 = sinc_interpolate(arr[1:end, 1:end], (1024, 1024));
    arr_ds = downsample(arr_interp, (128, 128))
    arr_ds2 = downsample(arr_interp, (128, 128))
    arr_ds23 = downsample(arr_interp2, (512, 512))
    arr_ds3 = downsample(arr_interp, (128, 128))

    @test ≈(arr_ds3, arr)
    @test ≈(arr_ds2, arr)
    @test ≈(arr_ds, arr)
    @test ≈(arr_ds23, arr_interp2)

end
