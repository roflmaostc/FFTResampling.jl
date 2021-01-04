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


@testset "make hermitian" begin
    a = FFTInterpolations.make_hermitian([1im 2.0; 3.0im 4.0im; 5.0 6.0im])

    ã = [0.0+0.5im  2.0+0.0im  2.5-0.0im; 0.0+1.5im  0.0+4.0im  0.0-1.5im; 2.5+0.0im  0.0+6.0im  0.0-0.5im]

    @test a == ã
    
    a = FFTInterpolations.make_hermitian([1.0; 2.0; 3.0im])
    ã = [1.0; 2.0; 3.0im]
    @test a == ã


    a = FFTInterpolations.make_hermitian([1.0+1.0im; 2.0; 3.0; 3.0im])
    ã = [0.5+0.5im; 2.0; 3.0; 3.0im; 0.5-0.5im]
    @test a == ã

    function test_sum(s)
        a = randn(s)
        ã = FFTInterpolations.make_hermitian(a)
        @test sum(a) ≈ sum(ã)
    end    
    test_sum((12, 3, 12, 11, 1))
    test_sum((123, 312))
    test_sum((1))
    test_sum((2))
    test_sum((100, 101))
    test_sum((101, 101))
    test_sum((101, 100))
    test_sum((101, 101))

    

    function test_symmetry(s)
        a = randn(s) 
        b = FFTInterpolations.make_hermitian(a)

        for k1 = 1:size(b)[1]
            for k2 = 1:size(b)[2]
                for k3 = 1:size(b)[3]
                    if k1 != 1 && k2 != 1 && k3 != 1 
                        continue
                    end
                    @test b[k1, k2, k3] ≈ 
                        conj(b[size(b)[1] + 1 - k1, size(b)[2] + 1 - k2, size(b)[3] + 1 - k3])
                end
            end
        end
    end            
    
    test_symmetry((2,2,2))
    test_symmetry((1,1,1))
    test_symmetry((4,4,4))
    test_symmetry((12,8,12))

    function test_symmetry_4D(s)
        a = randn(s) 
        b = FFTInterpolations.make_hermitian(a)

        for k1 = 1:size(b)[1]
            for k2 = 1:size(b)[2]
                for k3 = 1:size(b)[3]
                    for k4 = 1:size(b)[4]
                        if k1 != 1 && k2 != 1 && k3 != 1 && k4 != 1
                            continue
                        end
                        @test b[k1, k2, k3, k4] ≈ conj(b[size(b)[1] + 1 - k1, size(b)[2] + 1 - k2, size(b)[3] + 1 - k3, size(b)[4] + 1 - k4])
                    end
                end
            end
        end
    end            

    test_symmetry_4D((2,4, 6, 8))


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


