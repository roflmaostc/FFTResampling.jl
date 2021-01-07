
function center_test(x1, x2, x3, y1, y2, y3)
    arr1 = randn((x1, x2, x3))
    arr2 = zeros((y1, y2, y3))

    FFTInterpolations.center_set!(arr2, arr1)
    arr3 = FFTInterpolations.center_extract(arr2, (x1, x2, x3))
    @test arr1 == arr3
end

 # test center set and center extract methods
@testset "center methods" begin
    center_test(4, 4, 4, 6,7,4)
    center_test(5, 4, 4, 7, 8, 4)
    center_test(5, 4, 4, 8, 8, 8)
    center_test(6, 4, 4, 7, 8, 8)


    @test 1 == FFTInterpolations.center_pos(1)
    @test 2 == FFTInterpolations.center_pos(2)
    @test 2 == FFTInterpolations.center_pos(3)
    @test 3 == FFTInterpolations.center_pos(4)
    @test 3 == FFTInterpolations.center_pos(5)
    @test 513 == FFTInterpolations.center_pos(1024)

    @test FFTInterpolations.get_indices_around_center((5), (2)) == (2, 3)
    @test FFTInterpolations.get_indices_around_center((5), (3)) == (2, 4)
    @test FFTInterpolations.get_indices_around_center((4), (3)) == (2, 4)
    @test FFTInterpolations.get_indices_around_center((4), (2)) == (2, 3)
end




@testset "slice" begin
    
    x = randn((1,2,3,4))
    y = FFTInterpolations.slice(x, 2, 2)
    @test x[:, 2:2, :, :] == y

    x = randn((5,2,3,4))
    y = FFTInterpolations.slice(x, 1, 4)
    @test x[4:4, :, :, :] == y

    x = randn((5))
    y = FFTInterpolations.slice(x, 1, 5)
    @test x[5:5] == y

end


@testset "slice indices" begin
    x = randn((1,2,3))
    y = FFTInterpolations.slice_indices(x, 1, 1)
    @test y == [1:1, :, :]


    x = randn((20,4,20, 1, 2))
    y = FFTInterpolations.slice_indices(x, 2, 3)
    @test y == [:, 3:3, :, :, :]
end



@testset "reverse function" begin
    reverse_all = FFTInterpolations.reverse_all
    reverse_all_indices = FFTInterpolations.reverse_all_indices
    function test_reverse(x, y)
        x_rev = reverse_all(x)
        @test x == reverse_all(x_rev)
        @test x_rev == y
        
        @test x_rev == x[reverse_all_indices(x)...]
    end

    test_reverse([1,2,3], [3,2,1])
    test_reverse([1 2; 3 4], [4 3; 2 1])

    x = randn((2,3,4,1, 5,6))
    
    x_rev = x[end:-1:1, end:-1:1,end:-1:1,end:-1:1,end:-1:1, end:-1:1]
    test_reverse(x, x_rev)

    
    y = [3:-1:1, 5:-1:1, 1:-1:1, 20:-1:1]
    @test y == reverse_all_indices(randn((3,5,1,20)))
end



 # Just a reference implementation I don't want to delete
 # is slower, but works
@gg function make_hermitian_ref(arr::AbstractArray{T, N}) where {T, N}
    N2 = N - 1
    quote 
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
        for d = 1:$N
            # special case when singleton dimension 
            if size(arr_new)[d] == 1 || size(arr)[d] % 2 == 1
                continue
            end
            
            # the N2 dimensional loops, loop over this
            iter_nloops = map(x -> 1:x, collect(size_new))
            # delete at dimension d because we have only N2 dimensions to loop over
            deleteat!(iter_nloops, d)
            
            # @show iter_nloops
            @nloops $N2 i (k -> iter_nloops[k]) begin
                # indices to access the left slice
                # left slice is the existing slice which will be copied
                inds_l = convert(Array{Int},collect(@ntuple $N2 i))
                insert!(inds_l, d, 1)
                
                # right slice which will be copied from the left slice
                inds_r = similar(inds_l)
                
                # because of the hermitian property, we reverse the slices
                for j in 1:length(inds_r)
                    # in this case, we want to get the right hand slice
                    if j == d
                        inds_r[j] = size_new[j]
                    else
                        inds_r[j] = size_new[j] - inds_l[j] + 1
                    end
                end
                
                # assign the complex conjugate 
                arr_out[inds_l...] = 0.5 * arr_new[inds_l...]
                arr_out[inds_r...] = 0.5 * conj(arr_new[inds_l...])
            end
    
        end
        return arr_out
    end
end




@testset "make hermitian" begin
    a = FFTInterpolations.make_hermitian([1im 2.0; 3.0im 4.0im; 5.0 6.0im])

    ã = [0.0+0.5im  2.0+0.0im  2.5-0.0im; 0.0+1.5im  0.0+4.0im  0.0-1.5im; 2.5+0.0im  0.0+6.0im  0.0-0.5im]

    @test a ≈ ã
    
    a = FFTInterpolations.make_hermitian([1.0; 2.0; 3.0im])
    ã = [1.0; 2.0; 3.0im]
    @test a ≈ ã


    a = FFTInterpolations.make_hermitian([1.0+1.0im; 2.0; 3.0; 3.0im])
    ã = [0.5+0.5im; 2.0; 3.0; 3.0im; 0.5-0.5im]
    @test a ≈ ã

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
    
    function compare_ref(x)
        a = FFTInterpolations.make_hermitian(x)
        b = make_hermitian_ref(x)
        @test a ≈ b
    end


    compare_ref(randn((2, 1, 2,3)))
    compare_ref(randn((2, 1, 2,3, 2, 4, 8)))
    compare_ref(randn((1022, 4)))
    compare_ref(randn((1022, 1022)))
    compare_ref(randn((5, 4)))
end



@testset "test 1D dft" begin
    x = randn(1)
    @test fft(x) ≈ FFTInterpolations.dft_1D(x)
    
    x = randn(122)
    @test fft(x) ≈ FFTInterpolations.dft_1D(x)
    
    x = randn(123)
    @test fft(x) ≈ FFTInterpolations.dft_1D(x)
end
