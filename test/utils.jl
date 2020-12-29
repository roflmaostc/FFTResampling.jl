
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

