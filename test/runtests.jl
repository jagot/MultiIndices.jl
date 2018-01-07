using MultiIndices
using Base.Test

@testset "abstract" begin
    n = 1000
    𝕽¹ = Cartesian(n)

    𝕽̃¹ = 𝕽¹ ⊕ 𝕽¹ ⊕ 𝕽¹
    @test dims(𝕽̃¹) == dims(𝕽¹) == 1
    @test order(𝕽̃¹) == order(𝕽¹) == (1,)
    @test size(𝕽̃¹) == 3size(𝕽¹)[1] == 3n

    @test 𝕽̃¹[1] == 𝕽¹
    @test 𝕽̃¹[1,2] == 2
    # @test 𝕽̃¹[2,1] == 1 + n

    𝕽² = 𝕽¹ ⊗ 𝕽¹
    𝕽³a = 𝕽² ⊗ 𝕽¹
    𝕽³b = 𝕽¹ ⊗ 𝕽²

    @test dims(𝕽¹) == 1
    @test order(𝕽¹) == (1,)

    @test dims(𝕽²) == 2
    @test order(𝕽²) == (1,2)
    @test size(𝕽²) == (n, n)

    @test order(𝕽¹,𝕽¹) == (1,2)
    @test order(𝕽²,𝕽²) == (1,2,3,4)

    @test dims(𝕽³a) == dims(𝕽³b) == 3
    @test order(𝕽³a) == order(𝕽³b) == (1, 2, 3)
    @test size(𝕽³a) == size(𝕽³b) == (n, n, n)

    # # Actually, we would like this to be true, as well, since the
    # # constituent spaces are the same
    # @test 𝕽³a == 𝕽³b

    @test 𝕽²[2,1] == 2
    @test 𝕽²[1,2] == 1+n

    @test 𝕽³a[2,1,1] == 2
    @test 𝕽³a[1,2,1] == 1 + n
    @test 𝕽³a[1,1,2] == 1 + n^2
end

@testset "cartesian" begin
    𝖃𝖄 = CartesianXY(100,500)
    @test 𝖃𝖄[2,1] == 2
    @test 𝖃𝖄[1,2] == 101
    @test dims(𝖃𝖄) == 2
    @test order(𝖃𝖄) == (1,2)
    @test size(𝖃𝖄) == (100, 500)

    𝖄𝖃 = CartesianYX(100,500)
    @test 𝖄𝖃[2,1] == 501
    @test 𝖄𝖃[1,2] == 2
    @test dims(𝖄𝖃) == 2
    @test order(𝖄𝖃) == (2,1)
    @test size(𝖄𝖃) == (100, 500)
end

@testset "polar" begin
    𝔓 = Polar{(1,2)}(100,10)
    @test dims(𝔓) == 2
    @test order(𝔓) == (1,2)
    @test size(𝔓) == (100,10)
    @test 𝔓[2,1] == 2
    @test 𝔓[1,2] == 101
end

@testset "spin" begin
    𝔖 = Spin{5//2}()
    @test dims(𝔖) == 1
    @test order(𝔖) == (1,)
    @test size(𝔖) == (6,)
    @test 𝔖[4] == 4
end
