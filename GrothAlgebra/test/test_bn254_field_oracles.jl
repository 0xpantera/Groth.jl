using Test
using GrothAlgebra
using Random

serialize_field(x::Union{BN254Fq,BN254Fr}) = string(convert(BigInt, x))

@testset "BN254 field oracle vectors" begin
    big(s::AbstractString) = parse(BigInt, s)

    fq_a = bn254_fq(big("1234567890123456789012345678901234567890"))
    fq_b = bn254_fq(big("987654321098765432109876543210987654321"))
    fr_a = bn254_fr(big("1357913579135791357913579135791357913579"))
    fr_b = bn254_fr(big("2468024680246802468024680246802468024680"))

    @test serialize_field(fq_a * fq_b) == "15472953419057815038298011349517756066574095001981990885326662174623684880625"
    @test serialize_field(inv(fq_a)) == "21353693958014289373041527313235984653524774391338635633977575809931291938135"
    @test fq_a * inv(fq_a) == one(BN254Fq)

    @test serialize_field(fr_a * fr_b) == "2463067557993456352666791777416259958222174884297511614703693412527579300319"
    @test serialize_field(inv(fr_a)) == "13691281592788382201413824496711925257783763687512458266092765078913737141380"
    @test fr_a * inv(fr_a) == one(BN254Fr)
end

@testset "BN254 Montgomery oracle parity" begin
    rng = MersenneTwister(0xC0D3)

    function check_field_parity(::Type{F}, samples::Int) where {F<:Union{BN254Fq,BN254Fr}}
        p = prime(F)
        for _ in 1:samples
            a_raw = rand(rng, BigInt(0):p - 1)
            b_raw = rand(rng, BigInt(0):p - 1)
            a = F(a_raw, true)
            b = F(b_raw, true)

            @test convert(BigInt, a) == a_raw
            @test convert(BigInt, b) == b_raw

            @test convert(BigInt, a + b) == mod(a_raw + b_raw, p)
            @test convert(BigInt, a - b) == mod(a_raw - b_raw, p)
            @test convert(BigInt, a * b) == mod(a_raw * b_raw, p)
            @test convert(BigInt, a^5) == powermod(a_raw, 5, p)

            if !iszero(a)
                @test convert(BigInt, inv(a)) == invmod(a_raw, p)
            end
            if !iszero(b)
                @test convert(BigInt, a / b) == mod(a_raw * invmod(b_raw, p), p)
            end
        end
    end

    check_field_parity(BN254Fq, 64)
    check_field_parity(BN254Fr, 64)
end
