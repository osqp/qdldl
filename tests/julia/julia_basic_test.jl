# Unit test for QDLDL.jl

using QDLDL, Test, Random, Statistics, LinearAlgebra, SparseArrays

rng = Random.MersenneTwister(131123)


function generateSparsePosDefMatrix(rng,n::Int64,density::Float64)
      X = sprand(rng,n,n,density)
      X = 0.5(X+X')
      # make random  symmetric matrix pos def
      X = X+(2*n+1)*sparse(1.0I,n,n);
      return X
  end

dim = rand(rng,50:1:200);
density=rand(rng,0.1:0.1:0.4)
A = generateSparsePosDefMatrix(rng,dim,density)
xtrue = rand(rng,dim)
btrue = A*xtrue


F = QDLDL.qdldl(A)

@testset "Basic QDLDL.jl test" begin

  x1 = QDLDL.solve(F,btrue)
  @test norm(xtrue-x1,Inf) < 1e-10

  x2 = F\btrue
  @test norm(xtrue-x2,Inf) < 1e-10

  x3 = copy(btrue)
  QDLDL.solve!(F,x3)
  @test norm(xtrue-x3,Inf) < 1e-10
end
nothing