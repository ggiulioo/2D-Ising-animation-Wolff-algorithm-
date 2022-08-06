include("ising.jl")

function wolffmeasures!(M::Measures,x::Ising,count::Int)
    M.ene[count]=energy(x)
    M.mag[count]=abs(mean(x.spin))
    push!(M.conf,BitArray(map(spin2bool,x.spin)))
    M.time[count]=length(M.conf)
end

function region(x::Ising,S::Array,F::BitArray{1},J::Array{T,2};site::Int=rand(1:x.N)) where T<:Real
    empty!(S)
    push!(S,site)
    F[S[1]]=true
    for i ∈ S
        for j ∈ x.Λ.neig[i]
            p=1-exp(-2*x.β*abs(J[i,j]))
            if rand()<=p && !F[j] && x.spin[j]*x.spin[i]*J[i,j]>=0
                F[j]=true
                push!(S,j)
            end
        end
    end
    F[S].=false
end

function wolff(I::NTuple{N,Int},β::T;  # linear size of the lattice + temperature
                 h::Vector{T}=zeros(T,prod(I)), # external field
                 x₀::Vector{Int}=rand([-1,1],prod(I)), # initial configuration
                 J::Array{T,2}=ones(Float64,prod(I),prod(I)), # couplig matrix
                 nterm::Int=0,  # number of steps for thermalization
                 nsteps::Int=1000,  # number of measures
                 nsweep::Int=0, # number of steps between measurements
                 verbose::Bool=false) where {T<:Real,N}
    x=Ising(I,h,β,x₀)
    M=Measures(zeros(Int,nsteps),β,zeros(T,nsteps),zeros(T,nsteps),Vector{BitArray{1}}())
    S=[]
    F=BitArray{1}(zeros(Bool,x.N))
    for t ∈ 1:nterm
        region(x,S,F,J)
        x.spin[S]=-x.spin[S]
    end
    for t ∈ 1:nsteps
        for n ∈ 1:nsweep
            region(x,S,F,J)
            x.spin[S]=-x.spin[S]
        end
        region(x,S,F,J)
        x.spin[S]=-x.spin[S]
        wolffmeasures!(M,x,t)
    end
    return M
end

nothing
