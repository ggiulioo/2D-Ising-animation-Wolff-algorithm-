using Statistics, Random

include("lattice_p.jl")

mutable struct Ising{T<:AbstractFloat,D}
    N::Int               # number of spin
    β::T                 # inverse temperature
    Λ::Lattice{D}        # Lattice
    spin::Vector{Int}    # spin config
    h::Vector{T}         # external field
    H::Vector{T}         # total field
end

mutable struct Measures{T<:AbstractFloat}
    time::Vector{Int}
    β::T
    ene::Vector{T}
    mag::Vector{T}
    conf::Vector{BitArray{1}}
end

function totalfield(Λ::Lattice,h::Vector{Float64},spin::Vector{Int})
    H=zeros(Float64,prod(I))
    for i ∈ Λ.site
        H[i]=sum(spin[Λ.neig[i]])+h[i]
    end
    return H
end

function Ising(I::NTuple{D,Int},h::Vector{T},β::T,spin::Vector{Int})where T<:AbstractFloat where D
    N=prod(I)
    Λ=Lattice(I)
    H=totalfield(Λ,h,spin)
    return Ising(N,β,Λ,spin,h,H)
end

function energy(x::Ising)
    E=0.0
    for i ∈ x.N
        E-=0.5*x.spin[i]*(x.H[i]+x.h[i])
    end
    return E
end

function onemcstep!(x::Ising,site::Int)
    ΔE=2*x.H[site]*x.spin[site]
    r=rand()
    if exp(-x.β*ΔE)>r
        x.spin[site]=-x.spin[site]
        for i ∈ x.Λ.neig[site]
            x.H[i]=sum(x.spin[x.Λ.neig[i]])+x.h[site]
        end
    end
end

function onemcsweep!(x::Ising)
    perm=randperm(x.N)
    for i ∈ perm
        onemcstep!(x,i)
    end
end

spin2bool=x->x>0 #Funzione astratta, spin2bool(x)=true, se x>0, false, altrimenti

function measures!(M::Measures,x::Ising,count::Int)
    M.ene[count]=energy(x)
    M.mag[count]=mean(x.spin)
    push!(M.conf,BitArray(map(spin2bool,x.spin)))
    M.time[count]=length(M.conf)
end


function mcising(I::NTuple{N,Int},β::T;  # linear size of the lattice + temperature
                 h::Vector{T}=zeros(T,prod(I)), # external field
                 x₀::Vector{Int}=rand([-1,1],prod(I)), # initial configuration
                 nterm::Int=100,  # number of mcsweeps for thermalization
                 nmeas::Int=1000, # number of measurments
                 nsweep::Int=100, # number of mcsweeps between measurements
                 verbose::Bool=false) where {T<:Real,N}
    x=Ising(I,h,β,x₀)
    time=0
    #Termalizziamo
    for t ∈ 1:nterm
        onemcsweep!(x)
    end
    #Adesso misuriamo
    M=Measures(zeros(Int,nmeas),β,zeros(T,nmeas),zeros(T,nmeas),[BitArray{1}([])])
    pop!(M.conf)
    for n ∈ 1:nmeas
        for t ∈ 1:nsweep
            onemcsweep!(x)
        end
        measures!(M,x,n)
    end
    return M
end

function blockmeans(data::Vector{Measures{T}}; field::Symbol=:mag) where T<:AbstractFloat
    nblocks=Int(floor(log(2,length(data[1].time))))
    t,μ,σ=zeros(Float64,nblocks),zeros(Float64,nblocks),zeros(Float64,nblocks)
    n,m=Int(ceil(length(data[1].time)/2)+1),length(data[1].time)
    for i ∈ 1:nblocks
        blocks=[abs.(getfield(data[i],field)[n:m]) for i ∈ 1:length(data)]
        t[i]=(n-1+m) >> 1
        μ[i]=mean(mean.(blocks))
        σ[i]=std(mean.(blocks))/sqrt(length(n:m)-1)
        m=n-1
        n=Int(ceil(m/2))
    end
    return t,μ,σ
end

#Autocorrelation

function autocorrelation(data::Vector{Measures{T}},Δt::Int; field::Symbol=:mag) where T<:AbstractFloat
    ac=[acsingle(getfield(m,field),Δt) for m ∈ data]
    μ=mean(ac)
    σ=std(ac)/sqrt(length(data))
    return μ,σ
end

function acsingle(data::Vector{T},Δt::Int) where T
    c₀=acupdate(data,0)
    ac=zeros(typeof(c₀), Δt+1)
    ac[1]=c₀
    for t ∈ 1:Δt
        ac[t+1]=acupdate(data,t)
    end
    return ac./c₀
end

function acupdate(data::Vector{T}, Δt::Int) where T<:AbstractFloat
    a,b,c,cnt=zero(T),zero(T),zero(T),0
    for t in 1:length(data)-Δt
        a+=data[t]*data[t+Δt]
        b+=data[t]
        c+=data[t+Δt]
        cnt+=1
    end
    a,b,c=a/cnt,b/cnt,c/cnt
    return a-b*c
end

function acupdate(data::Vector{BitArray{1}}, Δt::Int)
    a,b,c,cnt=0,0,0,0
    N=length(data[1])
    for t in 1:length(data)-Δt
        a+=2*sum(data[t].==data[t+Δt]) - N
        b+=2*sum(data[t])-N
        c+=2*sum(data[t+Δt])-N
        cnt+=1
    end
    a,b,c=a/(N*cnt),b/(N*cnt),c/(N*cnt)
    return a-b*c
end


nothing
