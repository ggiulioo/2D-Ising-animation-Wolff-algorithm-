import Base.show
struct Lattice{D} 
    dims::NTuple{D,Int}
    site::UnitRange{Int}
    neig::Vector{Vector{Int}}
end

(Lattice(dims::NTuple{D,Int}) where D) = Lattice(dims,1:prod(dims),neighbors(dims))
Lattice(dims...) = Lattice(dims)

function show(io::IO, x::Lattice{D}) where D
    print(io,"$D-D lattice dims = $(x.dims[1])")
    for i in 2:length(x.dims)
        print(io,"×",x.dims[i])
    end
    print(io," PBC")
end


neighbors(I...) = neighbors(I)

function neighbors(I::NTuple{N,Int}) where N
    neigh = Vector{Vector{Int}}(undef,prod(I))
    vscra = zeros(Int,N)
    pm = (-1,1)
    site = 0
    @inbounds for i in CartesianIndices(I)
        cnt = 0
        site += 1
        neigh[site] = Vector{Int}(undef,2N)
        for pm1 in pm
            for k1=1:N
                cnt += 1
                for k2=1:N                                
                    vscra[k2] = mod1(i.I[k2] + pm1 * (k1 == k2), I[k2])
                end
                neigh[site][cnt] = mysub2ind(I,vscra)
            end
        end
    end   
    neigh
end

function mysub2ind(I::NTuple{N,Int},v::Vector{Int}) where N
    num = v[1]
    pdim  = 1
    @inbounds for i=1:N-1
        pdim *= I[i]
        num += pdim * (v[i+1]-1)
    end
    num
end

nothing
