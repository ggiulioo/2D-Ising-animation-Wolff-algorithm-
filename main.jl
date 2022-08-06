include("graphics_ising.jl")

I = (30,30)            #lattice dimensions
N = prod(I)		     #total spins
x₀ = rand([-1,1],N)    #initial configuration
h = zeros(Float64,N)   #external field
β = 10000.0            #beta, here very high
Ter_steps = 40         #termalisation steps

x = Ising(I,h,β,x₀)    #Ising model definition


for j ∈ 1:Ter_steps    #termalisation
    onemcsweep!(x)
end

plot()
anim= @animate for i ∈ 1:5 #animation
        onemcsweep!(x)
        commandplot!(I,x.spin)
    end
gif(anim, fps = 5)
