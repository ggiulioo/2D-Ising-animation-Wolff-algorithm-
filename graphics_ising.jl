using Plots
using LinearAlgebra
include("ising.jl")
include("wolff.jl")

function pos2indice(tupla, posizionex, posizioney)
    indice=(posizionex-1)*tupla[2]+posizioney
    return indice
end

"""
Help della funzione latticeplot(tupla, x).
Latticeplot esegue il grafico di una configurazione di spin x sul reticolo individuato dalla tupla.
"""
function latticeplot(tupla, x)
    s=plot(title="Regione di Spin, -1 rosso e 1 blu")
    for i ∈ (1:tupla[1]) #indice colonne x
        for j ∈ (1:tupla[2]) #indice righe y
            if x[pos2indice(tupla,i,j)]==1
                scatter!([i],[j], markercolor=:blue, markersize=7, leg=false)
            elseif x[pos2indice(tupla,i,j)]==-1
                scatter!([i],[j], markercolor=:red, markersize=7,leg=false)
            end
        end
    end
return display(s)
end
function commandplot!(tupla, x)
    for i ∈ (1:tupla[1]) #indice colonne x
        for j ∈ (1:tupla[2]) #indice righe y
            if x[pos2indice(tupla,i,j)]==1
                scatter!([i],[j], markercolor=:blue, markersize=7, leg=false)
            elseif x[pos2indice(tupla,i,j)]==-1
                scatter!([i],[j], markercolor=:red, markersize=7,leg=false)
            end
        end
    end
end
 nothing
