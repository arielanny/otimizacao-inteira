
## Gerar conjunto O (length)
coord_O = Int[] # Conjunto dos pontos p horizontais (length)
# n_max onde as iteracoes devem parar,
# o grid nao pode ser maior que isso pois a caixa nao cabe
# n_max_H = Comprimento(L)/(Lado da Caixa(i))
n_max_O = Vector{Int64}[]
# Cria um n_max_O para cada caixa
for n_i in 1:length(ary_posicoes[:,1])
    push!(n_max_O, Vector{Int64}(0:trunc(Int,D_P[1]/ary_posicoes[n_i,2])))
end
println(n_max_O)
ary_comb_n=collect(Iterators.product(n_max_O...))

ary_comb_n=collect(Iterators.product([[0, 1, 2, 3],[0, 1, 2, 3]]...))

# Posicoes que as caixas odem ocupar
# restricao para criacao do grid D_P(L) - max(d_c)
Max_O = D_P[1] - minimum(ary_posicoes[:,2])
# loop para montar o grid
for n_i in 1:length(ary_comb_n)
    coord_p=sum(collect(ary_comb_n[n_i]) .* ary_posicoes[:,2]) # produto vetorial
    if coord_p<=Max_O
        append!(coord_O,coord_p)
    end
end
unique!(coord_O)
println("Conjunto O de posições p (bottom-left) que pertmitem caixas, O=",coord_O)
# O=[0, 29, 43, 27, 54, 26, 53, 52] normal






n_max_O

for conjunto in n_max_O
    println(conjunto)
    for sub in conjunto
        println(sub)
    end
end