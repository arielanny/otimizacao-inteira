#Modelo escrito utilizando o JuMP, que é uma linguagem de modelagem específica para #otimização matemática incorporada a linguagem de programação Julia.

# Importa pacotes de otimização
using JuMP, Gurobi, CPLEX, Cbc
using DelimitedFiles

# Comprimento, largura e altura da caixa i (li,wi,hi)
c = [10  5  5;
    10  5  5;
    5  10  5;
    5  10  5;
    10  5  10;
    10  5  10]

n = size(c, 1) # Número de caixas
M = 10^3

# Volume da caixa i (vi)
v = prod.(eachrow(c))

# Origem do conteiner C0=(X0,Y0,Z0) e tamanho C=(X,Y,Z)
C0 = (20,20,20)
C = (12, 10, 15)

# Prepara o modelo
m = Model(Cbc.Optimizer)

# Declara as variáveis
@variable(m, x[i=1:n] >= 0)
@variable(m, y[i=1:n] >= 0)
@variable(m, z[i=1:n] >= 0)
@variable(m, u1[i=1:n, j=1:n], Bin)
@variable(m, u2[i=1:n, j=1:n], Bin)
@variable(m, u3[i=1:n, j=1:n], Bin)
@variable(m, p[i=1:n], Bin)

# Estabelece a função objetivo
@objective(m, Max, sum(v[i]*p[i] for i=1:n))

# Estabelece as restrições
@constraint(m, constraint_01[i=1:n], x[i] >= C0[1]*p[i])
@constraint(m, constraint_02[i=1:n], y[i] >= C0[2]*p[i])
@constraint(m, constraint_03[i=1:n], z[i] >= C0[3]*p[i])
@constraint(m, constraint_04[i=1:n], x[i] <= (C0[1]+C[1]) - c[i,1])
@constraint(m, constraint_05[i=1:n], y[i] <= (C0[2]+C[2]) - c[i,2])
@constraint(m, constraint_06[i=1:n], z[i] <= (C0[3]+C[3]) - c[i,3])
@constraint(m, constraint_07[i=1:n-1, j=i+1:n], x[j] + c[j,1] - x[i] <= M*(u2[i,j] + u3[i,j]))
@constraint(m, constraint_08[i=1:n-1, j=i+1:n], x[i] + c[i,1] - x[j] <= M*(u1[i,j] + u3[i,j]))
@constraint(m, constraint_09[i=1:n-1, j=i+1:n], y[j] + c[j,2] - y[i] <= M*(u1[i,j] + u2[i,j]))
@constraint(m, constraint_10[i=1:n-1, j=i+1:n], y[i] + c[i,2] - y[j] <= M*(2 - (u1[i,j] + u2[i,j])))
@constraint(m, constraint_11[i=1:n-1, j=i+1:n], z[j] + c[j,3] - z[i] <= M*(2 - (u2[i,j] + u3[i,j])))
@constraint(m, constraint_12[i=1:n-1, j=i+1:n], z[i] + c[i,3] - z[j] <= M*(2 - (u1[i,j] + u3[i,j])))
@constraint(m, constraint_13[i=1:n-1, j=i+1:n], 1 <= u1[i,j] + u2[i,j] + u3[i,j] <= 2)

# Resolve o problema
optimize!(m)