# Script para carregar instâncias teste
using CSV
using DataFrames



# lista de endereços
arquivos=["C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-20-L-C-50.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-20-U-C-50.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-20-U-C-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-20-U-R-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-40-D-R-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-40-F-R-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-40-U-R-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-60-D-R-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-60-F-C-90.csv",
"C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-60-F-R-90.csv"]

dados_caixas = CSV.File("C:/Users/mizus/OneDrive/Área de Trabalho/Palete Trabalho/Trabalho Interira/InstanciasSelecionadas/ep3-20-L-C-50.csv", header=0)
df_dados_caixas=DataFrame(dados_caixas)
rename!(df_dados_caixas,["etc","index","comprimento","altura","largura","volume","quantidade"])
D_P=df_dados_caixas[1,["index","altura","comprimento"]]

for arquivo in arquivos
    # Leitura de arquivos
    dados_caixas = CSV.File(arquivo, header=0)
    df_dados_caixas=DataFrame(dados_caixas)
    D_P=df_dados_caixas[1,["index","altura","comprimento"]]
    ary_posicoes=Array{Int64}(undef, 0, 4)
    rename!(df_dados_caixas,["etc","index","comprimento","altura","largura","volume","quantidade"])
    ary_posicoes=df_dados_caixas[2:end, ["index","comprimento","largura","altura"]]
    # Rodar Paletização

end