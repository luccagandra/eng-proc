from integration import IntegraçãoEnergética
import pandas as pd
import tkinter as tk
from tkinter import ttk

input = False
matrix = []

# Função chamada quando o botão é clicado
def obter_numeros():
    for i in range(4):
        linha = []
        for j in range(3):
            valor = float(entrada_matriz[i][j].get())
            linha.append(valor)
        matrix.append(linha)
    pass

if input == True:
    # Cria a janela da GUI
    janela = tk.Tk()
    janela.title("Inserir Matriz 4x3")

    # Cria rótulos para as colunas
    coluna_legendas = ["WCp, KW/°C", "To, °C", "Td, °C"]
    for j in range(3):
        coluna_legenda = tk.Label(janela, text=coluna_legendas[j])
        coluna_legenda.grid(row=0, column=j + 1)

    # Cria rótulos para as linhas e entradas para a matriz
    entrada_matriz = [[None] * 3 for _ in range(4)]
    linha_legendas = ["Q1", "Q2", "F1", "F2"]
    for i in range(4):
        linha_legenda = tk.Label(janela, text=linha_legendas[i])
        linha_legenda.grid(row=i + 1, column=0)
        for j in range(3):
            entrada_matriz[i][j] = tk.Entry(janela)
            entrada_matriz[i][j].grid(row=i + 1, column=j + 1)

    # Cria um botão para processar a matriz
    botao_processar = tk.Button(janela, text="Processar", command=obter_numeros)
    botao_processar.grid(row=5, columnspan=4)

    # Inicializa a GUI
    janela.mainloop()
else:
    matrix =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
            [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
            [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
            [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

    matriz_aula = [[3.0, 170.0, 60.0], # WCp, Q1_T0, Q1_Td
                [1.5, 150.0, 30.0], # WCp, Q2_T0, Q2_Td
                [2.0, 30.0, 140.0],  # WCp, F1_T0, F1_Td
                [4.0, 80.0, 140.0]] # WCp, F2_T0, F2_Td

matriz_escolhida = matrix
user_input = True

delta_T_min = 10
loop = IntegraçãoEnergética(matriz_escolhida, delta_T_min, Qx=None, Fx=None, user_input=user_input) #Qx=0 e Fx=1 dão erro
loop.criar_grafico()
loop.offer_demand()

print("\n","Matriz original: \n")
print(pd.DataFrame(matriz_escolhida, index=['']*len(matrix), columns=['']*len(matrix[0])),"\n")

# Quais as combinações possíveis pelo RPS?
loop.loop_RPS("QmTOxFmTO")
loop.completando_utilidades()
loop.plot_multiple()

# Das correntes que sobraram, pode haver outro loop
if loop.is_there_two_chains == True:
    print("\n","------------------------------","\n")

    print(f"A combinação fora do loop é: Q{loop.chains_not_in_the_loop[0]+1}xF{loop.chains_not_in_the_loop[1]+1}")
    matriz_loop2 = loop.matrix
    loop2 = IntegraçãoEnergética(matriz_loop2, delta_T_min, Qx=loop.chains_not_in_the_loop[0], Fx=loop.chains_not_in_the_loop[1], user_input=user_input)
    loop2.loop_RPS("QmTOxFmTO")
    loop2.completando_utilidades()
    loop2.plot_multiple()

    print("\n","------------------------------","\n")
    print(f"Custo cap total (duas redes): {round(loop.custo_cap,2)}+{round(loop2.custo_cap,2)}={round(loop.custo_cap+loop2.custo_cap,2)}")
    print(f"Custo útil total (duas redes): {round(loop.custo_util,2)}+{round(loop2.custo_util,2)}={round(loop.custo_util+loop2.custo_util,2)}")
    print(f"Custo total (duas redes): {round(loop.custo_util+loop2.custo_util+loop.custo_cap+loop2.custo_cap,2)}")



"""
QMTOxFMTO
QmTOxFmTO
"""