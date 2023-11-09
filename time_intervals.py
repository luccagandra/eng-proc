import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt

def display_table(data):
    root = tk.Tk()
    root.title("Table Display")

    columns = ["Intervalo"] + ["R(k-1)"] + ["Oferta"] + ["Demanda"] + ["Sk"]

    # Create a Treeview widget
    table = ttk.Treeview(root, columns=columns, show="headings")

    # Set column headings with center alignment
    for col in columns:
        table.heading(col, text=col, anchor="center")

    # Insert data into the table with center alignment
    for i, row in enumerate(data, start=1):
        values = [i] + row
        table.insert("", "end", values=values, tags=("centered",))

    # Center the cell values
    table.tag_configure("centered", anchor="center")

    # Pack the Treeview widget
    table.pack()

    root.mainloop()

def offer_demand(matriz, intervals):
    
    Rk = 0
    a1 = 1
    a2 = 1
    offer_demand = []

    for i in range(len(intervals)-1): 
        
        if matriz[0][2] <= (intervals[i] + intervals[i+1])/2 <= matriz[0][1]: # passa por arrow1
            a1 = matriz[0][0]
        else:
            a1 = 0

        if matriz[1][2] <= (intervals[i] + intervals[i+1])/2 <= matriz[1][1]: # passa por arrow2
            a2 = matriz[1][0]
        else:
            a2 = 0

        if matriz[2][1] <= (intervals[i] + intervals[i+1])/2 <= matriz[2][2]: # passa por arrow3
            a3 = matriz[2][0]
        else:
            a3 = 0

        if matriz[3][1] <= (intervals[i] + intervals[i+1])/2 <= matriz[3][2]: # passa por arrow4
            a4 = matriz[3][0]
        else:
            a4 = 0

        y = (intervals[i]-intervals[i+1])
        lines = [Rk, y*a1 +y*a2, y*a3 +y*a4, Rk + y*a1 +y*a2 - y*a3 - y*a4]

        if lines[3] < 0: # pinch
            Rk = 0 
        else:
            Rk = lines[3]
        
        offer_demand.append(lines)
    
    display_table(offer_demand)
    #print(offer_demand)
        
#offer_demand(matriz, intervals)

def draw_arrow(text,x, start, end, color_arrow):
    arrowprops = dict(arrowstyle='->', linestyle='-', linewidth=2, color=color_arrow)
    plt.annotate(text, xy=(x, end), xytext=(x, start), arrowprops=arrowprops, color=color_arrow)

def criar_grafico(matriz):
    y1 = [matriz[0][1], matriz[0][2], matriz[1][1], matriz[1][2]]  # Q1_To, Q1_Td, Q2_To, Q2_Td
    y2 = [matriz[2][1], matriz[2][2], matriz[3][1], matriz[3][2]]  # F1_To, F1_Td, F2_To, F2_Td

    x1 = [0, 0.5, 0.5, 1]  # Pontos de mudança nos degraus para y1
    x2 = [0.5, 1, 0, 0.5]  # Pontos de mudança nos degraus para y2
    
    y_ticks = []
    intervals = []

    for i in y1:
        plot_y1_step = [i,i,i-10, i-10]
        plt.step(x1, plot_y1_step)

        y_ticks.append(i)
        y_ticks.append(i-10)
        intervals.append(i-10)

    for i in y2:
        plot_y2_step = [i+10,i+10,i, i]
        plt.step(x1, plot_y2_step)
        y_ticks.append(i)
        y_ticks.append(i+10)
        intervals.append(i)

    intervals = list(set(intervals)) # sort and remove duplicates
    intervals.sort(reverse=True)
    #print(intervals)
    plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
    plt.xticks([])

    draw_arrow("",0.2,matriz[0][1],matriz[0][2], "red")
    draw_arrow("",0.4,matriz[1][1],matriz[1][2], "red")

    draw_arrow("",0.6,matriz[2][1],matriz[2][2], "blue")
    draw_arrow("",0.8,matriz[3][1],matriz[3][2], "blue")

    #plt.xlabel('Eixo X')
    #plt.ylabel('Eixo Y')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.show()

    offer_demand(matriz, intervals)

#criar_grafico(matriz)

# Função chamada quando o botão é clicado
def obter_numeros():
    matriz = []
    for i in range(4):
        linha = []
        for j in range(3):
            valor = float(entrada_matriz[i][j].get())
            linha.append(valor)
        matriz.append(linha)
    exibir_matriz(matriz)
    criar_grafico(matriz)

def exibir_matriz(matriz):
    for i in range(4):
        for j in range(3):
            print(matriz[i][j], end="\t")
        print()

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

