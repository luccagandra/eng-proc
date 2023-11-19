import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt

matriz =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
         [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
         [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
         [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

matriz_2 = [[3.0, 170.0, 60.0], # WCp, Q1_T0, Q1_Td
            [1.5, 150.0, 30.0], # WCp, Q2_T0, Q2_Td
            [2.0, 20.0, 135.0],  # WCp, F1_T0, F1_Td
            [4.0, 80.0, 140.0]] # WCp, F2_T0, F2_Td

matriz_trabalho = [[3.0, 170.0, 60.0], # WCp, Q1_T0, Q1_Td
            [1.5, 150.0, 30.0], # WCp, Q2_T0, Q2_Td
            [2.0, 30.0, 140.0],  # WCp, F1_T0, F1_Td
            [4.0, 80.0, 140.0]] # WCp, F2_T0, F2_Td

intervals_trabalho = [160.0, 140.0, 80.0, 50.0, 30.0, 20.0]

matriz_acima_pinch = [[3.0, 170.0, 90.0], # WCp, Q1_T0, Q1_Td
            [1.5, 150.0, 90.0], # WCp, Q2_T0, Q2_Td
            [2.0, 80.0, 140.0],  # WCp, F1_T0, F1_Td
            [4.0, 80.0, 140.0]] # WCp, F2_T0, F2_Td

matriz_abaixo_pinch = [[3.0, 90.0, 60.0], # WCp, Q1_T0, Q1_Td
            [1.5, 90.0, 30.0], # WCp, Q2_T0, Q2_Td
            [2.0, 30.0, 80.0]]  # WCp, F1_T0, F1_Td

intervals_2 = [160.0, 140.0, 135.0, 80.0, 50.0, 20.0]

def RPS(matriz, temp_min, first_hot, first_cold):
    count = 1

    while count <= 10:   
        chosen_hot = -1 
        chosen_cold = -1 

        if count != 1:
            hot_candidates = [(matriz[i][1], i) for i in range(2)]
            chosen_hot = 0 if hot_candidates[0][0] > hot_candidates[1][0] else 1
            QMT0 = matriz[chosen_hot].copy()

            max_value_cold = -1000
            
            for i in range(2, 3): # Troca para 3
                if QMT0[1] > (matriz[i][1] + temp_min) and matriz[i][1] > max_value_cold:
                    chosen_cold = i

            FMT0 = matriz[chosen_cold].copy()

            if chosen_cold == -1:
                break
        else:
            chosen_hot = first_hot
            chosen_cold = first_cold 
            QMT0 = matriz[chosen_hot].copy()
            FMT0 = matriz[chosen_cold].copy()

        print("Q", count, "Troca, corrente quente:", chosen_hot, QMT0)
        print("F", count, "Troca, corrente fria:", chosen_cold, FMT0)

        if QMT0[1] - FMT0[2] <= temp_min: # TEQ* - TSF < temp_min
            FMT0[2] = QMT0[1] - temp_min

        if QMT0[2] - FMT0[1] <= temp_min: # TSQ - TEF* < temp_min
            QMT0[2] = FMT0[1] + temp_min
        #print(QMT0[2])
        Q = min((QMT0[1] - QMT0[2])*QMT0[0], (FMT0[2] - FMT0[1])*FMT0[0])
        
        if Q == (QMT0[1] - QMT0[2])*QMT0[0]: # Q = oferta
            QMT0[2] = QMT0[2]
            FMT0[2] = round(FMT0[1] + Q/FMT0[0],1)
        
        if Q == (FMT0[2] - FMT0[1])*FMT0[0]: # Q = demanda
            QMT0[2] = round(QMT0[1] - Q/QMT0[0],1)
            FMT0[2] = FMT0[2]
        
        #print("QMT0: ", QMT0)
        #print("FMT0: ", FMT0)
        
        matriz[chosen_hot][1] = QMT0[2] # Entradas são as saídas do sistema anterior
        matriz[chosen_cold][1] = FMT0[2]
        
        print("RPS:", matriz)
        #display_table(matriz)
        
        count += 1

#RPS(matriz_abaixo_pinch, 10, 0, 2) # Q1xF1

RPS(matriz_acima_pinch, 10, 0, 3) # Q1xF2
#RPS(matriz_acima_pinch, 10, 1, 2) # Q2xF1
#RPS(matriz_acima_pinch, 10, 1, 3) # Q2xF2

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

def is_pair_overlapping(pair1, pair2):
    # Sort the pairs to ensure proper comparison
    pair1 = sorted(pair1)
    pair2 = sorted(pair2)

    # Check if pair1 is completely within pair2
    return pair1[0] >= pair2[0] and pair1[1] <= pair2[1]

def offer_demand(matriz, intervals, temp_min):
    
    Rk = 0
    a1 = 1
    a2 = 1
    offer_demand = []

    for i in range(len(intervals)-1): 
        
        if is_pair_overlapping((intervals[i]+temp_min, intervals[i+1]+temp_min), (matriz[0][2], matriz[0][1])): # passa por arrow1
            a1 = matriz[0][0]
        else:
            a1 = 0

        if is_pair_overlapping((intervals[i]+temp_min, intervals[i+1]+temp_min), (matriz[1][2], matriz[1][1])): # passa por arrow2
            a2 = matriz[1][0]
        else:
            a2 = 0

        if is_pair_overlapping((intervals[i], intervals[i+1]), (matriz[2][2], matriz[2][1])): # passa por arrow3
            a3 = matriz[2][0]
        else:
            a3 = 0

        if is_pair_overlapping((intervals[i], intervals[i+1]), (matriz[3][2], matriz[3][1])): # passa por arrow4
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
    
    #print(matriz)
    #print(offer_demand)
    display_table(offer_demand)
        
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

    print("Intervals", intervals)
    offer_demand(matriz, intervals, 10)

criar_grafico(matriz_trabalho)

# Função chamada quando o botão é clicado
def obter_numeros():
    matriz = []
    for i in range(4):
        linha = []
        for j in range(3):
            valor = float(entrada_matriz[i][j].get())
            linha.append(valor)
        matriz.append(linha)

    #exibir_matriz(matriz)
    criar_grafico(matriz)
    RPS(matriz, 10)

def exibir_matriz(matriz):
    for i in range(4):
        for j in range(3):
            print(matriz[i][j], end="\t")


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

#offer_demand(matriz_2, intervals_2, 10)

#RPS(matriz_2, 10)
