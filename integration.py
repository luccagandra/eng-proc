import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

def plot_single_transfer(x, y, number, values,type, offset=0.1, single=True, first=False):
    """
    Plot a point in a 2D Cartesian space, display a straight line segment around it,
    plot a number slightly offset from the center, and add labels at the ends of the segments.

    Parameters:
    - x: x-coordinate of the point
    - y: y-coordinate of the point
    - number: The number to be plotted slightly offset from the center
    - offset_x: Offset in the x-direction for the number (default is 0.1)
    - offset_y: Offset in the y-direction for the number (default is 0.1)
    - values: Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
    """

    global custo_cap
    global custo_util
    global custo_total

    if first==True:
        custo_cap = 0
        custo_util = 0
        custo_total = 0

    # Plotting the point
    plt.scatter(x, y, color='black', label='Point')

    # Displaying a straight line segment around the point
    segment_length = 1  # Length of the segment in both x and y directions

    if type == 'reta': # Usando trocador de integração
        plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
        plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')

        if values[2] == None:
            text_hot = ""
        else:
            text_hot = f'Q{values[0]} = {values[2]}'

        if values[3] == None:
            text_cold = ""
        else:
            text_cold = f'F{values[1]} = {values[3]}'

        # Adding labels 
        plt.text(x - segment_length/2, y + offset, text_hot, color='red', ha='right', va='center', fontsize=10)
        plt.text(x + offset, y + segment_length/2, text_cold, color='blue', ha='center', va='bottom', fontsize=10)

        plt.text(x + segment_length/2, y + offset, f'{values[4]}', color='red', ha='left', va='center', fontsize=10)
        plt.text(x + offset, y - segment_length/2, f'{values[5]}', color='blue', ha='center', va='top', fontsize=10)

        # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, WCp_Q, WCp_F
        # Cálculo da Oferta e Demanda
        oferta =  values[6] * (values[2] - values[4])
        demanda = values[7] * (values[5] - values[3])

        Q = min(oferta, demanda)

        U = 0.75

        custo_cap += custo_do_trocador(values, 'logarítmico', Q, U)

    elif type == 'hot': # Usando água (resfriador)
        plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='blue')
        plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
        
        plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, values[5], color='blue', ha='left', va='center', fontsize=10)
        plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, values[3], color='blue', ha='right', va='center', fontsize=10)
        plt.text(x + segment_length/2, y + offset, values[4], color='red', ha='left', va='center', fontsize=10)

        Q =  values[6] * (values[2] - values[4])

        U = 0.75

        custo_cap += custo_do_trocador(values, 'logarítmico', Q, U)

        # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor

        vazao = calculo_da_vazao(Q, 'água')
        custo_util += custo_utilidades(vazao, 'água')

    elif type == 'cold': # Usando vapor (aquecedor)
        plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='red')
        plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')

        plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, values[2], color='red', ha='left', va='center', fontsize=10)
        plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, values[4], color='red', ha='right', va='center', fontsize=10)
        plt.text(x + offset, y - segment_length/2, values[5], color='blue', ha='center', va='top', fontsize=10)

        Q = values[7] * (values[5] - values[3])
        U = 1

        custo_cap += custo_do_trocador(values, 'logarítmico', Q, U)

        # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor

        vazao = calculo_da_vazao(Q, 'vapor')
        custo_util += custo_utilidades(vazao, 'vapor')

    # Plotting the number slightly offset from the center
    plt.text(x + offset/2, y + offset/2, str(number), color='black', ha='center', va='center', fontsize=12, fontweight='bold')

    custo_cap = round(custo_cap,1)
    custo_util = round(custo_util,1)
    custo_total = round(custo_cap + custo_util,1)
    fig_size = plt.gcf().get_size_inches()

    if single == True:
        # Adding labels, title, legend, and grid for better visualization
        plt.title('Rede Final')
        plt.grid(False)
        plt.axis('off')

        text1 = f"Ccap = {custo_cap}"
        text2 = f"Cutil = {custo_util}"
        text3 = f"Ctotal = {custo_total}"

        plt.text(x, 1.2, text1, ha='right', va='top', color='red', fontsize=10)
        plt.text(x, 1, text2, ha='right', va='top', color='green', fontsize=10)
        plt.text(x, 0.8, text3, ha='right', va='top', color='blue', fontsize=10)

        #print(custo_cap, custo_util, custo_total)
        plt.show()

def calculo_da_vazao(Q, tipo):
    if tipo == 'água':
        W = Q/((0.00116)*(50-30))
    if tipo == 'vapor':
        W = Q/0.48

    return W

def custo_do_trocador(values, tipo, Q, U):

    # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F

    delta_1 = values[2] - values[5] # TEQ - TFS
    delta_2 = values[4] - values[3] # TSQ - TEF

    if tipo == 'aritmético':
        area = Q / (U*(delta_1+delta_2)/2)
    
    if tipo == 'logarítmico':
        if delta_1 == delta_2: # Não podemos passar na fórmula abaixo pq gera uma indeterminação
            area = delta_1
        else:
            area = Q / (U*(delta_1-delta_2)/np.log(delta_1/delta_2))

    custo_do_trocador = 130 * (math.pow(area, (65/100)))

    return custo_do_trocador

def custo_utilidades(consumo, tipo):
    if tipo == 'água':
        custo_unitario = 0.00005 # $/kg
    if tipo == 'vapor':
        custo_unitario = 0.0015 # $/kg

    return 8500*(custo_unitario*consumo)

def plot_multiple_transfers(plot_list, last_matrix):
    # Plotting the point
    plt.scatter(1, 1, color='red', label='Point')

    x,y = 0,0
    plotted_out = []
    values_to_clear = []

    for i in range(len(last_matrix)):
        if last_matrix[i][1] == last_matrix[i][2]:
            values_to_clear.append(last_matrix[i][1])

    for i in range(len(plot_list)):

        values = plot_list[i]

        if i == 0:
            plotted_out.append([values[4], 0,0, 'hot']) # TSQ, plotted at 0,0
            plotted_out.append([values[5], 0,0, 'cold']) # TSF, plotted at 0,0

            if len(plot_list)-1 == 0:
                plot_single_transfer(0,0,(i+1),values,'reta', single=False, first=True)

                for a in range(len(plotted_out)): # Completando com utilidades
                    
                    if a+1 == len(plotted_out):
                        single = True
                    else:
                        single = False  
                    
                    if plotted_out[a][3] == 'hot': # água de 30 a 50
                        
                        for k in range(len(last_matrix)):
                            if last_matrix[k][1] == plotted_out[a][0]:
                                meta = last_matrix[k][2]
                                WCp_Q = last_matrix[k][0]

                        values = [0,0,plotted_out[a][0],30,meta,50,WCp_Q,0] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F
                        plot_single_transfer(plotted_out[a][1]+1, plotted_out[a][2], (a+i+2),values,'hot', single=single)

                    if plotted_out[a][3] == 'cold': # vapor de 250 a 250

                        for h in range(len(last_matrix)):
                            if last_matrix[h][1] == plotted_out[a][0]:
                                meta = last_matrix[h][2]
                                WCp_F = last_matrix[h][0]

                        values = [0,0,250,plotted_out[a][0],250,meta, 0, WCp_F] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
                        plot_single_transfer(plotted_out[a][1], plotted_out[a][2]-1, (a+i+2),values,'cold', single=single)

                # Quando uma troca n pode ser realizada na rede, ela é realizada individualmente.

            plot_array = np.array(plot_list)

            if not 1 in plot_array[:, [0]]: # Não tem Q1 na rede, a corrente deve ser resfriada
                values = [0,0,last_matrix[0][1],30,last_matrix[0][2],50,last_matrix[0][0],0] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F
                plot_single_transfer(0, 0, 1,values,'hot', single=True)            

            if not 2 in plot_array[:, [0]]: # Não tem Q2 na rede, a corrente deve ser resfriada
                values = [0,0,last_matrix[1][1],30,last_matrix[1][2],50,last_matrix[1][0],0] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F
                plot_single_transfer(0, 0, 1,values,'hot', single=True)
    
            if not 1 in plot_array[:, [1]]: # Não tem F1 na rede, a corrente deve ser resfriada
                values = [0,0,250,last_matrix[2][1],250,last_matrix[2][2],0,last_matrix[2][0]] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F
                plot_single_transfer(0, 0, 1,values,'cold', single=True)

            if not 2 in plot_array[:, [1]]: # Não tem F2 na rede, a corrente deve ser resfriada
                values = [0,0,250,last_matrix[3][1],250,last_matrix[3][2],0,last_matrix[3][0]] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor, Wcp_Q, Wcp_F
                plot_single_transfer(0, 0, 1,values,'cold', single=True)

            else:
                plot_single_transfer(0,0,(i+1),values,'reta', single=False, first=True)
            prev_value_x = 0
            prev_value_y = 0
            prev_values = values

        else:
            if values[3] == prev_values[5]: # Conexão por corrente fria
                x = prev_value_x
                y = prev_value_y - 1
                values_to_clear.append(values[3])

            if values[2] == prev_values[4]: # Conexão por corrente quente
                x = prev_value_x + 1
                y = prev_value_y
                values_to_clear.append(values[2])

            prev_value_x = x
            prev_value_y = y
            prev_values = values

            if i == len(plot_list)-1: # Última iteração
                plotted_out.append([values[4], x,y, 'hot']) 
                plotted_out.append([values[5], x,y, 'cold']) 

                plot_single_transfer(x,y,(i+1),values,'reta', single=False) # Terão elementos livres!

                #print(plotted_out)
                #print(values_to_clear)

                for j in range(len(plotted_out)-1, -1, -1):  
                    if plotted_out[j][0] in values_to_clear:
                        plotted_out.pop(j)
                
                for a in range(len(plotted_out)): # Completando com utilidades
                    
                    if a+1 == len(plotted_out):
                        single = True
                    else:
                        single = False  
                    
                    if plotted_out[a][3] == 'hot': # água de 30 a 50
                        
                        for k in range(len(last_matrix)):
                            if last_matrix[k][1] == plotted_out[a][0]:
                                meta = last_matrix[k][2]
                                WCp_Q = last_matrix[k][0]

                        values = [0,0,plotted_out[a][0],30,meta,50,WCp_Q,0] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
                        plot_single_transfer(plotted_out[a][1]+1, plotted_out[a][2], (a+i+2),values,'hot', single=single)

                    if plotted_out[a][3] == 'cold': # vapor de 250 a 250

                        for h in range(len(last_matrix)):
                            if last_matrix[h][1] == plotted_out[a][0]:
                                meta = last_matrix[h][2]
                                WCp_F = last_matrix[h][0]

                        values = [0,0,250,plotted_out[a][0],250,meta, 0, WCp_F] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
                        plot_single_transfer(plotted_out[a][1], plotted_out[a][2]-1, (a+i+2),values,'cold', single=single)

            else:
                plotted_out.append([values[4], x,y, 'hot']) 
                plotted_out.append([values[5], x,y, 'cold']) 
                plot_single_transfer(x,y,(i+1),values,'reta', single=False)
            

def perform_transform(matrix, Qx, Fx, delta_T_min):

    QMT0 = matrix[Qx-1].copy()
    FMT0 = matrix[Fx+1].copy()

    """
        TEQ*=ToQ
           ↓
TEF*=ToF  →o→ TSF=TdF
           ↓
        TSQ=TdQ
    """

    # WCp
    WCp_Q = QMT0[0]
    WCp_F = FMT0[0]

    # Temperaturas de entrada
    TEQ_fixado = QMT0[1] # TEQ* = ToQ
    TEF_fixado = FMT0[1] # TEF* = ToF

    # Temperaturas de saída
    TSQ_meta = QMT0[2] # TSQ = TdQ
    TSF_meta = FMT0[2] # TSF = TdF

    # Se entrada quente - saída fria é menor que ΔT(min)
    if TEQ_fixado - TSF_meta < delta_T_min:
        TSF_meta = TEQ_fixado - 10

    # Se saída quente - entrada fria é menor que ΔT(min)
    if TSQ_meta - TEF_fixado < delta_T_min:
        TSQ_meta = TEF_fixado + 10
    
    # Cálculo da Oferta e Demanda
    oferta =  WCp_Q * (TEQ_fixado - TSQ_meta)
    demanda = WCp_F * (TSF_meta - TEF_fixado)

    # Q

    Q = min(oferta, demanda)

    # Se Q = oferta, TSF = entrada fria + Q / WCp_F
    if Q == oferta:
        TSQ_meta = TSQ_meta # Temperatura de saída quente confirmada
        TSF_meta = TEF_fixado + round(Q/WCp_F,1) 
    
    # Se Q = demanda, TSQ = entrada quente - Q / WCp_F
    if Q == demanda:
        TSF_meta = TSF_meta
        TSQ_meta = TEQ_fixado - round(Q/WCp_Q,1) 

    # TSF_meta e TSQ_meta são atualizados
    QMT0[1] = TSQ_meta
    FMT0[1] = TSF_meta

    new_matrix = matrix.copy()

    new_matrix[Qx-1] = QMT0
    new_matrix[Fx+1] = FMT0

    return new_matrix

def combinations(matrix):
    # Temperaturas de entrada das corrents quentes 
    # ToQ[0] = Q1, ToQ[1] = Q2
    ToQ = (matrix[0][1], matrix[1][1])

    # Temperaturas de entrada das corrents frias
    # ToF[0] = F1, ToF[1] = F2   
    ToF = (matrix[2][1], matrix[3][1])

    index_max_ToQ = np.argmax(ToQ)
    index_min_ToQ = np.argmin(ToQ)
    index_max_ToF = np.argmax(ToF)
    index_min_ToF = np.argmin(ToF)

    possible_combinations = []
    
    max_combination = (index_max_ToQ,index_max_ToF)
    min_combination = (index_min_ToQ,index_min_ToF)
    mix_1 = (index_min_ToQ,index_max_ToF)
    mix_2 = (index_max_ToQ,index_min_ToF)
    #print("Matrix: ", matrix)
    #print("ToQ: ", ToQ)
    #print("ToF: ", ToF)

    # Combinações válidas para Q1
    for i in range(len(ToQ)):
        if (ToQ[0] > ToF[i]) and (matrix[0][1] != matrix[0][2]) and (matrix[2+i][1] != matrix[2+i][2]):
            possible_combinations.append((0,i))

    # Combinações válidas para Q2
    for i in range(len(ToQ)):
        if (ToQ[1] > ToF[i]) and (matrix[1][1] != matrix[1][2]) and (matrix[2+i][1] != matrix[2+i][2]):
            possible_combinations.append((1,i))

    #print(possible_combinations)
    
    QMTOxFMTO = max_combination if max_combination in possible_combinations else None
    QmTOxFmTO = min_combination if min_combination in possible_combinations else None
    mix_1 = mix_1 if mix_1 in possible_combinations else None
    mix_2 = mix_2 if mix_2 in possible_combinations else None
    
    combinations_ranked = QMTOxFMTO, QmTOxFmTO, mix_1, mix_2
    
    """
    print("Combinações possíveis: ", possible_combinations)
    print("QMTOxFMTO:","Q", max_combination[0], "F", max_combination[1])
    print("QmTOxFmTO:","Q", min_combination[0], "F", min_combination[1])
    """
    
    QMTOxFMTO = (max_combination[0],max_combination[1])
    QmTOxFmTO = (min_combination[0],min_combination[1])

    if len(possible_combinations) == 0:
        QMTOxFMTO = None
        QmTOxFmTO = None
    
    return combinations_ranked # QMTOxFMTO, QmTOxFmTO, QmTOxFMTO, QMTOxFmTO


def perform_RPS(matrix, Qx=None, Fx=None, delta_T_min=10, plot=False):

    comb = (0,0)
    count = 1
    plot_multiple = []
    plot_single = [0,0,0,0,0,0,0,0]

    #while type != None:
    while count<10:
        if count == 1: 
            new_matrix = matrix
            print("\n","Matriz original: ",count,"\n")
            l1, l2 = len(new_matrix), len(new_matrix[0])
            print(pd.DataFrame(new_matrix, index=['']*l1, columns=['']*l2),"\n")
            last_comb = 0

        valid = combinations(new_matrix)
        
        if all(v is None for v in valid): # Se não houver combinações válidas o loop para
            break

        comb = next(item for item in valid if item is not None) # Retorna a primeira combinação existente.
        
        if comb == last_comb:
            break
        
        last_comb = comb

        if Qx != None and Fx != None and count == 1: # Caso o user passe uma combinação inicial
            Q_x = Qx+1
            F_x = Fx+1
        elif (Qx != None and Fx == None) or (Qx == None and Fx != None):
            print("A função requer tanto Qx como Fx para funcionar, caso passe uma combinação inicial.")
        else: # Caso o user não passe uma combinação inicial -> QMTOxFMTO default
            Q_x = comb[0]+1
            F_x = comb[1]+1

        # Qx, Fx, Qx_valor, Fx_valor
        plot_single = [Q_x, F_x, new_matrix[Q_x-1][1], new_matrix[F_x+1][1],0,0,0,0] # Sabemos as entradas, mas ainda não as saídas
    
        new_matrix = perform_transform(new_matrix, Q_x, F_x, delta_T_min)
        
        print("------------------------------")
        print("\n","Número da troca: ",count,"\n")
        print("Q", comb[0]+1, "F", comb[1]+1,"\n", sep='')

        plot_single[4] = new_matrix[Q_x-1][1] # Saída quente -> TSQ_valor
        plot_single[5] = new_matrix[F_x+1][1] # Saída fria -> TSF_valor
        plot_single[6] = new_matrix[Q_x-1][0] # Wcp_Q
        plot_single[7] = new_matrix[F_x+1][0] # Wcp_F

        l1, l2 = len(new_matrix), len(new_matrix[0])

        #print(plot_single)
        #print("Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor")
        print(pd.DataFrame(new_matrix, index=['']*l1, columns=['']*l2),"\n")

        plot_multiple.append(plot_single)

        count+=1
    
    if plot == True:
        return plot_multiple, new_matrix

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
    
    # Chamar RPS
    plot_multiple, last_matrix =  perform_RPS(matriz, plot=True)

    plot_multiple_transfers(plot_multiple, last_matrix)

def exibir_matriz(matriz):
    for i in range(4):
        for j in range(3):
            print(matriz[i][j], end="\t")

input_user = False

matriz =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
         [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
         [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
         [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

if input_user == True:
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
else:
    plot_multiple, last_matrix =  perform_RPS(matriz, Qx=0, Fx=1, plot=True)
    print(plot_multiple)
    print(last_matrix)
    plot_multiple_transfers(plot_multiple, last_matrix)