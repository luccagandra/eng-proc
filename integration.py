import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

class IntegraçãoEnergética:
    def __init__(self, matrix, delta_T_min=10, Qx=None, Fx=None, U_tc=0.75, U_re=0.75, U_aq=1.00,custo_unit_agua=0.00005,custo_unit_vapor=0.0015,temp_vapor=(250,250),temp_agua=(30,50),  user_input=False):
        # Define propriedades e variáveis que serão usadas nos cálculos
        self.matrix = matrix
        self.user_input = user_input
        self.prev_matrix = [[0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0], 
                            [0.0, 0.0, 0.0]]
        self.custo_cap = 0
        self.custo_util = 0
        self.custo_min_vapor = 0
        self.custo_min_agua = 0
        self.Q = 0
        self.last_comb = (-10,-10)
        self.is_there_two_chains = False
        self.chains_not_in_the_loop = (-1,-1)
        self.plot = []
        self.iterations = 0
        self.coordinates = [0,0]
        self.F_coords= [[0,0],
                        [0,0]]
        self.Q_coords = [[0,0],
                        [0,0]]
        self.comb_history = []
        self.used = [0,0,0,0]
        self.intervals = []

        self.Qx_init = Qx
        self.Fx_init = Fx
        self.temp_vapor = temp_vapor # T_e, T_s - Mudar caso a corrente de vapor apresente temperaturas diferentes
        self.temp_agua = temp_agua # T_e, T_s - Mudar caso a corrente de água apresente temperaturas diferentes
        self.delta_T_min = delta_T_min
        self.custo_unit_agua = custo_unit_agua # Mudar caso o custo da água for diferente
        self.custo_unit_vapor = custo_unit_vapor # Mudar caso o custo do vapor for diferente
        self.U_tc = U_tc # Mudar caso o coeficiente de troca de calor global do trocador de integração for diferente
        self.U_re = U_re # Mudar caso o coeficiente de troca de calor global do resfriador for diferente
        self.U_aq = U_aq # Mudar caso o coeficiente de troca de calor global do aquecedor for diferente

    def perform_transform(self, Qx, Fx):
        # Realiza os cálculos de um trocador de calor
        self.atualizar_prev_matrix(self.matrix)

        QMT0 = self.matrix[Qx-1].copy()
        FMT0 = self.matrix[Fx+1].copy()

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
        if TEQ_fixado - TSF_meta < self.delta_T_min:
            TSF_meta = TEQ_fixado - self.delta_T_min

        # Se saída quente - entrada fria é menor que ΔT(min)
        if TSQ_meta - TEF_fixado < self.delta_T_min:
            TSQ_meta = TEF_fixado + self.delta_T_min

        oferta =  WCp_Q * (TEQ_fixado - TSQ_meta)
        demanda = WCp_F * (TSF_meta - TEF_fixado)

        # Q

        self.Q = min(oferta, demanda)

        if self.Q == oferta:
            TSQ_meta = TSQ_meta # Temperatura de saída quente confirmada
            TSF_meta = TEF_fixado + round(self.Q/WCp_F,1) 
        
        # Se Q = demanda, TSQ = entrada quente - Q / WCp_F
        if self.Q == demanda:
            TSF_meta = TSF_meta
            TSQ_meta = TEQ_fixado - round(self.Q/WCp_Q,1) 

        # TSF_meta e TSQ_meta são atualizados
        QMT0[1] = TSQ_meta
        FMT0[1] = TSF_meta

        new_matrix = self.matrix.copy()

        new_matrix[Qx-1] = QMT0
        new_matrix[Fx+1] = FMT0

        self.atualizar_matrix(matrix)

        return new_matrix

    def next_combination(self, tipo):
        # Determina quais as trocas são válidas e quais não são
        # Seleciona a troca a ser realizada no caso da aplicação dum um método
        
        # Temperaturas de entrada das corrents quentes 
        # ToQ[0] = Q1, ToQ[1] = Q2
        ToQ = (self.matrix[0][1], self.matrix[1][1])

        # Temperaturas de entrada das corrents frias
        # ToF[0] = F1, ToF[1] = F2   
        ToF = (self.matrix[2][1], self.matrix[3][1])

        index_max_ToQ = np.argmax(ToQ)
        index_min_ToQ = np.argmin(ToQ)
        index_max_ToF = np.argmax(ToF)
        index_min_ToF = np.argmin(ToF)

        max_combination = (index_max_ToQ,index_max_ToF)
        min_combination = (index_min_ToQ,index_min_ToF)
        mix_1 = (index_min_ToQ,index_max_ToF)
        mix_2 = (index_max_ToQ,index_min_ToF)

        possible_combinations = []

        # Combinações válidas para Q1
        for i in range(len(ToQ)): #(ToQ[0] > ToF[i]+10)
            if (ToQ[0] > ToF[i]) and (self.matrix[0][1] != self.matrix[0][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append((0,i))

        # Combinações válidas para Q2
        for i in range(len(ToQ)): #(ToQ[1] > ToF[i]+10)
            if (ToQ[1] > ToF[i]) and (self.matrix[1][1] != self.matrix[1][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append((1,i))

        QMTOxFMTO = max_combination if max_combination in possible_combinations else None
        QmTOxFmTO = min_combination if min_combination in possible_combinations else None
        QmTOxFMTO = mix_1 if mix_1 in possible_combinations else None
        QMTOxFmTO = mix_2 if mix_2 in possible_combinations else None

        if tipo == 'QMTOxFMTO':
            combinations_ranked = QMTOxFMTO, QMTOxFmTO,QmTOxFMTO, QmTOxFmTO,

        if tipo == 'QmTOxFmTO':
            combinations_ranked = QmTOxFmTO, QmTOxFMTO, QMTOxFmTO, QMTOxFMTO

        if tipo == 'QmTOxFMTO':
            combinations_ranked = QmTOxFMTO, QmTOxFmTO, QMTOxFmTO, QMTOxFMTO

        if tipo == 'QMTOxFmTO':
            combinations_ranked = QMTOxFmTO, QMTOxFMTO, QmTOxFMTO, QmTOxFmTO

        """
        print("Combinações possíveis: ", possible_combinations)
        print("QMTOxFMTO:","Q", max_combination[0], "F", max_combination[1])
        print("QmTOxFmTO:","Q", min_combination[0], "F", min_combination[1])
        """

        if len(possible_combinations) == 0:
            QMTOxFMTO = None
            QmTOxFmTO = None

        are_all_none = all(item is None for item in combinations_ranked)

        if are_all_none:
            comb = None
        else:   
            comb = next(item for item in combinations_ranked if item is not None) # Retorna a primeira combinação existente.

        return comb 

    def valid_combinations(self):
        # Chama o método anterior e as organisa numa lista
        
        ToQ = (self.matrix[0][1], self.matrix[1][1])

        # Temperaturas de entrada das corrents frias
        # ToF[0] = F1, ToF[1] = F2   
        ToF = (self.matrix[2][1], self.matrix[3][1])

        possible_combinations = []

        # Combinações válidas para Q1
        for i in range(len(ToQ)): #(ToQ[0] > ToF[i]+10)
            if (ToQ[0] > ToF[i]) and (self.matrix[0][1] != self.matrix[0][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append([0,i])

        # Combinações válidas para Q2
        for i in range(len(ToQ)): #(ToQ[1] > ToF[i]+10)
            if (ToQ[1] > ToF[i]) and (self.matrix[1][1] != self.matrix[1][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append([1,i])
        
        return possible_combinations

    def atualizar_matrix(self, matrix):
        self.matrix = matrix

    def atualizar_prev_matrix(self, prev_matrix):
        self.prev_matrix = prev_matrix

    def draw_arrow(self, text,x, start, end, color_arrow):
        arrowprops = dict(arrowstyle='->', linestyle='-', linewidth=2, color=color_arrow)
        plt.annotate(text, xy=(x, end), xytext=(x, start), arrowprops=arrowprops, color=color_arrow)

    def criar_grafico(self):
        y1 = [self.matrix[0][1], self.matrix[0][2], self.matrix[1][1], self.matrix[1][2]]  # Q1_To, Q1_Td, Q2_To, Q2_Td
        y2 = [self.matrix[2][1], self.matrix[2][2], self.matrix[3][1], self.matrix[3][2]]  # F1_To, F1_Td, F2_To, F2_Td

        x1 = [0, 0.5, 0.5, 1]  # Pontos de mudança nos degraus para y1
        
        y_ticks = []

        for i in y1:
            plot_y1_step = [i,i,i-10, i-10]
            plt.step(x1, plot_y1_step)

            y_ticks.append(i)
            y_ticks.append(i-10)
            self.intervals.append(i-10)

        for i in y2:
            plot_y2_step = [i+10,i+10,i, i]
            plt.step(x1, plot_y2_step)
            y_ticks.append(i)
            y_ticks.append(i+10)
            self.intervals.append(i)

        self.intervals = list(set(self.intervals)) # sort and remove duplicates
        self.intervals.sort(reverse=True)
        #print(intervals)
        plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
        plt.xticks([])

        self.draw_arrow("",0.2,self.matrix[0][1],self.matrix[0][2], "red")
        self.draw_arrow("",0.4,self.matrix[1][1],self.matrix[1][2], "red")

        self.draw_arrow("",0.6,self.matrix[2][1],self.matrix[2][2], "blue")
        self.draw_arrow("",0.8,self.matrix[3][1],self.matrix[3][2], "blue")

        plt.grid(linestyle = '--')
        plt.show()

    def display_table(self, data):
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

    def is_pair_overlapping(self, pair1, pair2):
        # Sort the pairs to ensure proper comparison
        pair1 = sorted(pair1)
        pair2 = sorted(pair2)

        # Check if pair1 is completely within pair2
        return pair1[0] >= pair2[0] and pair1[1] <= pair2[1]
    
    def offer_demand(self):
        # Calcula o método do transbordo para obter o custo mínimo de utilidades
        
        Rk = 0
        a1 = 1
        a2 = 1
        offer_demand = []

        for i in range(len(self.intervals)-1): 
            
            if self.is_pair_overlapping((self.intervals[i]+self.delta_T_min, self.intervals[i+1]+self.delta_T_min), (self.matrix[0][2], self.matrix[0][1])): # passa por arrow1
                a1 = self.matrix[0][0]
            else:
                a1 = 0

            if self.is_pair_overlapping((self.intervals[i]+self.delta_T_min, self.intervals[i+1]+self.delta_T_min), (self.matrix[1][2], self.matrix[1][1])): # passa por arrow2
                a2 = self.matrix[1][0]
            else:
                a2 = 0

            if self.is_pair_overlapping((self.intervals[i], self.intervals[i+1]), (self.matrix[2][2], self.matrix[2][1])): # passa por arrow3
                a3 = self.matrix[2][0]
            else:
                a3 = 0

            if self.is_pair_overlapping((self.intervals[i], self.intervals[i+1]), (self.matrix[3][2], self.matrix[3][1])): # passa por arrow4
                a4 = self.matrix[3][0]
            else:
                a4 = 0
            
            y = (self.intervals[i]-self.intervals[i+1])
            lines = [Rk, y*a1 +y*a2, y*a3 +y*a4, Rk + y*a1 +y*a2 - y*a3 - y*a4]

            if lines[3] < 0: # pinch
                Rk = 0 
                self.custo_util_min(abs(lines[3]), self.U_aq, "aquecedor")
                print(f"Custo mínimo de utilidades(vapor):  {abs(round(lines[3],2))}kW, ou {round(self.custo_min_vapor,2)}$")
            else:
                Rk = lines[3]

                if i == len(self.intervals)-2:
                    self.custo_util_min(abs(lines[3]), self.U_re, "resfriador")
                    print(f"Custo mínimo de utilidades(água): {abs(round(lines[3],2))}kW, ou {round(self.custo_min_agua,2)}$")
            
            offer_demand.append(lines)

        print(f"Custo mínimo de utilidades(total): {round(self.custo_min_agua+self.custo_min_vapor,2)}$")
        self.display_table(offer_demand)
    
    def loop_RPS(self, tipo):
        # Aplica o método do RPS
        
        while True:
            if self.user_input == False:
                if self.iterations == 0 and self.Qx_init != None and self.Fx_init != None:
                    comb = (self.Qx_init, self.Fx_init)
                else:    
                    comb = self.next_combination(tipo)

                if comb == None: # Primeira condição de quebra -> Deve haver combinações possíveis
                    break

                if comb[0] != self.last_comb[0] and comb[1] != self.last_comb[1] and self.iterations != 0:
                    self.is_there_two_chains = True
                    self.chains_not_in_the_loop = (comb[0],comb[1])
                    break

                if comb in self.comb_history: # Segunda condição de quebra -> As combinações não podem ser iguais
                    break

            else:
                if self.iterations == 0 and self.Qx_init != None and self.Fx_init != None:
                    comb = (self.Qx_init, self.Fx_init)
                else: 
                    valid = self.valid_combinations()
                    
                    are_all_none = all(item is None for item in valid) # Primeira condição de quebra -> Deve haver combinações possíveis

                    if are_all_none:
                        print("\n","------------------------------","\n")
                        print("NÃO HÁ MAIS COMBINAÇÕES POSSÍVEIS. COMPLETANDO COM UTILIDADES.") 
                        break
                    else:
                        for i in range(len(valid)):
                            print(f"Combinação válida: Q{valid[i][0]+1}xF{valid[i][1]+1}")
                    
                        comb = [0,0]

                        comb[0] = int(input("Escolha a corrente quente: Q")) -1
                        comb[1] = int(input("Escolha a corrente fria: F")) -1

                        if comb[0] != self.last_comb[0] and comb[1] != self.last_comb[1] and self.iterations != 0: # Segunda condição de quebra -> As combinações não podem ser iguais
                            print("\n","------------------------------","\n")
                            print("COMBINAÇÃO FORA DA REDE. COMPLETANDO COM UTILIDADES A PRIMEIRA E CRIANDO OUTRA REDE.") 

                            self.is_there_two_chains = True
                            self.chains_not_in_the_loop = (comb[0],comb[1])
                            break
            
            self.used[comb[0]] = True
            self.used[comb[1]+2] = True

            new_matrix = self.perform_transform(comb[0]+1, comb[1]+1)
            self.atualizar_matrix(new_matrix)

            self.print_comb(comb, new_matrix)

            # Com a nova matriz, devemos calcular o custo do trocador
            
            TEQ = self.prev_matrix[comb[0]][1]
            TSQ = new_matrix[comb[0]][1]
            TEF = self.prev_matrix[comb[1]+2][1]
            TSF = new_matrix[comb[1]+2][1]

            self.atualizar_custo_cap(TEQ, TSF, TSQ, TEF, self.U_tc)

            Qx = comb[0]+1
            Fx = comb[1]+1

            self.plot.append([TEQ, TSQ, TEF, TSF, Qx, Fx]) #self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx])

            self.comb_history.append(comb)
            self.last_comb = comb

            self.iterations += 1

    def completando_utilidades(self):

        matrix_arr = np.array(self.matrix)
        
        for i in range(2): # Utilidades para Q1 e Q2
            """
            Obs: Completar com utilidades é simplesmente fazer uma troca de calor como
            qualquer outra, com a diferença de que o aquecimento ouresfriamento não será
            com uma corrente específica, mas sim com uma de água. Dessa forma, a restrição
            de delta_T_min continua valendo.
            """ 

            # Se entrada quente - saída fria é menor que ΔT(min)
            """
            Exemplo: a meta de resfriamento é 0 graus, mas a temp de entrada e saída da água é, 30 e 50.
            Logo, só podemos resfriar a corrente até 30. A meta de saída da corrente é mudada.
            """

            WCp = matrix_arr[i,0]
            TEQ = matrix_arr[i,1]
            TSQ = matrix_arr[i,2]
            TEF = self.temp_agua[0]
            TSF = self.temp_agua[1]

            self.Q = WCp * (TEQ - TSQ)
            
            if TSQ - self.temp_agua[1] < self.delta_T_min and self.used[i] == True: # Ex: 50 - 50 < 10, TSQ = 50+10
                print(f"Q{i} terá nova meta de resfriamento: {self.temp_agua[1] + self.delta_T_min}ºC, ao invés de {TSQ}ºC, por causa de ΔT(min)={self.delta_T_min}.")
                TSQ = self.temp_agua[1] + self.delta_T_min 

            # Se ainda não chegou na meta definida
            if TEQ > TSQ and self.used[i] == True:
                self.iterations += 1

                Qx = i+1
                Fx = None
                
                self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx]) #self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx])

                self.print_utilidades_resfriamento(f"{i+1}", TSQ, matrix_arr[i])
                self.atualizar_custo_cap(TEQ, TSF, TSQ, TEF, self.U_re)
                self.atualizar_custo_util(WCp, TEQ, TSQ, self.U_re, "resfriador")
        
        for i in range(2): # Utilidades para F1 e F2
            WCp = matrix_arr[i+2,0]
            TEF = matrix_arr[i+2,1]
            TSF = matrix_arr[i+2,2]
            TEQ = self.temp_vapor[0]
            TSQ = self.temp_vapor[1]

            self.Q = WCp * (TSF - TEF)

            if self.temp_vapor[1] - TSF < self.delta_T_min and self.used[i+2] == True: # Ex: 250 - 250 < 10, TSF = 250+10
                print(f"F{i+1} terá nova meta de resfriamento: {self.temp_vapor[1] + self.delta_T_min}ºC, ao invés de {TSF}ºC, por causa de ΔT(min)={self.delta_T_min}.")
                TSF = self.temp_vapor[1] - self.delta_T_min
            
            # Se ainda não chegou na meta definida
            if TEF < TSF and self.used[i+2] == True: 
                self.iterations += 1

                Qx = None
                Fx = i+1

                self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx]) #self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx])

                self.print_utilidades_aquecimento(f"{i+1}", TSF, matrix_arr[i+2])
                self.atualizar_custo_cap(TEQ, TSF, TSQ, TEF, self.U_aq)
                self.atualizar_custo_util(WCp, TSF, TEF, self.U_aq, "aquecedor")

    def print_comb(self, comb, new_matrix):
        
        print("\n","------------------------------","\n")

        print(f"(TROCADOR DE CALOR) Combinação feita: Q{comb[0]+1}xF{comb[1]+1}")

        print(pd.DataFrame(new_matrix, index=['']*len(new_matrix), columns=['']*len(new_matrix[0])),"\n")

    def print_utilidades_resfriamento(self, id, TSQ, corrente):
        
        print("\n","------------------------------","\n")

        print(f"(RESFRIADOR) Resfriamento feito em: Q{id} \n")

        print(pd.DataFrame(self.matrix, index=['']*len(self.matrix), columns=['']*len(self.matrix[0])),"\n")

        print(corrente, "\n")

        print(f"Resfriamento até {TSQ}ºC com água de T(in)={self.temp_agua[0]}ºC e T(out)={self.temp_agua[1]}ºC \n")

    def print_utilidades_aquecimento(self, id, TSF, corrente):
        
        print("\n","------------------------------","\n")

        print(f"(AQUECEDOR) Aquecimento feito em: F{id}")

        print(pd.DataFrame(self.matrix, index=['']*len(self.matrix), columns=['']*len(self.matrix[0])),"\n")

        print(corrente, "\n")

        print(f"Aquecimento até {TSF}ºC com vapor de T(in)={self.temp_vapor[0]}ºC e T(out)={self.temp_vapor[1]}ºC")

    def atualizar_custo_cap(self, TEQ, TSF, TSQ, TEF, U):

        delta_1 = TEQ - TSF 
        delta_2 = TSQ - TEF 

        print("TEQ:",TEQ,"TSQ:",TSQ,"TEF:",TEF,"TSF:",TSF)

        #if tipo == 'aritmético':
        #    area = Q / (U*(delta_1+delta_2)/2)
        
        #if tipo == 'logarítmico':
        if delta_1 == delta_2: # Não podemos passar na fórmula abaixo pq gera uma indeterminação
            area = delta_1
        else:
            area = self.Q / (U*(delta_1-delta_2)/np.log(delta_1/delta_2))

        # Mudar caso o modelo do custo de capital for diferente
        custo_do_trocador = 130 * (math.pow(area, (65/100))) # EQUAÇÃO DO MODELO DE CUSTO DO TROCADOR -> SLIDE

        self.custo_cap += custo_do_trocador
        print("Novo custo cap:", round(self.custo_cap,2))

    def atualizar_custo_util(self, WCp, T_in, T_out, U, tipo):

        Q = WCp*(T_in - T_out)

        # W = vazão

        if tipo == "resfriador":

            W = Q/((0.00116)*(self.temp_agua[1] - self.temp_agua[0])) # 0.00116 = Cp_água
            # Mudar caso o modelo do custo de utilidades for diferente
            self.custo_util += 8500*(self.custo_unit_agua*W)

        if tipo == "aquecedor":

            W = Q/0.48 # 0.48 = Cp_vapor
            # Mudar caso o modelo do custo de utilidades for diferente
            self.custo_util += 8500*(self.custo_unit_vapor*W)
        
        print("Novo custo útil:", round(self.custo_util,2),"\n")

    def custo_util_min(self, Q, U, tipo):

        if tipo == "resfriador":

            W = Q/((0.00116)*(self.temp_agua[1] - self.temp_agua[0])) # 0.00116 = Cp_água
            # Mudar caso o modelo do custo de utilidades for diferente
            self.custo_min_agua = 8500*(self.custo_unit_agua*W)

        if tipo == "aquecedor":

            W = Q/0.48 # 0.48 = Cp_vapor
            # Mudar caso o modelo do custo de utilidades for diferente
            self.custo_min_vapor = 8500*(self.custo_unit_vapor*W)

    def plot_single(self, plot_list, id_plot, last):

        TEQ = plot_list[0] # self.plot.append([TEQ, TSQ, TEF,TSF, Qx, Fx])
        TSQ = plot_list[1]
        TEF = plot_list[2]
        TSF = plot_list[3]
        Qx = plot_list[4]
        Fx = plot_list[5]
        segment_length = 1
        offset = 0.1
        x = self.coordinates[0]
        y = self.coordinates[1]
        
        if (Qx and Fx) != None: # Trocador
            plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
            plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')            

            text_hot = f'Q{Qx} = {TEQ}'
            text_cold = f'F{Fx} = {TEF}'

            # Adding labels 
            plt.text(x - segment_length/2, y + offset, text_hot, color='red', ha='right', va='center', fontsize=10)
            plt.text(x + offset, y + segment_length/2, text_cold, color='blue', ha='center', va='bottom', fontsize=10)

            plt.text(x + segment_length/2, y + offset, f'{TSQ}', color='red', ha='left', va='center', fontsize=10)
            plt.text(x + offset, y - segment_length/2, f'{TSF}', color='blue', ha='center', va='top', fontsize=10)

        if Fx == None: # Resfriador 
            plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='blue')
            plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
            
            plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, TSF, color='blue', ha='left', va='center', fontsize=10)
            plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, TEF, color='blue', ha='right', va='center', fontsize=10)
            plt.text(x + segment_length/2, y + offset, TSQ, color='red', ha='left', va='center', fontsize=10)

        if Qx == None: # Aquecedor 
            plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='red')
            plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')

            plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, TEQ, color='red', ha='left', va='center', fontsize=10)
            plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, TSQ, color='red', ha='right', va='center', fontsize=10)
            plt.text(x + offset, y - segment_length/2, TSF, color='blue', ha='center', va='top', fontsize=10)
        
        # Plotting the number slightly offset from the center
        plt.text(x + offset/2, y + offset/2, str(id_plot), color='black', ha='center', va='center', fontsize=12, fontweight='bold')

        if last == True:
            # Adding labels, title, legend, and grid for better visualization
            plt.title('Rede Final')
            plt.grid(False)
            plt.axis('off')

            text1 = f"Ccap = {round(self.custo_cap,2)}"
            text2 = f"Cutil = {round(self.custo_util,2)}"
            text3 = f"Ctotal = {round(self.custo_cap + self.custo_util,2)}"

            plt.text(2.5, 0.8, text1, ha='right', va='top', color='red', fontsize=10)
            plt.text(2.5, 0.6, text2, ha='right', va='top', color='green', fontsize=10)
            plt.text(2.5, 0.4, text3, ha='right', va='top', color='blue', fontsize=10)         

            #print(custo_cap, custo_util, custo_total)
            plt.show()

    def plot_multiple(self):
        for i in range(len(self.plot)): # TEQ, TSQ, TEF,TSF, Qx, Fx
            if i == 0:
                self.plot_single(self.plot[0], 1, False) # plot_single(self, plot_list, id_plot, last)
            else:
                
                self.update_coordinates(self.plot[i-1], self.plot[i])
                
                last = False

                if i == len(self.plot)-1:
                    last = True
                
                self.plot_single(self.plot[i], i+1, last)

    def update_coordinates(self, prev, current):
        
        Qx_prev = prev[4]
        Qx_current = current[4]
        Fx_prev = prev[5]
        Fx_current = current[5]

        if Qx_prev == Qx_current and Qx_current != None and Fx_current != None: # Qx_prev = Qx_current -> conexão trocador por corrente quente
            self.F_coords[Fx_prev-1][0] = self.coordinates[0]
            self.coordinates[0] += 1
            self.F_coords[Fx_current-1][0] = self.coordinates[0]
            self.F_coords[Fx_current-1][1] = self.coordinates[1]
            self.Q_coords[Qx_current-1][0] = self.coordinates[0]
            self.Q_coords[Qx_current-1][1] = self.coordinates[1]

        if Fx_prev == Fx_current and Qx_current != None and Fx_current != None: # Fx_prev = Fx_current -> conexão trocador por corrente fria
            self.Q_coords[Qx_prev-1][1] = self.coordinates[1]
            self.coordinates[1] -= 1
            self.Q_coords[Qx_current-1][1] = self.coordinates[1]
            self.Q_coords[Qx_current-1][0] = self.coordinates[0]
            self.F_coords[Fx_current-1][1] = self.coordinates[1]
            self.F_coords[Fx_current-1][0] = self.coordinates[0]

        if Fx_current == None: #resfriador
            self.coordinates[0] = self.Q_coords[Qx_current-1][0] + 1
            self.coordinates[1] = self.Q_coords[Qx_current-1][1]      

        if Qx_current == None: #aquecedor
            self.coordinates[0] = self.F_coords[Fx_current-1][0]
            self.coordinates[1] = self.F_coords[Fx_current-1][1] - 1

matrix =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
         [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
         [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
         [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

matriz_aula = [[3.0, 170.0, 60.0], # WCp, Q1_T0, Q1_Td
               [1.5, 150.0, 30.0], # WCp, Q2_T0, Q2_Td
               [2.0, 30.0, 140.0],  # WCp, F1_T0, F1_Td
               [4.0, 80.0, 140.0]] # WCp, F2_T0, F2_Td