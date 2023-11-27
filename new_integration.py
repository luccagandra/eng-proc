import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

class IntegraçãoEnergética:
    def __init__(self, matrix, delta_T_min):
        self.matrix = matrix
        self.delta_T_min = delta_T_min
        self.custo_cap = 0
        self.custo_util = 0
        self.Q = 0

    def perform_transform(self, Qx, Fx):

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
            TSF_meta = TEQ_fixado - 10

        # Se saída quente - entrada fria é menor que ΔT(min)
        if TSQ_meta - TEF_fixado < self.delta_T_min:
            TSQ_meta = TEF_fixado + 10

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

        return new_matrix

    def next_combination(self, tipo):
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

        possible_combinations = []
        
        max_combination = (index_max_ToQ,index_max_ToF)
        min_combination = (index_min_ToQ,index_min_ToF)
        mix_1 = (index_min_ToQ,index_max_ToF)
        mix_2 = (index_max_ToQ,index_min_ToF)
        #print("self.matrix: ", self.matrix)
        #print("ToQ: ", ToQ)
        #print("ToF: ", ToF)

        # Combinações válidas para Q1
        for i in range(len(ToQ)):
            if (ToQ[0] > ToF[i]+10) and (self.matrix[0][1] != self.matrix[0][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append((0,i))

        # Combinações válidas para Q2
        for i in range(len(ToQ)):
            if (ToQ[1] > ToF[i]+10) and (self.matrix[1][1] != self.matrix[1][2]) and (self.matrix[2+i][1] != self.matrix[2+i][2]):
                possible_combinations.append((1,i))

        QMTOxFMTO = max_combination if max_combination in possible_combinations else None
        QmTOxFmTO = min_combination if min_combination in possible_combinations else None
        QmTOxFMTO = mix_1 if mix_1 in possible_combinations else None
        QMTOxFmTO = mix_2 if mix_2 in possible_combinations else None
        
        if tipo == 'QMTOxFMTO':
            combinations_ranked = QMTOxFMTO, QmTOxFmTO, QmTOxFMTO, QMTOxFmTO

        if tipo == 'QmTOxFmTO':
            combinations_ranked = QmTOxFmTO, QMTOxFMTO, QmTOxFMTO, QMTOxFmTO

        if tipo == 'QmTOxFMTO':
            combinations_ranked = QmTOxFMTO, QmTOxFmTO, QMTOxFMTO, QMTOxFmTO

        if tipo == 'QMTOxFmTO':
            combinations_ranked = QMTOxFmTO, QmTOxFmTO, QMTOxFMTO, QmTOxFMTO

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
        
        comb = next(item for item in combinations_ranked if item is not None) # Retorna a primeira combinação existente.

        return comb 

    def print_comb(self, comb, new_matrix):
        
        print("\n","------------------------------","\n")

        print(f"Combinação feita: Q{comb[0]+1}xF{comb[1]+1}")

        print(pd.DataFrame(new_matrix, index=['']*len(new_matrix), columns=['']*len(new_matrix[0])),"\n")

    def atualizar_custo_cap(self, TEQ, TSF, TSQ, TEF, U):

        delta_1 = TEQ - TSF 
        delta_2 = TSQ - TEF 

        #if tipo == 'aritmético':
        #    area = Q / (U*(delta_1+delta_2)/2)
        
        #if tipo == 'logarítmico':
        if delta_1 == delta_2: # Não podemos passar na fórmula abaixo pq gera uma indeterminação
            area = delta_1
        else:
            area = self.Q / (U*(delta_1-delta_2)/np.log(delta_1/delta_2))

        custo_do_trocador = 130 * (math.pow(area, (65/100)))

        self.custo_cap += custo_do_trocador
        print("Novo custo cap:", round(self.custo_cap,2))

    def atualizar_matrix(self, matrix):
        self.matrix = matrix

matrix =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
         [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
         [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
         [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

delta_T_min = 10
last_comb = 0

loop = IntegraçãoEnergética(matrix, delta_T_min)

print("\n","Matriz original: \n")
print(pd.DataFrame(matrix, index=['']*len(matrix), columns=['']*len(matrix[0])),"\n")

# Quais as combinações possíveis?

while True:
    comb = loop.next_combination("QMTOxFMTO")

    if comb == None: # Primeira condição de quebra -> Deve haver combinações possíveis
        break

    if comb == last_comb: # Segunda condição de quebra -> As combinações não podem ser iguais
        break
    last_comb = comb

    new_matrix = loop.perform_transform(comb[0]+1, comb[1]+1)
    loop.atualizar_matrix(new_matrix)

    loop.print_comb(comb, new_matrix)

    # Com a nova matriz, devemos calcular o custo do trocador
    
    TEQ = matrix[comb[0]][1]
    TSQ = new_matrix[comb[0]][2]
    TEF = matrix[comb[1]+2][1]
    TSF = new_matrix[comb[1]+2][2]

    print("TEQ:",TEQ,"TSF:",TSF,"TSQ:",TSQ,"TEF:",TEF)

    loop.atualizar_custo_cap(TEQ, TSF, TSQ, TEF, 0.75) # U = 0.75 pois é trocador de calor

