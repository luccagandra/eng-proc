import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")  # Defina o backend para TkAgg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

matrix =[[10.0, 180.0, 90.0], # WCp, Q1_T0, Q1_Td
         [2.0, 250.0, 140.0], # WCp, Q2_T0, Q2_Td
         [5.0, 60.0, 150.0],  # WCp, F1_T0, F1_Td
         [7.0, 100.0, 220.0]] # WCp, F2_T0, F2_Td

"""
exemplo = perform_transform(matrix, 2, 2, 10)
exemplo2 = perform_transform(exemplo, 1, 2, 10)
exemplo3 = perform_transform(exemplo2, 1, 1, 10)

print(exemplo)
print(exemplo2)
print(exemplo3)
"""

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
    
    # Combinações válidas para Q1
    for i in range(len(ToQ)):
        if (ToQ[0] > ToF[i]) and (matrix[0][1] != matrix[0][2]):
            possible_combinations.append((0,i))

    # Combinações válidas para Q2
    for i in range(len(ToQ)):
        if (ToQ[1] > ToF[i]) and (matrix[1][1] != matrix[1][2]):
            possible_combinations.append((1,i))
    
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

def perform_RPS(matrix, Qx=None, Fx=None, delta_T_min=10):

    comb = (0,0)
    count = 1

    #while type != None:
    while count<1000:
        if count == 1: 
            new_matrix = matrix
            print("Matriz original: ",count,"\n")
            l1, l2 = len(new_matrix), len(new_matrix[0])
            print(pd.DataFrame(new_matrix, index=['']*l1, columns=['']*l2),"\n")

        valid = combinations(new_matrix)

        if all(v is None for v in valid): # Se não houver combinações válidas o loop para
            break

        comb = next(item for item in valid if item is not None) # Retorna a primeira combinação existente.
        
        if Qx != None and Fx != None: # Caso o user passe uma combinação inicial
            Q_x = Qx
            F_x = Fx
        elif (Qx != None and Fx == None) or (Qx == None and Fx != None):
            print("A função requer tanto Qx como Fx para funcionar, caso passe uma combinação inicial.")
        else: # Caso o user não passe uma combinação inicial -> QMTOxFMTO default
            Q_x = comb[0]+1
            F_x = comb[1]+1

        new_matrix = perform_transform(new_matrix, Q_x, F_x, delta_T_min)
        
        print("------------------------------")
        print("Número da troca: ",count,"\n")
        print("Q", comb[0]+1, "F", comb[1]+1,"\n", sep='')

        l1, l2 = len(new_matrix), len(new_matrix[0])
        print(pd.DataFrame(new_matrix, index=['']*l1, columns=['']*l2),"\n")
        
        count+=1

perform_RPS(matrix)