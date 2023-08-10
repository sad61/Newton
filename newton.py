import inspect, re
import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
from prettytable import PrettyTable as pt

def f(M):
    IMPULSO_DO_FOGUETE = 283.5
    GRAVIDADE_DA_TERRA = 9.81
    VELOCIDADE_DE_EXAUSTAO = IMPULSO_DO_FOGUETE * GRAVIDADE_DA_TERRA
    MASSA_SEM_COMBUSTIVEL = 137000
    ANGULO_DO_FOGUETE = 80
    TEMPO_DE_QUEIMA = 161
    VELOCIDADE_DE_ESCAPE_DA_TERRA = 11200
    FRACAO_DA_VELOCIDADE_DE_ESCAPE = 0.5

    c = VELOCIDADE_DE_EXAUSTAO
    m = MASSA_SEM_COMBUSTIVEL
    theta = np.radians(ANGULO_DO_FOGUETE)
    t = TEMPO_DE_QUEIMA
    g = GRAVIDADE_DA_TERRA
    v = VELOCIDADE_DE_ESCAPE_DA_TERRA * FRACAO_DA_VELOCIDADE_DE_ESCAPE

    return ((c * np.log(M / m) * np.sin(theta)) - (t * g)) - v

def df(M):
    IMPULSO_DO_FOGUETE = 283.5
    GRAVIDADE_DA_TERRA = 9.81
    VELOCIDADE_DE_EXAUSTAO = IMPULSO_DO_FOGUETE * GRAVIDADE_DA_TERRA
    ANGULO_DO_FOGUETE = 80
    
    c = VELOCIDADE_DE_EXAUSTAO
    theta = np.radians(ANGULO_DO_FOGUETE)
    
    return c*np.sin(theta)/M

def newton(M0,f,df,tol,N):

    table = pt()
    
    table.field_names = ['i','M','f(M)','df(M)','e']
 
    for i in range(0,N):
        
        M = M0 - f(M0) / df(M0) 
        
        e = abs(M-M0)/abs(M)

        table.add_row([i, np.round(M,8), np.round(f(M),8), np.round(df(M),4), f'{e:e}'])    
        
        if e < tol:
            break
        M0 = M                
    table.align = 'c';  
    print(table)
       
    if i == N:
        print(f'Solução não obtida em {N:d} iterações')
    else:
        print(f'Solução obtida: M = {M:.6f}')

    return M

val = newton(1000000, f, df, 1e-4, 30)