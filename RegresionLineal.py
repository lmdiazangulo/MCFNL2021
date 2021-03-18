#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:54:19 2021

@author: dla
"""

def RegresionLineal(x, y):
    '''Funcion que realiza la regresión lineal de la lista de números X e Y
    de igual longitud (ambos listas) y devuelve, si la ecuación de la recta 
    de mejor ajuste es Y = mX + n: m, \Delta{m}, n, \Delta{n}. También 
    devuelve el coeficiente de correlación r (a veces llamado r^2)'''
    
    #0) Modules needed
    import math
    import pandas as pd   
    from statsmodels.formula.api import ols 
    
    #-1) De prueba, valores x e y:
    #x = [1, 2]
    #y = [3, 4]
    
    #1) Sumatorios necesarios
    suma_x = sum(x)
    suma_y = sum(y)
    N = len(x)                                  #longitud de los vectores
        #1.1) Calculo de sumatorios
        
    suma_x2 = sum([x_value **2 for x_value in x])
                        #creation of the array x^2, and then sum its elements
    suma_y2 = sum([y_value **2 for y_value in y])
                        #same as above but w\ y^2
    suma_xy = sum([x[i] * y[i] for i in range(0,N)])
                        #same as above but w\  xy
    
    '''Long way to calculate the sums (matlab influenced)
    suma_x2 = 0                                     #inicializaicon
    suma_y2 = 0                                     #inicializaicon
    suma_xy = 0                                     #inicializaicon              
    
        """
        Way to calculate the sum x^2 and sum y^2, but not sum x*y (Matlab style)
    for element in x:               #calc of sum x^2
        suma_x2 = suma_x2 + element**2
    
    for element in y:                #calc of sum y^2
        suma_y2 = suma_y2 + element**2
        """

    for value in range(0,N):
        suma_x2 = suma_x2 + x[value]**2
        suma_y2 = suma_y2 + y[value]**2
        suma_xy = suma_xy + x[value]*y[value]
    '''
    
    #2) Realizacion del ajuste
    data = pd.DataFrame({'X': x, 'Y': y})  
                                        #creacion de un panda con los datos
    ajuste= ols("Y ~ X", data).fit()                                #ajuste
                    #es fundamental lo de "Y ~ X". Si lo pones al contrario
                    #hace lo contrario!!
                    
    slope = ajuste.params[1]                         #pendiente
    intercept = ajuste.params[0]               #ordenada origen
    r = ajuste.rsquared                     #coef correlacion r
    
    #3) Calculo de errores del ajuste
    delta_slope = 3 * math.sqrt(slope**2/(N-2) * (1/r**2 - 1))
    delta_intercept = 3 * math.sqrt(suma_x2 * delta_slope**2 / N)
            
    #4) Return of values
        #the values will be returned in a dictionary indicating what is each
        #value
    values = {'Slope' : slope, 'Intercept' : intercept, 'r' : r, 
              '\Delta{slope}' : delta_slope, 
              '\Delta{intercept}' : delta_intercept}
    return values