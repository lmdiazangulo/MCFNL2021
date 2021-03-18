#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 10:28:13 2021

@author: dla


SCRIPT PRACTICA ALFA MASTER, PASADO DE MATLAB A PYTHON

"""

#   %reset                                      #to delete all the variables

#########################1), Data loading #####################3
#The files to load are in .mca. Have find out ot read the data. They key was
#finding the appropiate codification, which was find in the GitHub mca reader.
#the standar, utf-8 do not fits here. The problem was with "Board Temp: 37\B0C",
#line 4151

with open('Am_con_Al_no_vacio_no_pol.mca', encoding = 'iso-8859-15') as file_object:    #file_object stores
            # the object representing the file, which is opened with
            #the open function
            
        Am_Al_mal_data = file_object.read()  
                #This contains strings (have to be converted to numbers using int()
                #and \n, so the \n (salto de linea) have to be removed
        
        Am_Al_mal_data = Am_Al_mal_data.split('\n')         #creating a list
                #split turns a string into a list. Between brackets goes the 
                #sepparator.
                
        while '' in Am_Al_mal_data:   #removal of empty lines
            Am_Al_mal_data.remove('')    
            
        #We have everything of the .mca file in a list. It contains
        #the counts, but also info in the form of text. We can choose the numbers
        #easily since they are wlawys from column 13 to 4108:
        Am_Al_mal = Am_Al_mal_data[12:4108]    #the counts (12 since the first
        #column is 0, not 1!!!)
        
        Am_Al_mal = list(map(int, Am_Al_mal))       #converting it to numbers
        channels = list(range(1, len(Am_Al_mal)+1))     #number of chanels
                #list(range(1,5)) gives [1,2,3,4]
                #same channels as if you use the .txt files
 
                       
with open('Am_con_Al.mca', encoding = 'iso-8859-15') as file_object:
    Am_Al_data = file_object.read()
    Am_Al_data = Am_Al_data.split('\n')
    
    while '' in Am_Al_data:   #removal of empty lines
        Am_Al_data.remove('')
    
    Am_Al = Am_Al_data[12:4108]                     #choosing only the counts    
    Am_Al = list(map(int, Am_Al))       #converting it to numbers
   
     
with open('americioAireNopol.mca', encoding = 'iso-8859-15') as file_object:
    Am_mal_data = file_object.read()
    Am_mal_data = Am_mal_data.split('\n')   
    
    while '' in Am_mal_data:   #removal of empty lines
        Am_mal_data.remove('')

    Am_mal = Am_mal_data[12:4108]                     #choosing only the counts    
    Am_mal = list(map(int, Am_mal))       #converting it to numbers    


with open('americioVacio_pol.mca', encoding = 'iso-8859-15') as file_object:
    Am_data = file_object.read()
    Am_data = Am_data.split('\n')   
    
    while '' in Am_data:   #removal of empty lines
        Am_data.remove('')
    
    Am = Am_data[12:4108]                     #choosing only the counts    
    Am = list(map(int, Am))       #converting it to numbers  
    
    
with open('muestra_problema_aire_no_pol.mca', encoding = 'iso-8859-15') as file_object:
    problema_mal_data = file_object.read()
    problema_mal_data = problema_mal_data.split('\n')   
    
    while '' in problema_mal_data:   #removal of empty lines
        problema_mal_data.remove('')

    problema_mal = problema_mal_data[12:4108]         #choosing only the counts        
    problema_mal= list(map(int, problema_mal))       #converting it to numbers
   
    
with open("muestra_problema.mca", encoding = 'iso-8859-15') as file_object:
    problema_data = file_object.read()
    problema_data = problema_data.split('\n')
    
    while '' in problema_data:   #removal of empty lines
        problema_data.remove('')

    problema = problema_data[12:4108]         #choosing only the counts            
    problema = list(map(int, problema))       #converting it to numbers  
    
    
with open('pulser.mca', encoding = 'iso-8859-15') as file_object:
    pulser_data = file_object.read()
    pulser_data = pulser_data.split('\n')
    
    while '' in pulser_data:   #removal of empty lines
        pulser_data.remove('')

    pulser = pulser_data[12:4108]         #choosing only the counts                
    pulser = list(map(int, pulser))       #converting it to numbers
    
    
    
    
'''##@@@@Panda to read@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import pandas as pd   #to read files

Am_Al_mal_panda = pd.read_csv('Am_con_Al_no_vacio_no_pol.txt')
    #the format is data_frame, something usual in data science
Am_Al_mal_panda.rename(columns={'0':'cuentas'})
################################################################'''


'''#@@@@@Reading the .txt manually created form the .mca@@@@@@@@@@@@@@22@@
with open('Am_con_Al_no_vacio_no_pol.txt') as file_object:   #file_object stores
            # the object representing the file (pi.txt), which is opened with
            #the open function
            
        Am_Al_mal = file_object.read()
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        Am_Al_mal = Am_Al_mal.split('\n')         #creating a list
                #split turns a string into a list. Between brackets goes the 
                #sepparator.
                
        while '' in Am_Al_mal:   #removal of empty lines
            Am_Al_mal.remove('')    
            
        Am_Al_mal = list(map(int, Am_Al_mal))       #converting it to numbers
        channels = list(range(1, len(Am_Al_mal)+1))     #number of chanels
                #list(range(1,5)) gives [1,2,3,4]        
###########################################################
'''

###################1.1) Representacion###################

    #i) Plot normal
import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
        
plt.plot(channels,Am_Al_mal)

# Set chart title and label axes.
plt.title("Am con capa Al aire no polarización", fontsize=24)           #title
plt.xlabel("Canal", fontsize=14)                        #xlabel
plt.ylabel("Cuentas", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis

plt.close("all")                    #to close all the plots
    
    #ii) Bar plot       
plt.bar(channels,Am, width=1.0)
        #width 1.0 to remove the gaps between bars

# Set chart title and label axes.
plt.title("Am", fontsize=24)          #title
plt.xlabel("Canal", fontsize=14)                                       #xlabel
plt.ylabel("Cuentas", fontsize=14)                                     #ylabel
plt.tick_params(axis='both', labelsize=14)            #size of tick labels  
plt.xlim(0,2500)                                            #limits of x axis
plt.ylim(0,200)                                             #limits of y axis

    #iii) Panda    
#Am_Al_mal_panda.plot()

    #iV) Pygal (what the book crash course uses) 
            #not recommended becasue the x axis has to be introduced manually


################2) Calibración canal-energía##########################

#from sklearn.linear_model import LinearRegression #a way to Linear Regression
    #discarted because it does not give error of the coefficients, just
    #residuals
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)

unidad_rara = [1.59, 2.09, 3.48, 4.97, 5.96, 6.95, 7.94, 8.94] 
        #valores unidades raras de los pulsos (voltaje, o relacionado)
canal_unidad_rara = [649, 854, 1427, 2037, 2445, 2849, 3255, 3658]
    #los canales de los amximos de cuentas para cada pulso

'''
#Reshape to do an instance of Linear Regression
u_r = np.array(unidad_rara).reshape(-1,1)
canal_u_r = np.array(canal_unidad_rara).reshape(-1,1)
    
ajuste_c_u_rara = LinearRegression().fit(canal_u_r,
                                                  u_r) #this is an instance
    #of the class LineaRegression that contains our linear regression
r_square = ajuste_c_u_rara.score(canal_u_r, u_r)           #r^2

#u_r = m * pendiente + ord_or
pendiente = ajuste_c_u_rara.coef_              
ord_origen = ajuste_c_u_rara.intercept_
ajuste_c_u_rara.coef_

#Esta por ver esto, pq da residuos, q no errores
#################
'''

'''OTRA FORMA DE HACER EL AJUSTE, PERO QUE NO DA ERRORES
data = pd.DataFrame({'canal': canal_unidad_rara, 'unidad_rara': unidad_rara, })
            #creacion de un panda con los datos (primero x luego y)
ajuste_u_r_canal_u_r = ols(" unidad_rara ~ canal_unidad_rara", data).fit()  
                #model
                
print(ajuste_u_r_canal_u_r.summary())
'''

#Version con mi funcion (testeada por supuesto)
import RegresionLineal

ajuste = RegresionLineal.RegresionLineal(canal_unidad_rara, unidad_rara)

#############3) Calibracion canal-energía#############################
'''hora le debes meter un valor de enrgía, para variar así b, y obtener ya
energia=a+canal+c, siendo c distinto de b. 
O sea, tu con el pulser le
metes pulsos en el ampli. Barres todo el espectro, dando pulsos en
distintos canales. Asi hallarias la proporcionalidad entre el voltaje de
los pulsos y los canales. Y sabes que el votlaje stá relacionado con la
energia, pues si se deposita una energia dad se crea un pulso de un
voltaje dado.
 O sea, metes ahi las energias y canal de un dato sabido, de un espectro,
 y sacas la cte c: E=a*canal+c==> c=E-a*canal
'''

c = 5.4856 - ajuste['Slope'] * 2212;        #nueva ordenada en origen [MeV]

#Entonces, si tengo el canal, para pasarlo a energía debo hacer:
    #E = slope * canal+c. 
    #La forma rapida de hacerlo es con una comprehensive list:

Energies = [channel *ajuste['Slope'] + c  for channel in channels] #MeV

    
    
#################4) Plot para energías################

plt.bar(Energies,Am, width = 0.0025)
        #width modified to match matlab plots. Default = 0.8

# Set chart title and label axes. (fontsize changed to match the space that 
    #is saved)
plt.title("Am", fontsize=24)          #title
plt.xlabel("E (MeV)", fontsize=12)                                    #xlabel
plt.ylabel("Cuentas", fontsize=12)                                    #ylabel
plt.tick_params(axis='both', labelsize=12)            #size of tick labels  
plt.grid(True)                                              #show grid
plt.xlim(5.35,5.55)                                         #limits of x axis
#plt.ylim(0,200)                                            #limits of y axis
#plt.yscale('log')                                          #y axis in logscale
plt.savefig('espectro_Am.eps', format='eps')
plt.savefig('espectro_Am.png', format='png')
plt.savefig('espectro_Am.pdf', format='pdf')



#############5) Gaussian fit to calc FWHM######################
#Source: http://emilygraceripka.com/blog/16

import math 
#from scipy.stats import norm               ##norm.fit() fit to gaussian
import scipy

#In this section I will carry out a fit of a peak of one spectra


 
#Definition of the function to use to fit the data
def gaussian(x, a, b, c):
    return a * np.exp(- (x-b)**2 / (2 * c**2))   

        #this is a gaussian function (more general than normal distribution)
        
        #if using math.exp the fit gives error: 
        #(only size-1 arrays can be converted to Python scalars)
        #with numpy.exp everything fine.       
        #
        #FRIENDSHIP ENDED WITH math.exp(), NOW numpy.exp() IS MY 
        #BEST FRIEND
        #

#Creation of the array needed to do the fit
x_data = np.array(Energies[2195:2236])      
y_data = np.array(Am[2195:2236])


'''Guesses for the fit (from 
https://stackoverflow.com/questions/19206332/gaussian-fit-for-python)
n = len(x_data)                          #the number of data
mean = sum(x_data*y_data)/n                   #note this correction
sigma = sum(y_data*(x_data-mean)**2)/n        #note this correction

gaussian_fit = scipy.optimize.curve_fit(gaussian, x_data, y_data,
                                        p0=[max(y_data),mean,sigma])

VERY BAD IDEA, DO NOT CONVERGE. wITHOUT GUESSES, IT FINDS THE APPROPIATE FIT!
Have commented on that post that guesses are not needed
'''

gaussian_fit = scipy.optimize.curve_fit(gaussian, x_data, y_data)

opt_values = gaussian_fit[0]   #optimal values of the function to fit the data
cov_of_opt_val = gaussian_fit[1]            #covariances of the optimal values
    #the diagonal are the variance of the parameter to estimate.
    
a = opt_values[0]  
b = opt_values[1]
cc = opt_values[2]  
        #similar values as the ones given by the fit function in Matlab :)

perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error (el 
                                                #error de toa la via vamos)
    #source: 
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html


sigma = cc/np.sqrt(2)                   #standard deviation of the gaussian fit
Delta_sigma = perr[2]/np.sqrt(2)        #error of the standar deviation
print('sigma: ' + str(sigma) + ' +/- ' + str(Delta_sigma) + ' MeV')

FWHM = 2 * np.sqrt(2 * np.exp(2)) * sigma                   #FWHM of the peak
Delta_FWHM = 2 * np.sqrt(2 * np.exp(2)) * Delta_sigma     #error of the FWHM
print('FWHM: ' + str(FWHM) + ' +/- ' + str(Delta_FWHM) + ' MeV')

##Plot
plt.plot(x_data, y_data, label = 'data')        #original data
plt.plot(x_data, gaussian(x_data, a, b, cc), 'ro', label = 'fit')
plt.title("Ajuste Gaussiano a pico de Am", fontsize=24)          #title
plt.xlabel("E (MeV)", fontsize=12)                                    #xlabel
plt.ylabel("Cuentas", fontsize=12)                                    #ylabel
plt.tick_params(axis='both', labelsize=12)            #size of tick labels  
plt.grid(True)                                              #show grid
#plt.xlim(5.35,5.55)                                         #limits of x axis

###

#Bueno, estas son las cosas fundamentales que harás analizando espectros,
#asi q ya me doy por satisfecho xD. He probado una vez cada una.


########################################################################
########################################################################
########################################################################
'''ABOUT SIZE
np.size(Energies[2195:2236])
#arrays:
xx = [1, 2, 3, 4, 5] #horizontal
xxx = np.array([[1],[2],[3]])   #vertical (muestra colorines en el explorer)
'''