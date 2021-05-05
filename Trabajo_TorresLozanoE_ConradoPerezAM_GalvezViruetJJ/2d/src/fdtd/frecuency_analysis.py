import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.fft import rfft, fftfreq
from fdtd.common import X, Y, L, U

def Fourier_trans(measures,data,data_input):
    """
    Calculates the Fourier transform of input data, and the reflection and transmission coefficients, total and 
    as a function of frequency.
    | Inputs:
    | - measures: Object with ports' powers.
    | - data: Object with simulation times.
    | - data_input: Object with ports' positions.
    | Output, touple with three objects, contains:
    | - ports: Contains information about ports: distance to the source, time of signals and measures in those times,
    | frequencies and Fourier transform.
    | - [R,T]: A list with the reflexion and transmission coefficients as functions of frecuency.
    | - [R_t,T_t]: Total reflection and transmission coefficients.

    """
    # Power values
    measures = copy.deepcopy(measures)
    
    # Times
    times = np.array(list(map(lambda x: x*10**9, copy.deepcopy(data[0]["time"]))))
    Ntimes = len(times)

    # Ports positions
    sourceposition = data_input["coordinates"][data_input["elements"][data_input["sources"][0]["elemId"]][0]][0]
    ports = {"0": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_inc"]["elemId"]][0]][0])],\
        "1": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_refl"]["elemId"]][0]][0])],\
        "2": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_trans"]["elemId"]][0]][0])],}
    
    # Port times
    ports["0"].append([0,15])    # Firsts times at port 0 in ns
    ports["1"].append([18,33])   # Second times at port 1 in ns
    ports["2"].append([20,35])   # Firsts times at port 3 in ns

    for i in ports:
        ports[i].append([[],[]]) 


    for i in ports:
        for j in range(0,Ntimes): 
            if (ports[i][1][0] < times[j] and ports[i][1][1] > times[j]):
                ports[i][2][1].append(abs(measures.Ports(int(i))[j]))
                ports[i][2][0].append(times[j])
            else:
                ports[i][2][1].append(0)
                ports[i][2][0].append(times[j])

    # Nmeasures = len(ports[i][2][0])
    Nmeasures = Ntimes

    # Transmission and reflection coefficients
    R_t = sum(ports["2"][2][1])/sum(ports["0"][2][1])
    T_t = sum(ports["1"][2][1])/sum(ports["0"][2][1])

    # Frecuencies and Fourier transform.
    timestep = times[1] - times[0]
    for i in ports.values():
        frequencies = (fftfreq(Nmeasures)/timestep)[:Nmeasures//2]
        transform = np.abs((rfft(i[2][1])))[1:]
        i.append([frequencies, transform] )

    # Coefficients as functions of frequency.
    R = ports["2"][3][1]/ports["0"][3][1]
    T = ports["1"][3][1]/ports["0"][3][1]

    return (ports,[R,T],[R_t,T_t])

