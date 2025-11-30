# -*- coding: utf-8 -*-
import os
import ctypes

# Load dll and functions
PyRainflowDLL = ctypes.PyDLL(os.path.join(os.path.dirname(__file__),'PyRainflow.dll'))

PyRainflowDLL.PyGetPeaks.argtypes = (ctypes.py_object, )
PyRainflowDLL.PyGetPeaks.restype = ctypes.py_object

PyRainflowDLL.PyRainflowCounting.argtypes = (ctypes.py_object, ctypes.c_int, ctypes.c_double, ctypes.c_double)
PyRainflowDLL.PyRainflowCounting.restype = ctypes.py_object

PyRainflowDLL.PyRainflowCountingRandomOrder.argtypes = (ctypes.py_object, ctypes.c_double, ctypes.c_double, ctypes.c_bool, ctypes.c_bool)
PyRainflowDLL.PyRainflowCountingRandomOrder.restype = ctypes.py_object


# Python wrapper functions
def PyGetPeaks(values):
    if not isinstance(values, list):
        raise TypeError('Expected a list object')
        
    global PyRainflowDLL

    return PyRainflowDLL.PyGetPeaks(values)


def PyRainflowCounting(values, algorithm='4-points', tolerance=1e-3, cutoff=1e-1):
    if not isinstance(values, list):
        raise TypeError('Expected a list object')
    elif len(values) < 2:
        return [],[],[]
    elif not isinstance(values[0], float):
        raise TypeError('Expected a list of floats')

    global PyRainflowDLL

    alg_enum = 0
    if algorithm.lower() == '3-points, non-periodic':
        algealg_enumnum = 1
    elif algorithm.lower() == '3-points':
        alg_enum = 2

    sigmaDelta, sigmaMean, cycleCount = PyRainflowDLL.PyRainflowCounting(values, alg_enum, tolerance, cutoff)

    return sigma_delta, sigma_mean, cycle_count


def PyRainflowCountingRandomOrder(data, tolerance=1e-3, cutoff=1e-1):
    if not isinstance(data, list):
        raise TypeError('Expected a list object')
    elif len(data) == 0:
        return [],[],[]
    elif not isinstance(data[0], list):
        raise TypeError('Expected a list of lists -> [[int, list], [int, list], ...]')
    elif not isinstance(data[0][1][0], float):
        raise TypeError('Expected a list of floats')
    
    global PyRainflowDLL
    randomise = True
    verify = False

    sigmaDelta, sigmaMean, cycleCount = PyRainflowDLL.PyRainflowCountingRandomOrder(data, tolerance, cutoff, randomise, verify)

    return sigma_delta, sigma_mean, cycle_count