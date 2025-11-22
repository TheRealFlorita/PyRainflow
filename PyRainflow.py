# -*- coding: utf-8 -*-
import os
import ctypes

# Load dll and functions
PyRainflowDLL = ctypes.PyDLL(os.path.join(os.path.dirname(__file__),'PyRainflow.dll'))

PyRainflowDLL.PyFilterPeaks.argtypes = (ctypes.py_object, )
PyRainflowDLL.PyFilterPeaks.restype = ctypes.py_object

PyRainflowDLL.PyRainflowCounting.argtypes = (ctypes.py_object, ctypes.c_int, ctypes.c_double, ctypes.c_double)
PyRainflowDLL.PyRainflowCounting.restype = ctypes.py_object

PyRainflowDLL.PyRainflowCountingRandomOrder.argtypes = (ctypes.py_object, ctypes.c_double, ctypes.c_double, ctypes.c_bool, ctypes.c_bool)
PyRainflowDLL.PyRainflowCountingRandomOrder.restype = ctypes.py_object


# Python wrapper functions
def PyFilterPeaks(values):
    if not isinstance(values, list):
        raise TypeError('Expected a list object')
        
    global PyRainflowDLL

    return PyRainflowDLL.PyFilterPeaks(values)


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

    output = PyRainflowDLL.PyRainflowCounting(values, alg_enum, tolerance, cutoff)

    sigma_delta = [item[0] for item in output]
    sigma_mean =  [item[1] for item in output]
    cycle_count = [item[2] for item in output]
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

    output = PyRainflowDLL.PyRainflowCountingRandomOrder(data, tolerance, cutoff, randomise, verify)

    sigma_delta = [item[0] for item in output]
    sigma_mean =  [item[1] for item in output]
    cycle_count = [item[2] for item in output]
    return sigma_delta, sigma_mean, cycle_count