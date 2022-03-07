import csv
import numpy as np

# returns line 1 of a file in integer format
def line1(filename):

    with open(filename, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        return int(data[0][0])

# compile data from csv to array (skipping line 1)
def compile(DATA, arr):
    
    with open(DATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(1, len(data)):
            arr[i - 1] = float(data[i][0])
    return

# compile data from csv to complex array
def complex_compile(filename, arr):
    
    with open(filename, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(1, len(data)):
            arr[i - 1] = complex(float(data[i][1]), float(data[i][2]))

    return

# compile parameter data from csv to list
def param_compile(filename, arr):

    with open(filename, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(1, len(data)):
            arr.append([float(data[i][0]), float(data[i][1]), float(data[i][2])])

# old compile method, used for quick testing when compile() doesn't work (TODO: remove)
def oldcompile(DATA, arr):
    
    with open(DATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(len(data)):
            arr[i] = float(data[i][0])
    return