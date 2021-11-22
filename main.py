import S4 as S4
import sys, csv, math, os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from tqdm import tqdm
from pathlib import Path

PATH = ""

# flags for Reflection and Transmission simulations
R = 1
T = 1

# thickness of air above layer
tBuff = 10.00
tPath = 50.00

# random initial number just for random initialization purposes
initial = 1.00

# default wvl values (wavelength in nm)
wvlVal = 1550
wvlL, wvlH, wvlS = wvlVal - 2, wvlVal + 2, 1/100
wvlTotal = int((wvlH - wvlL) / wvlS + 1)

# default t range values (thickness in microns)
tL, tH = 0.1, 0.5
tTotal = int((tH - tL) * 1000 + 1)

# default r/a range values
raL, raH, raS = 0.01, 0.3, 0.01
raTotal = round((raH - raL)/raS) + 1

# set default a range values (somewhat arbitrary)
aL, aH, aS = 0.5, round((wvlVal - 300) / 1000, 3), 0.001
aTotal = int((aH - aL) / aS)

# selected t values for Reflection and Transmission (t is thickness)
tListR = []
tListT = []

# selected r/a ratio and a values for Reflection and Transmission (r is radius, a is lattice constant)
raListR = []
raListT = []

# selected t, r, a parameters for Reflection and Transmission
paramListR = []
paramListT = []

# create S4 simulation object and initializes everything
def init(mat, eps, a, num = 50):
    S = S4.New(Lattice = ((a, 0), (0, a)), NumBasis = num)

    S.AddMaterial(Name = str(mat), Epsilon = eps)
    S.AddMaterial(Name = "SiO2", Epsilon = 2.1904 + 0j)
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)

    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    S.AddLayer(Name = str(mat), Thickness = initial, Material = str(mat))
    # S.SetRegionCircle(Layer = str(mat), Material = str(mat), Center = (0, 0), Radius = 0) # initialize
    S.AddLayer(Name = 'SOGBelow', Thickness = 0.00, Material = 'SiO2')

    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    return S

# set radius and thickness for a simulation object
def config(S, mat, r, t):

    S.SetLayerThickness(Layer = str(mat), Thickness = t)
    S.SetRegionCircle(Layer = str(mat), Material = 'Vacuum', Center = (0, 0), Radius = r)

    return

# compile data from csv to complex array (for incident and reference simulations)
def complex_compile(RDATA, arr):
    
    with open(RDATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(len(data)):
            arr[i] = complex(float(data[i][1]), float(data[i][2]))

    return

# run incident simulation
def incSim(wvlL, wvlH, wvlTotal, RDATA):

    output = open(os.path.join(RDATA, 'incident.txt'), 'w+')

    S = S4.New(Lattice = ((1, 0), (0, 1)), NumBasis = 50)
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)

    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    S.AddLayer(Name = 'SOGBelow', Thickness = 0.00, Material = 'Vacuum')
    
    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    for wvl in np.linspace(wvlL, wvlH, wvlTotal):
        freq = 1000 / wvl
        S.SetFrequency(freq)

        (E, _) = S.GetFields(0, 0, tBuff)
        output.write(str(round(wvl, 1)) + ',' + str(E[1].real) + ',' + str(E[1].imag) + '\n')

        print("incident", round(wvl, 1), '\n', end = " ")
        sys.stdout.flush()

    output.close()

    return

# run reference simulation
def refSim(wvlL, wvlH, wvlTotal, RDATA):

    output = open(os.path.join(RDATA, 'reference.txt'), 'w+')

    S = S4.New(Lattice = ((1, 0), (0, 1)), NumBasis = 50)
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)
        
    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    S.AddLayer(Name = 'SOGBelow', Thickness = 0.00, Material = 'Vacuum')
    
    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    for wvl in np.linspace(wvlL, wvlH, wvlTotal):
        freq = 1000 / wvl
        S.SetFrequency(freq)

        (E, _) = S.GetFields(0, 0, tBuff + 2 * tPath)
        output.write(str(round(wvl, 1)) + ',' + str(E[1].real) + ',' + str(E[1].imag) + '\n')

        print("reference", round(wvl, 1), '\n', end = " ")
        sys.stdout.flush()

    output.close()

    return

# compile data from csv to array
def compile(DATA, arr):
    
    with open(DATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(len(data)):
            arr[i] = float(data[i][0])
    
    return

# find viable t values (for both R and T)
def tSim(mat, eps, wvlVal, tL, tH, tTotal, RDATA):

    S = init(mat, eps, initial)
    filenameR = os.path.join(RDATA, "tR.txt")
    filenameT = os.path.join(RDATA, "tT.txt")

    out_tR = open(filenameR, "w+")
    out_tT = open(filenameT, "w+")
    for t in np.linspace(tL, tH, tTotal):
        t = round(t, 3)

        S.SetLayerThickness(Layer = str(mat), Thickness = t)

        freq = 1000 / wvlVal
        S.SetFrequency(freq)

        _, pyntbacktop = S.GetPowerFlux('AirBuff', 0)
        pyntfowbot, _ = S.GetPowerFlux('SOGBelow', 0)

        out_tR.write(str(-pyntbacktop.real) + '\n')
        out_tT.write(str(pyntfowbot.real) + '\n')

        print("thickness", t, '\n', end = " ")
        sys.stdout.flush()
    out_tR.close()
    out_tT.close()

    reflection = np.zeros(tTotal)
    transmission = np.zeros(tTotal)
    compile(filenameR, reflection)
    compile(filenameT, transmission)

    print("Finding peak values for t...\n")

    # find peak values for t
    rPeaks, _ = find_peaks(reflection, prominence = 0.1)
    tPeaks, _ = find_peaks(transmission, prominence = 0.1)

    print("Updating list of potential t values...\n")

    # update list with thickness values to try
    global tListR, tListT
    tListR.append([tL + rPeaks[i] / 1000 for i in range(rPeaks.size)])
    tListT.append([tL + tPeaks[i] / 1000 for i in range(tPeaks.size)])
    tListR = np.ndarray.tolist(np.array(np.squeeze(np.array(tListR)), ndmin = 1)) # remove extra (useless) dimension
    tListT = np.ndarray.tolist(np.array(np.squeeze(np.array(tListT)), ndmin = 1))

    return

# find viable r/a and a values ("state" parameter is a flag--1 for Reflection and 0 for Transmission)
def raSim(mat, eps, wvlL, wvlH, wvlTotal, t, ra, a, r, state, paramList, RDATA):

    if state: filename1 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

    S = init(mat, eps, a)
    out_ra_power = open(filename1, "w+")
    out_ra_field = open(filename2, "w+")
    for wvl in np.linspace(wvlL, wvlH, wvlTotal):
        S.SetRegionCircle(Layer = str(mat), Material = 'Vacuum', Center = (0, 0), Radius = r)
        S.SetLayerThickness(Layer = str(mat), Thickness = t)

        freq = 1000 / wvl
        S.SetFrequency(freq)

        _, pyntbacktop = S.GetPowerFlux('AirBuff', 0)
        pyntfowtop, _ = S.GetPowerFlux('SOGBelow', 0)
        
        if state: out_ra_power.write(str(-pyntbacktop.real) + '\n')
        if not state: out_ra_power.write(str(pyntfowtop.real) + '\n')

        (E, _) = S.GetFields(0.5 * a, 0.5 * a, tBuff)
        out_ra_field.write(str(round(wvl, 1)) + ',' + str(E[1].real) + ',' + str(E[1].imag) + '\n')
        
        print("thickness: ", round(t, 3), "\tra: ", round(ra, 3), "\ta: ", round(a, 3), "\twvl: ", round(wvl, 3), end = "\n")
        sys.stdout.flush()
    out_ra_power.close()
    out_ra_field.close()

    raVals = np.zeros(wvlTotal)            
    compile(filename1, raVals)
    fderiv = np.diff(raVals); fderiv /= wvlS # first derivative
    sderiv = np.diff(fderiv); sderiv /= wvlS # second derivative

    # find peak points (gives indices of those points)
    fderivPeaks, _ = find_peaks(fderiv, prominence = 0.1)
    fderivTroughs, _ = find_peaks(-fderiv, prominence = 0.1)
    sderivPeaks, _ = find_peaks(sderiv, prominence = 0.1)
    sderivTroughs, _ = find_peaks(-sderiv, prominence = 0.1)

    if fderiv[fderivPeaks].size: # if it has notable features (not just a flat/curved line or anything)
            if (fderivPeaks.size == 1) and (fderivTroughs.size == 1) and (sderivPeaks.size == 2) and (sderivTroughs.size == 1): # fits conditions for "good" mode part 1 (based on observations made on first/second derivative graphs)
                if (fderivPeaks[0] < fderivTroughs[0]) and (sderivPeaks[0] < sderivTroughs[0]) and (sderivTroughs[0] < sderivPeaks[1]): # verify shape is correct part 2
                    paramList.append([t, ra, a])

    return

# find good phase shift parameters (see raSim for state parameter info)
def phaseSim(t, ra, a, state, paramList, RDATA):

    if state: filename1 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

    # compile data
    vals = np.zeros(wvlTotal); compile(filename1, vals)
    reference = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, 'reference.txt'), reference)
    incident = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, 'incident.txt'), incident)
    z0 = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, filename2), z0)
    reflected = z0 - incident
    # phase = 

    ratioR = (reference/reflected).real
    ratioI = (reference/reflected).imag
    phase = np.zeros(ratioR.size)
    for i in range(phase.size): phase[i] = math.atan2(ratioI[i], ratioR[i]) + 2 * math.pi if ratioI[i] < 0 else math.atan2(ratioI[i], ratioR[i] )
    phase /= math.pi

    # append to paramList
    for i in range(phase.size): 
        if round(phase[i], 1) == 2.0: paramList.append([t, ra, a]); return

    return

def phase(t, ra, a, state, RDATA):

    if state: filename1 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(RDATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(RDATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

    # compile data
    vals = np.zeros(wvlTotal); compile(filename1, vals)
    reference = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, 'reference.txt'), reference)
    incident = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, 'incident.txt'), incident)
    z0 = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, filename2), z0)
    reflected = z0 - incident

    ratioR = (reference/reflected).real
    ratioI = (reference/reflected).imag
    phase = np.zeros(ratioR.size)
    for i in range(phase.size): phase[i] = math.atan2(ratioI[i], ratioR[i]) + 2 * math.pi if ratioI[i] < 0 else math.atan2(ratioI[i], ratioR[i] )
    phase /= math.pi

    return phase

def main():

    # store time taken info
    Path('logs').mkdir(parents = True, exist_ok = True)
    time_logs = open('logs/time_logs.txt', 'w+')
    start_time = time.perf_counter()

    PATH = os.getcwd() # automatically grab current working directory

    # get input data
    DATA = os.path.join(PATH, str(input("Input folder name (data files will be saved in this folder): ")))
    if not os.path.isdir(DATA): os.mkdir(DATA) # create data directory
    print(f"\nData files will be saved in : '{DATA}'\n\n")

    # raw data folder indside data folder to reduce clutter
    RDATA = os.path.join(DATA, "raw_data")
    if not os.path.isdir(RDATA): os.mkdir(RDATA)

    wvlVal = int(input("Input target wavelength in nm: "))
    print(f"\nTarget wavelength is {wvlVal} nm\n")

    R = int(input("Run Reflection simulation? (1 for Yes, 0 for No): "))
    T = int(input("Run Transmission simulation? (1 for Yes, 0 for No): "))

    mat = str(input("Input material name: ")) # TODO: change to a list for multilayer
    print(f"\nMaterial name is '{mat}'\n\n")

    # TODO: possibly prompt for layer name

    re = float(input("Input epsilon real component: "))
    im = float(input("Input epsilon imaginary component: "))
    eps = complex(re, im)
    print(f"\nMaterial epsilon is {eps}\n\n")

    # print("Input n value (refractive index of the PC layer material): ")
    # n = float(input())
    # print(f"n value is {n}")

    print("Input thickness range in microns (use default recommended values if no preference):")
    tL = float(input("Lower bound (ex. 0.1): "))
    tH = float(input("Upper bound (ex. 0.5): "))
    tTotal = int((tH - tL) * 1000 + 1)
    print(f"\nLower bound is {tL} microns and upper bound is {tH} microns, for a total of {tTotal} iterations\n\n")
    time.sleep(3)

    """
    # default values used for testing:
    DATA = os.path.join(PATH, "data")
    if not os.path.isdir(DATA): os.mkdir(DATA)
    RDATA = os.path.join(DATA, "raw_data")
    if not os.path.isdir(RDATA): os.mkdir(RDATA)
    mat = 'Silicon'
    re = 12.1104
    im = 0
    eps = complex(re, im)
    """

    # update wvl values
    global wvlL, wvlH, wvlS, wvlTotal
    wvlL, wvlH, wvlS = wvlVal - 2, wvlVal + 2, 1/100
    wvlTotal = int((wvlH - wvlL) / wvlS + 1)

    # update a values
    global aL, aH, aS, aTotal
    aL, aH, aS = 0.5, round((wvlVal - 300) / 1000, 3), 0.001
    aTotal = int((aH - aL) / aS)

    # create parameter output file
    R_selected_parameters = open(os.path.join(DATA, "R_selected_parameters.txt"), "w+")
    T_selected_parameters = open(os.path.join(DATA, "T_selected_parameters.txt"), "w+")

    # update time logs
    curr_time = time.perf_counter()
    time_logs.write(f"input + intial setup took {int((curr_time - start_time) / 3600)} hours, {int(((curr_time - start_time) % 3600) / 60)} minutes, and {int(((curr_time - start_time) % 3600) % 60)} seconds" + '\n')
    prev_time = curr_time
    
    print("Starting Incident and Reference Simulations...\n")

    # run incident and reference simulations + unload data
    incSim(wvlL, wvlH, wvlTotal, RDATA)
    refSim(wvlL, wvlH, wvlTotal, RDATA)
    incident = np.zeros(wvlTotal, dtype=np.complex_)
    reference = np.zeros(wvlTotal, dtype=np.complex_)
    complex_compile(os.path.join(RDATA, "incident.txt"), incident)
    complex_compile(os.path.join(RDATA, "reference.txt"), reference)

    # update time logs
    curr_time = time.perf_counter()
    time_logs.write(f"incident + reference simulations took {int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds" + '\n')
    prev_time = curr_time

    print("Starting Thickness Simulation...\n")

    # run t value simulation + unload data
    tSim(mat, eps, wvlVal, tL, tH, tTotal, RDATA)
    reflection = np.zeros(tTotal)
    transmission = np.zeros(tTotal)
    compile(os.path.join(RDATA, "tR.txt"), reflection)
    compile(os.path.join(RDATA, "tT.txt"), transmission)

    # update time logs
    curr_time = time.perf_counter()
    time_logs.write(f"thickness simulation took {int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds" + '\n')
    prev_time = curr_time

    print("Starting Radius and Lattice Constant Simulations...\n")

    if R:
        # run r/a ratio + a value simulation for Reflection
        for t in tListR:
            t = round(t, 3)

            for ra in np.linspace(raL, raH, raTotal):
                ra = round(ra, 3)

                for a in np.linspace(aL, aH, aTotal):
                    a = round(a, 3)
                    r = round(ra * a, 3)

                    raSim(mat, eps, wvlL, wvlH, wvlTotal, t, ra, a, r, 1, raListR, RDATA)

        curr_time = time.perf_counter()
        time_logs.write(f"ra + a simulation (Reflection) took {int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds" + '\n')
        prev_time = curr_time

    if T:
        # run r/a ratio + a value simulation for Transmission
        for t in tListT:
            t = round(t, 3)

            for ra in np.linspace(raL, raH, raTotal):
                ra = round(ra, 3)

                for a in np.linspace(aL, aH, aTotal):
                    a = round(a, 3)
                    r = round(ra * a, 3)

                    raSim(mat, eps, wvlL, wvlH, wvlTotal, t, ra, a, r, 0, raListT, DATA)

        curr_time = time.perf_counter()
        time_logs.write(f"ra + a simulation (Transmission) took {int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds" + '\n')
        prev_time = curr_time

    print("Checking Phase Values...\n")

    if R:
        # calculate and check phase shift values for possible parameters (for Reflection)
        for t, ra, a in raListR: phaseSim(t, ra, a, 1, paramListR, RDATA)
    
    if T:
        # calculate and check phase shift values for possible parameters (for Transmission)
        for t, ra, a in raListT: phaseSim(t, ra, a, 0, paramListT, RDATA)
    
    # update time logs
    curr_time = time.perf_counter()
    time_logs.write(f"phase simulation/calculation took {int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds" + '\n')
    prev_time = curr_time

    print("Saving Selected Parameters...\n")

    if R:
        # save output parameters
        R_selected_parameters.write(f"Selected parameters (Reflection, {wvlVal} nm): \n")
        for t, ra, a in paramListR:
            R_selected_parameters.write(f"t: {t}\tra: {ra}\ta: {a}\n")
        
        if len(paramListR) == 0: R_selected_parameters.write("No parameter combinations found with 2pi phase shift.\n")
        
        R_selected_parameters.close()
    
    if T:
        # save output parameters
        T_selected_parameters.write(f"Selected parameters (Transmission, {wvlVal} nm): \n")
        for t, ra, a in paramListT:
            T_selected_parameters.write(f"t: {t}\tra: {ra}\ta: {a}\n")
        
        if len(paramListT) == 0: T_selected_parameters.write("No parameter combinations found with 2pi phase shift.\n")
        
        T_selected_parameters.close()

    print("Graphing Selected Parameters...\n")

    if R:

        # graph the selected parameters for user
        for t, ra, a in paramListR:

            # compile values
            rVals = np.zeros(wvlTotal)
            with open(os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}" + "_power.txt"), newline = '') as csvfile:
                data = np.array(list(csv.reader(csvfile)))
                for i in range(len(data)):
                    rVals[i] = float(data[i][0])
            fderiv = np.diff(rVals); fderiv /= wvlS
            sderiv = np.diff(fderiv); sderiv /= wvlS
            ph = phase(t, ra, a, 1, RDATA)

            # # peak/trough coords
            # fderivPeaks, _ = find_peaks(fderiv, prominence = 0.1)
            # fderivTroughs, _ = find_peaks(-fderiv, prominence = 0.1)
            # sderivPeaks, _ = find_peaks(sderiv, prominence = 0.1)
            # sderivTroughs, _ = find_peaks(-sderiv, prominence = 0.1)
            
            # # second derivative
            # fig, ax1 = plt.subplots()
            # ax1.plot(np.linspace(wvlL, wvlH - 2 * wvlS, wvlTotal - 2), sderiv)
            # ax1.plot(np.round(wvlL + sderivPeaks * wvlS, 3), sderiv[sderivPeaks], 'xr')
            # ax1.plot(np.round(wvlL + sderivTroughs * wvlS, 3), sderiv[sderivTroughs], 'xr')
            # ax1.set_xlabel("Wavelength (nm)")
            # ax1.set_ylabel("Reflectance")
            # ax1.set_title(f"second: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")

            # # first derivative
            # fig, ax2 = plt.subplots()
            # ax2.plot(np.linspace(wvlL, wvlH - wvlS, wvlTotal - 1), fderiv)
            # ax2.plot(np.round(wvlL + fderivPeaks * wvlS, 3), fderiv[fderivPeaks], 'xr')
            # ax2.plot(np.round(wvlL + fderivTroughs * wvlS, 3), fderiv[fderivTroughs], 'xr')
            # ax2.set_xlabel("Wavelength (nm)")
            # ax2.set_ylabel("Reflectance")
            # ax2.set_title(f"first: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")

            # original
            fig, ax3 = plt.subplots()
            ax3.plot(np.linspace(wvlL, wvlH, wvlTotal), rVals)
            ax3.set_title(f"R: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")
            ax3.set_xlabel("Wavelength (nm) [solid blue]")
            ax3.set_ylabel("Reflectance")
            ax3.set_ylim([0, 2])
            ax = ax3.twinx()
            ax.plot(np.linspace(wvlL, wvlH, wvlTotal), ph, 'y:')
            ax.set_ylabel("Phase Shift (π) [dotted yellow]")
            ax.set_ylim([0, 2])

            plt.tight_layout()
            plt.show()

    if T:

        # graph the selected parameters for user
        for t, ra, a in paramListT:

            # compile values
            tVals = np.zeros(wvlTotal)
            with open(os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}" + "_power.txt"), newline = '') as csvfile:
                data = np.array(list(csv.reader(csvfile)))
                for i in range(len(data)):
                    tVals[i] = float(data[i][0])
            fderiv = np.diff(tVals); fderiv /= wvlS
            sderiv = np.diff(fderiv); sderiv /= wvlS
            ph = phase(t, ra, a, 1, RDATA)

            # # peak/trough coords
            # fderivPeaks, _ = find_peaks(fderiv, prominence = 0.1)
            # fderivTroughs, _ = find_peaks(-fderiv, prominence = 0.1)
            # sderivPeaks, _ = find_peaks(sderiv, prominence = 0.1)
            # sderivTroughs, _ = find_peaks(-sderiv, prominence = 0.1)
            
            # # second derivative
            # fig, ax1 = plt.subplots()
            # ax1.plot(np.linspace(wvlL, wvlH - 2 * wvlS, wvlTotal - 2), sderiv)
            # ax1.plot(np.round(wvlL + sderivPeaks * wvlS, 3), sderiv[sderivPeaks], 'xr')
            # ax1.plot(np.round(wvlL + sderivTroughs * wvlS, 3), sderiv[sderivTroughs], 'xr')
            # ax1.set_xlabel("Wavelength (nm)")
            # ax1.set_ylabel("Transmittance")
            # ax1.set_title(f"second: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")

            # # first derivative
            # fig, ax2 = plt.subplots()
            # ax2.plot(np.linspace(wvlL, wvlH - wvlS, wvlTotal - 1), fderiv)
            # ax2.plot(np.round(wvlL + fderivPeaks * wvlS, 3), fderiv[fderivPeaks], 'xr')
            # ax2.plot(np.round(wvlL + fderivTroughs * wvlS, 3), fderiv[fderivTroughs], 'xr')
            # ax2.set_xlabel("Wavelength (nm)")
            # ax2.set_ylabel("Transmittance")
            # ax2.set_title(f"first: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")

            # original
            fig, ax3 = plt.subplots()
            ax3.plot(np.linspace(wvlL, wvlH, wvlTotal), tVals)
            ax3.set_title(f"T: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")
            ax3.set_xlabel("Wavelength (nm) [solid blue]")
            ax3.set_ylabel("Transmittance")
            ax3.set_ylim([0, 2])
            ax = ax3.twinx()
            ax.plot(np.linspace(wvlL, wvlH, wvlTotal), ph, 'y:')
            ax.set_ylabel("Phase Shift (π) [dotted yellow]")
            ax.set_ylim([0, 2])

            plt.tight_layout()
            plt.show()

    # save final time taken for whole simulation process
    curr_time = time.perf_counter()
    time_logs.write(f"total time taken: {round(curr_time - start_time, 3)} seconds ({int((curr_time - prev_time) / 3600)} hours, {int(((curr_time - prev_time) % 3600) / 60)} minutes, and {int(((curr_time - prev_time) % 3600) % 60)} seconds)")
    time_logs.close()

if __name__ == "__main__":
    main()