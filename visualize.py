# alternatively run from command line and use sys.argv
import csv, math, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

PATH = ""

# R and T flags
STATE = 1

wvlVal = 940
wvlL, wvlH, wvlS = wvlVal - 2, wvlVal + 2, 1/100

wvlTotal = int((wvlH - wvlL) / wvlS + 1)
tL, tH = 0.1, 0.5
tTotal = int((tH - tL) * 1000 + 1)

raL, raH, raS = 0.01, 0.3, 0.01
raTotal = round((raH - raL)/raS) + 1

tList = []
raList = []
paramList = []

# compile data from csv to complex array (for incident and reference simulations)
def complex_compile(DATA, arr):
    
    with open(DATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(len(data)):
            arr[i] = complex(float(data[i][1]), float(data[i][2]))

    return

# compile data from csv to array
def compile(DATA, arr):
    
    with open(DATA, newline = '') as csvfile:
        data = np.array(list(csv.reader(csvfile)))
        for i in range(len(data)):
            arr[i] = float(data[i][0])
    
    return

# find viable r/a and a values (given state--R is 1 and T is 0)
def raSim(wvlTotal, t, ra, a, r, state, paramList, DATA):

    if state: filename1 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

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

# find good phase shift parameters
def phaseSim(t, ra, a, state, paramList, DATA):

    if state: filename1 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

    # compile data
    vals = np.zeros(wvlTotal); compile(filename1, vals)
    reference = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, 'reference.txt'), reference)
    incident = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, 'incident.txt'), incident)
    z0 = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, filename2), z0)
    reflected = z0 - incident

    ratioR = (reference/reflected).real
    ratioI = (reference/reflected).imag
    phase = np.zeros(ratioR.size)
    for i in range(phase.size): phase[i] = math.atan2(ratioI[i], ratioR[i]) + 2 * math.pi if ratioI[i] < 0 else math.atan2(ratioI[i], ratioR[i] )
    phase /= math.pi

    # append to paramList
    for i in range(phase.size): 
        if round(phase[i], 1) == 2.0: paramList.append([t, ra, a]); return

    return

def phase(t, ra, a, state, DATA):

    if state: filename1 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if not state: filename1 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_power.txt")
    if state: filename2 = os.path.join(DATA, f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")
    if not state: filename2 = os.path.join(DATA, f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}_field.txt")

    # compile data
    vals = np.zeros(wvlTotal); compile(filename1, vals)
    reference = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, 'reference.txt'), reference)
    incident = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, 'incident.txt'), incident)
    z0 = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(DATA, filename2), z0)
    reflected = z0 - incident

    ratioR = (reference/reflected).real
    ratioI = (reference/reflected).imag
    phase = np.zeros(ratioR.size)
    for i in range(phase.size): phase[i] = math.atan2(ratioI[i], ratioR[i]) + 2 * math.pi if ratioI[i] < 0 else math.atan2(ratioI[i], ratioR[i] )
    phase /= math.pi

    return phase

def main():

    global wvlVal, wvlL, wvlH, wvlS, wvlTotal, STATE

    PATH = os.getcwd()
    DATA = os.path.join(PATH, str(input("data folder name: ")))
    RDATA = os.path.join(DATA, "raw_data")
    wvlVal = int(input("wavelength value (in nm): "))
    STATE = int(input("run visualization for reflection data or transmission data? (1 for Reflection, 0 for Transmission): "))

    # default wvl values (nm)
    wvlL, wvlH, wvlS = wvlVal - 2, wvlVal + 2, 1/100
    wvlTotal = int((wvlH - wvlL) / wvlS + 1)

    # set default a range values (somewhat arbitrary)
    aL, aH, aS = 0.5, round((wvlVal - 300) / 1000, 3), 0.001
    aTotal = int((aH - aL) / aS)
    
    print("Starting Incident and Reference Compilation...\n")

    # run incident and reference simulations + unload data
    incident = np.zeros(wvlTotal, dtype=np.complex_)
    reference = np.zeros(wvlTotal, dtype=np.complex_)
    complex_compile(os.path.join(RDATA, "incident.txt"), incident)
    complex_compile(os.path.join(RDATA, "reference.txt"), reference)

    print("Starting Thickness Compilation...\n")

    # run t value simulation + unload data
    reflection = np.zeros(tTotal)
    transmission = np.zeros(tTotal)
    compile(os.path.join(RDATA, "tR.txt"), reflection)
    compile(os.path.join(RDATA, "tT.txt"), transmission)

    # find peak values for t
    global tList
    if STATE: peaks, _ = find_peaks(reflection, prominence = 0.1)
    else: peaks, _ = find_peaks(transmission, prominence = 0.1)
    tList.append([tL + peaks[i] / 1000 for i in range(peaks.size)])
    tList = np.ndarray.tolist(np.array(np.squeeze(np.array(tList)), ndmin = 1))

    print("Starting Radius and Lattice Constant Compilation\n")

    for t in tList:
        t = round(t, 3)

        for ra in np.linspace(raL, raH, raTotal):
            ra = round(ra, 3)

            for a in np.linspace(aL, aH, aTotal):
                a = round(a, 3)
                r = round(ra * a, 3)

                raSim(wvlTotal, t, ra, a, r, STATE, raList, RDATA)

    print("Checking Phase Values…\n")

    # calculate and check phase shift values for possible parameters
    for t, ra, a in raList: phaseSim(t, ra, a, STATE, paramList, RDATA)

    print("Graphing Selected Parameters…\n")

    # graph the selected parameters for user
    for t, ra, a in paramList:

        filename = ""
        if STATE: filename = f"R_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}" + "_power.txt"
        else: filename = f"T_t{int(round(t, 3) * 1000)}_ra" + "{:.2f}".format(round(ra, 3)) + f"_a{int(round(a, 3) * 1000)}" + "_power.txt"

        # compile values
        rVals = np.zeros(wvlTotal)
        with open(os.path.join(RDATA, filename), newline = '') as csvfile:
            data = np.array(list(csv.reader(csvfile)))
            for i in range(len(data)):
                rVals[i] = float(data[i][0])
        fderiv = np.diff(rVals); fderiv /= wvlS
        sderiv = np.diff(fderiv); sderiv /= wvlS
        ph = phase(t, ra, a, STATE, DATA)

        fig, ax3 = plt.subplots()
        ax3.plot(np.linspace(wvlL, wvlH, wvlTotal), rVals)
        ax3.set_title(f"original: t = {int(round(t, 3) * 1000)}, ra = {round(ra, 3)}, a = {int(round(a, 3) * 1000)}")
        ax3.set_xlabel("Wavelength (nm)")
        ax3.set_ylabel("Transmittance")
        ax3.set_ylim([0, 2])
        ax = ax3.twinx()
        ax.plot(np.linspace(wvlL, wvlH, wvlTotal), ph, 'y:')
        ax.set_ylabel("Phase Shift")
        ax.set_ylim([0, 2])

        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()