from dataclasses import fields
from aux_functions.sim import *
from aux_functions.calc import *
from scipy.signal import find_peaks
import math
import matplotlib.pyplot as plt

# R3calc
F_PEAK_PROM = 0.08 # fderiv peak prominence
F_TROUGH_PROM = 0.08 # fderiv trough prominence
VAL_PEAK_PROM = 0.03 # R/T peak prominence
F_PAIR_X_DIST = 5 # max x-dist between peak/trough pair
F_PAIR_Y_DIST = 5 # max y-diff between peak/trough pair
F_NEIGHBOR_X_DIST = 8 # min x-dist between pair/another spike
FLAT_LENGTH = 20 # flat section x-dist
FLAT_ALIGNMENT = 0.3 # median y-diff between flat diffs
MIN_AMPLITUDE = 0.8 # min amplitude (% R/T)

# # scanning range of r/a + amplitude filter, must be less than or equal to WINDOW_RANGE (see main.py)
# SCANNING_RANGE = 0 # note: variable is just initialized here; if you want to change (add) this feature, see first few lines of R3() below

# selected [t, r/a, a]
paramList = []

def R3calc1(wvlS, wvlTotal, filename):

    # t = int(filename[-27:-23])
    # ra = float(filename[-20:-15])
    # a = int(filename[-13:-10])

    # fderiv: shape (spike pattern)
    vals = np.zeros(wvlTotal)            
    compile(filename, vals)
    fderiv = np.diff(vals); fderiv /= wvlS # first derivative

    # spike indices
    fderivPeaks, _ = find_peaks(fderiv, prominence = F_PEAK_PROM)
    fderivTroughs, _ = find_peaks(-fderiv, prominence = F_TROUGH_PROM)
    iPairs = [] # selected pair indices
    to_delete = [] # pairs to delete

    if not (fderivPeaks.size and fderivTroughs.size): return # no peaks and/or no troughs

    # pair peaks/troughs
    i2 = 0
    for i1 in range(fderivPeaks.size):

        if i2 > fderivTroughs.size: break

        while i2 < fderivTroughs.size and fderivTroughs[i2] < fderivPeaks[i1]: i2 += 1
        if i2 < fderivTroughs.size: iPairs.append([i1, i2])
    
    for iPeak, iTrough in iPairs:

        # make sure peak is close to trough and each is about the same distance away from zero (somewhat symmetrical shape)
        if not (fderivTroughs[iTrough] - fderivPeaks[iPeak] < F_PAIR_X_DIST and abs(fderiv[fderivPeaks[iPeak]]) - abs(fderiv[fderivTroughs[iTrough]]) < F_PAIR_Y_DIST): # approximate cutoff
            to_delete.append([iPeak, iTrough])

        # # remove peak/trough pair it's within [0, SCANNING_RANGE] or [wvlTotal - SCANNING_RANGE, wvlTotal] (instead of checking through whole window range)
        # elif fderivPeaks[iPeak] < int(SCANNING_RANGE / wvlS) or fderivTroughs[iTrough] > wvlTotal - int(SCANNING_RANGE / wvlS):
        #     to_delete.append([iPeak, iTrough])
    
    for peak, trough in to_delete: iPairs.remove([peak, trough]) # remove bad pairs
    to_delete.clear()

    # check for close trough before peak or close peak after trough
    for iPeak, iTrough in iPairs:
        
        # check if index in range too
        if (iTrough > 0 and fderivPeaks[iPeak] - fderivTroughs[iTrough - 1] < F_NEIGHBOR_X_DIST) or (iPeak < len(fderivPeaks) - 1 and fderivPeaks[iPeak + 1] - fderivTroughs[iTrough] < F_NEIGHBOR_X_DIST):
            to_delete.append([iPeak, iTrough])
    
    for peak, trough in to_delete: iPairs.remove([peak, trough]) # remove bad pairs

    # transform to np array and check if nonempty
    if not len(iPairs): return # no good pairs
    
    # vals: alignment + amplitude
    valPeaks, _ = find_peaks(vals, prominence = VAL_PEAK_PROM) # potential modes
    if not valPeaks.size: return
    iPeaks = [] # selected peaks

    for iVal in valPeaks:

        for iPeak, iTrough in iPairs:

            fPeak = fderivPeaks[iPeak] # index of peak points
            fTrough = fderivTroughs[iTrough] # index of trough points
            
            # check if peak is within range of this analyzed fderiv section
            if fPeak <= iVal and iVal <= fTrough:

                # index of bounds of flat sections
                endpointLL, endpointLH, endpointRL, endpointRH = fPeak, fPeak, fTrough, fTrough

                # left side
                i1 = fPeak
                while i1 >= 0 and endpointLL == fPeak:

                    if fderiv[i1] <= 0.1: # dips down into flat section
                        if endpointLH == fPeak: endpointLH = i1
                        else: endpointLL == i1 # curves up over threshold before dropping down enough; just take this value as next endpoint

                    if endpointLH != fPeak and fderiv[i1] < -0.1: # out of flat section range now
                        if endpointLL == fPeak: endpointLL = i1 # technically should go with i1 + 1, but i1 isn't too diff and will prevent possible indexing errors (TODO: maybe fix if necessary)

                    i1 -= 1

                # right side
                i2 = fTrough
                while i2 < fderiv.size and endpointRH == fTrough:

                    if fderiv[i2] >= -0.1:
                        if endpointRL == fTrough: endpointRL = i2
                        else: endpointRH == i2

                    if endpointRL != fTrough and fderiv[i2] > 0.1:
                        if endpointRH == fTrough: endpointRH = i2

                    i2 += 1
                
                if endpointLH == fPeak or endpointRL == fTrough: continue # curve never flattens down enough within scanned range of one/both of the "flat sections"
                if endpointLL == fPeak: endpointLL = 0 # flat section continues until end of scanned range
                if endpointRH == fTrough: endpointRH = fderiv.size - 1 # flat section continues until end of scanned range
                
                # flat length
                if endpointLL != 0 and endpointLL - endpointLH < FLAT_LENGTH: continue
                if endpointRH != fderiv.size - 1 and endpointRH - endpointRL < FLAT_LENGTH: continue

                # median-based flat alignment
                medianL = round((endpointLH - endpointLL) / 2)
                medianR = round((endpointRH - endpointRL) / 2)
                if vals[medianR] - vals[medianL] < FLAT_ALIGNMENT:

                    # amplitude filter
                    if np.round(min(vals[medianR], vals[medianL]), 1) >= MIN_AMPLITUDE:
                        iPeaks.append(iVal)
                        break
    
    iPeaks = np.array(iPeaks)
    if not iPeaks.size: return # no values found
    
    # save peak indices to graph later
    if not os.path.exists(filename[:-9] + "iPeaks.txt"):
        iPeaks_files = open(filename[:-9] + "iPeaks.txt", "w+")
        iPeaks_files.write(f"{iPeaks.size}\n")
        for i in range(iPeaks.size): iPeaks_files.write(f"{iPeaks[i]}\n")
    
    return

def Pcalc(wvlTotal, filename, RDATA):

    # t = int(filename[-41:-37])
    # ra = float(filename[-34:-29])
    # a = int(filename[-4:-1])

    # compile data
    reference = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, "reference.txt"), reference)
    incident = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, "incident.txt"), incident)
    z0 = np.zeros(wvlTotal, dtype=np.complex_); complex_compile(os.path.join(RDATA, filename + "field.txt"), z0)
    reflected = z0 - incident

    ratioR = (reference / reflected).real
    ratioI = (reference / reflected).imag
    phase = np.zeros(ratioR.size)
    for i in range(phase.size): phase[i] = math.atan2(ratioI[i], ratioR[i]) + 2 * math.pi if ratioI[i] < 0 else math.atan2(ratioI[i], ratioR[i])
    phase /= math.pi

    if not os.path.exists(filename + "phase.txt"):
        phase_files = open(filename + "phase.txt", "w+")
        phase_files.write(f"{phase.size}\n")
        for i in range(phase.size): phase_files.write(f"{phase[i]}\n")

    # plotname = os.path.join(R3PLOTS, filename + ".jpg")

    # global paramList
    # for i in range(phase.size): 
    #     if round(phase[i], 1) == 2.0: 
    #         paramList.append([t1, t2, ra1, ra2, a])

    #         plotname = os.path.join(PLOTS, filename + ".jpg")
            
    #         break

    # vals = np.zeros(wvlTotal)
    # compile(os.path.join(RDATA, filename + "power.txt"), vals)
    # peaks = np.zeros(line1(filename + "iPeaks.txt"))
    # compile(os.path.join(RDATA, filename + "iPeaks.txt"), peaks)
    
    # # graph
    # ax = plt.subplot()
    # ax.plot(np.linspace(wvlL, wvlH, wvlTotal), vals)
    # ax.plot(wvlL + int(peaks) * wvlS, vals[int(peaks)], 'xr')
    # ax.set_title(f"{CHAR}: t = {round(t1, 4)}, {round(t2, 4)}, ra = {round(ra1, 4)}, {round(ra2, 4)}, a = {round(a, 4)}")
    # ax.set_xlabel("Wavelength (nm)")
    # ax.set_ylabel(f"{CHAR} (%)")
    # ax.set_ylim([0, 1])

    # # phase shift graph (on same sest of axes as the one above)
    # ax = ax.twinx()
    # ax.plot(np.linspace(wvlL, wvlH, wvlTotal), phase, 'y:')
    # ax.set_ylabel("Phase Shift (π) [dotted yellow]")
    # ax.set_ylim([0, 2])

    # # save plot...
    # if not os.path.exists(plotname): plt.savefig(plotname, bbox_inches = "tight") # ...if it doesn't already exist

    # # clear plot
    # plt.clf()

    return

# bundle for 2 functions above
def R3calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, filenames, RDATA):
    cnt = 0
    for fname in filenames:
        filename = os.path.join(RDATA, fname)
        R3calc(wvlS, wvlTotal, filename)
        cnt += 1

    filenames.clear()
    filenames += [os.path.join(RDATA, fname[:-10]) for fname in sorted(os.listdir(RDATA)) if os.path.join(RDATA, fname).endswith("iPeaks.txt")]
    total = len(filenames)
    for filename in filenames:
        fieldSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, filename, cnt, total, RDATA)
        Pcalc(wvlL, wvlH, wvlS, wvlTotal, filename, RDATA)

def R3(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, DATA, RDATA):
    
    root_filename = os.path.join(DATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_")
    
    T_sp_filename = root_filename + "selected_parameters_T.txt"
    R_sp_filename = root_filename + "selected_parameters_R.txt"
    A_sp_filename = root_filename + "selected_parameters_A.txt"
    # if os.path.exists(sp_filename): return

    # unzip filenames to analyze
    T_file = open(root_filename + "final_ra_selected_T.txt", "r")
    T_filenames = T_file.read().splitlines()
    R3calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, T_filenames, RDATA)

    R_file = open(root_filename + "final_ra_selected_R.txt", "r")
    R_filenames = R_file.read().splitlines()
    R3calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, R_filenames, RDATA)

    A_file = open(root_filename + "final_ra_selected_A.txt", "r")
    A_filenames = A_file.read().splitlines()
    R3calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, A_filenames, RDATA)
    
    # create parameter output file
    selected_parameters = open(T_sp_filename, "w+")
    selected_parameters.write("selected [t, r/a ratio, a] values (separated by commas below):\n")
    for t, ra, a in paramList: selected_parameters.write(f"{t},{ra},{a}\n")
    selected_parameters.close()

    selected_parameters = open(R_sp_filename, "w+")
    selected_parameters.write("selected [t, r/a ratio, a] values (separated by commas below):\n")
    for t, ra, a in paramList: selected_parameters.write(f"{t},{ra},{a}\n")
    selected_parameters.close()

    selected_parameters = open(A_sp_filename, "w+")
    selected_parameters.write("selected [t, r/a ratio, a] values (separated by commas below):\n")
    for t, ra, a in paramList: selected_parameters.write(f"{t},{ra},{a}\n")
    selected_parameters.close()

    return