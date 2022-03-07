from aux_functions.sim import *
from aux_functions.calc import *
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# R2calc1
F_PEAK_PROM = 0.08 # fderiv peak prominence
F_TROUGH_PROM = 0.08 # fderiv trough prominence
F_PAIR_X_DIST = 5 # max x-dist between peak/trough pair
F_PAIR_Y_DIST = 5 # max y-diff between peak/trough pair
F_NEIGHBOR_X_DIST = 8 # min x-dist between pair/another spike

# R2calc2 & R2calc3
VAL_PEAK_PROM = 0.1 # R/T peak prominence
VAL_TROUGH_PROM = 0.1 # R/T trough prominence
VAL_PAIR_X_DIST = 5 # max x-dist between a peak/trough pair
VAL_NEIGHBOR_X_DIST = 8 # max x-dist between pair/another spike

# chosen filenames
T_fList = []
R_fList = []
A_fList = []

# checks for symmetric spike shape (just fderiv test done here, rest will be done in R3)
def R2calc1(wvlS, wvlTotal, filename):

    vals = np.zeros(wvlTotal)            
    compile(filename, vals)
    fderiv = np.diff(vals); fderiv /= wvlS # first derivative

    # find peak points (gives indices of those points)
    fderivPeaks, _ = find_peaks(fderiv, prominence = F_PEAK_PROM)
    fderivTroughs, _ = find_peaks(-fderiv, prominence = F_TROUGH_PROM)
    iPairs = [] # indices of selected peak/trough pairs
    to_delete = [] # pairs to delete (used for filtering iPairs)

    if not (fderivPeaks.size and fderivTroughs.size): return False # no peaks and/or no troughs

    # pair peaks with troughs
    i2 = 0
    for i1 in range(fderivPeaks.size):

        if i2 > fderivTroughs.size: break

        while i2 < fderivTroughs.size and fderivTroughs[i2] < fderivPeaks[i1]: i2 += 1
        if i2 < fderivTroughs.size: iPairs.append([i1, i2])
    
    # check each pair for possibilities
    for iPeak, iTrough in iPairs:

        # make sure peak is close to trough and each is about the same distance away from zero (somewhat symmetrical shape)
        if not (fderivTroughs[iTrough] - fderivPeaks[iPeak] < F_PAIR_X_DIST and abs(fderiv[fderivPeaks[iPeak]]) - abs(fderiv[fderivTroughs[iTrough]]) < F_PAIR_Y_DIST): 
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

    return len(iPairs) # True if nonempty

# checks for up-down shape (+ finer iteration)
def R2calc2(wvlTotal, filename):

    vals = np.zeros(wvlTotal)            
    compile(filename, vals)

    # find peak points (gives indices of those points)
    fderivPeaks, _ = find_peaks(vals, prominence = F_PEAK_PROM)
    fderivTroughs, _ = find_peaks(-vals, prominence = F_TROUGH_PROM)
    iPairs = [] # indices of selected peak/trough pairs
    to_delete = [] # pairs to delete (used for filtering iPairs)

    if not (fderivPeaks.size and fderivTroughs.size): return False # no peaks and/or no troughs

    # systematically pair off each peak with the next trough in range and store
    i2 = 0
    for i1 in range(fderivPeaks.size):

        if i2 > fderivTroughs.size: break

        while i2 < fderivTroughs.size and fderivTroughs[i2] < fderivPeaks[i1]: i2 += 1
        if i2 < fderivTroughs.size: iPairs.append([i1, i2])
    
    # check each pair for possibilities
    for iPeak, iTrough in iPairs:

        # make sure peak is close to trough and each is about the same distance away from zero (somewhat symmetrical shape)
        if not (fderivTroughs[iTrough] - fderivPeaks[iPeak] < VAL_PAIR_X_DIST): 
            to_delete.append([iPeak, iTrough])

        # # remove peak/trough pair it's within [0, SCANNING_RANGE] or [wvlTotal - SCANNING_RANGE, wvlTotal] (instead of checking through whole window range)
        # elif fderivPeaks[iPeak] < int(SCANNING_RANGE / wvlS) or fderivTroughs[iTrough] > wvlTotal - int(SCANNING_RANGE / wvlS):
        #     to_delete.append([iPeak, iTrough])
    
    for peak, trough in to_delete: iPairs.remove([peak, trough]) # remove bad pairs
    to_delete.clear()

    # check for close trough before peak or close peak after trough
    for iPeak, iTrough in iPairs:
        
        # check if index in range too
        if (iTrough > 0 and fderivPeaks[iPeak] - fderivTroughs[iTrough - 1] < VAL_NEIGHBOR_X_DIST) or (iPeak < len(fderivPeaks) - 1 and fderivPeaks[iPeak + 1] - fderivTroughs[iTrough] < VAL_NEIGHBOR_X_DIST): 
            to_delete.append([iPeak, iTrough])
    
    for peak, trough in to_delete: iPairs.remove([peak, trough]) # remove bad pairs

    return len(iPairs) # True if nonempty

# checks for down-up shape (- finer iteration)
def R2calc3(wvlTotal, filename):

    vals = np.zeros(wvlTotal)            
    compile(filename, vals)

    # find peak points (gives indices of those points)
    fderivPeaks, _ = find_peaks(vals, prominence = F_PEAK_PROM)
    fderivTroughs, _ = find_peaks(-vals, prominence = F_TROUGH_PROM)
    iPairs = [] # indices of selected peak/trough pairs
    to_delete = [] # pairs to delete (used for filtering iPairs)

    if not (fderivPeaks.size and fderivTroughs.size): return False # no peaks and/or no troughs

    # systematically pair off each peak with the next trough in range and store
    i2 = 0
    for i1 in range(fderivTroughs.size):

        if i2 > fderivPeaks.size: break

        while i2 < fderivPeaks.size and fderivPeaks[i2] < fderivTroughs[i1]: i2 += 1
        if i2 < fderivPeaks.size: iPairs.append([i1, i2])
    
    # check each pair for possibilities
    for iTrough, iPeak in iPairs:

        # make sure peak is close to trough and each is about the same distance away from zero (somewhat symmetrical shape)
        if not (fderivPeaks[iPeak] - fderivTroughs[iTrough] < VAL_PAIR_X_DIST): 
            to_delete.append([iTrough, iPeak])

        # # remove peak/trough pair it's within [0, SCANNING_RANGE] or [wvlTotal - SCANNING_RANGE, wvlTotal] (instead of checking through whole window range)
        # elif fderivPeaks[iPeak] < int(SCANNING_RANGE / wvlS) or fderivTroughs[iTrough] > wvlTotal - int(SCANNING_RANGE / wvlS):
        #     to_delete.append([iPeak, iTrough])
    
    for trough, peak in to_delete: iPairs.remove([trough, peak]) # remove bad pairs
    to_delete.clear()

    # check for close trough before peak or close peak after trough
    for iTrough, iPeak in iPairs:
        
        # check if index in range too
        if (iPeak > 0 and fderivTroughs[iTrough] - fderivPeaks[iPeak - 1] < VAL_NEIGHBOR_X_DIST) or (iTrough < len(fderivTroughs) - 1 and fderivTroughs[iTrough + 1] - fderivPeaks[iPeak] < VAL_NEIGHBOR_X_DIST): 
            to_delete.append([iTrough, iPeak])
    
    for trough, peak in to_delete: iPairs.remove([trough, peak]) # remove bad pairs

    return len(iPairs) # True if nonempty

# wraps together above three calc functions
def R2calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, raS, filenames, fList, RDATA):

    total = len(filenames)
    cnt = 0

    for fname in filenames:

        root_filename = os.path.join(RDATA, fname)

        # TODO: verify this
        t = int(fname[-23:-19])
        raVal = float(fname[-16:-11])
        a = int(fname[-9:-6])
        CHAR = str(fname[-5:-4])
    
        if R2calc1(wvlS, wvlTotal, root_filename): fList.append(fname)

        if R2calc2(wvlTotal, root_filename):

            for ra in np.linspace(raVal, raH, int((raH - raVal) / raS)):
                
                powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, round(ra * a, 4), a, cnt, total, "R2", RDATA)

                filename = os.path.join(RDATA, str(int(round(t, 4) * 1000)).zfill(4) + "_ra" + "{:.3f}".format(round(ra, 4)) + f"_a{int(round(a, 4) * 1000)}_{CHAR}.txt")

                if R2calc1(wvlS, wvlTotal, filename): fList.append(filename)

                if not R2calc2(wvlTotal, filename): break
        
        if R2calc3(wvlTotal, root_filename):

            for ra in np.linspace(raVal, raH, int((raH - raVal) / raS)):

                powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, round(ra * a, 4), a, cnt, total, "R2", RDATA)

                filename = os.path.join(RDATA, str(int(round(t, 4) * 1000)).zfill(4) + "_ra" + "{:.3f}".format(round(ra, 4)) + f"_a{int(round(a, 4) * 1000)}_{CHAR}.txt")

                if R2calc1(wvlS, wvlTotal, filename): fList.append(filename)

                if not R2calc3(wvlTotal, filename): break
        
        cnt += 1
    
    return

def R2(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, DATA, RDATA):

    raS = 0.001
    
    root_filename = os.path.join(DATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_")
    T_f_ra_filename = root_filename + "final_ra_selected_T.txt"
    R_f_ra_filename = root_filename + "final_ra_selected_R.txt"
    A_f_ra_filename = root_filename + "final_ra_selected_A.txt"
    # if os.path.exists(A_f_ra_filename): return

    # unzip filenames to analyze
    T_i_ra_filename = root_filename + "initial_ra_selected_T.txt"
    T_i_ra_filenames = open(T_i_ra_filename, "r")
    T_filenames = T_i_ra_filenames.read().splitlines()
    R2calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, raS, T_filenames, T_fList, RDATA)

    R_i_ra_filename = root_filename + "initial_ra_selected_R.txt"
    R_i_ra_filenames = open(R_i_ra_filename, "r")
    R_filenames = R_i_ra_filenames.read().splitlines()
    R2calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, raS, R_filenames, R_fList, RDATA)

    A_i_ra_filename = root_filename + "initial_ra_selected_A.txt"
    A_i_ra_filenames = open(A_i_ra_filename, "r")
    A_filenames = A_i_ra_filenames.read().splitlines()
    R2calc(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, raS, A_filenames, A_fList, RDATA)
    
    # save chosen filenames
    selected_files = open(T_f_ra_filename, "w+")
    for filename in T_fList: selected_files.write(f"{filename}\n")
    selected_files.close()

    selected_files = open(R_f_ra_filename, "w+")
    for filename in R_fList: selected_files.write(f"{filename}\n")
    selected_files.close()

    selected_files = open(A_f_ra_filename, "w+")
    for filename in A_fList: selected_files.write(f"{filename}\n")
    selected_files.close()

    # selected_files = open(final_ra_filename, "r")
    # fnames = selected_files.read().splitlines()
    # for fname in fnames:

    #     t1 = int(fname[-50:-46])
    #     raVal1 = float(fname[-43:-38])
    #     t2 = int(fname[-27:-23])
    #     raVal2 = float(fname[-20:-15])
    #     a = int(fname[-13:-10])
        
    #     filename = os.path.join(RDATA, fname)

    #     vals = np.zeros(wvlTotal)
    #     compile(filename, vals)

    #     # graph
    #     ax = plt.subplot()
    #     ax.plot(np.linspace(wvlL, wvlH, wvlTotal), vals)
    #     # ax.plot(np.linspace(wvlL, wvlH - wvlS, wvlTotal - 1), fderiv)
    #     ax.set_title(f"{CHAR}: t = {round(t1, 4)}, {round(t2, 4)}, ra = {round(ra, 4)}, {round(ra, 4)}, a = {round(a, 4)}")
    #     ax.set_xlabel("Wavelength (nm)")
    #     ax.set_ylabel(f"{CHAR} (%)")
    #     ax.set_ylim([0, 1])

    #     # save plot...
    #     plotname = os.path.join(RPLOTS, fname[:-4] + ".jpg")
    #     if not os.path.exists(plotname): plt.savefig(plotname, bbox_inches = "tight") # ...if it doesn't already exist

    #     # clear plot
    #     plt.clf()
    # plt.close()

    return