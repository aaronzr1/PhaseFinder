from logging import root
from aux_functions.sim import *
from aux_functions.calc import *
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# R1calc constants
F_PEAK_PROM = 0.08 # fderiv peak prominence
F_TROUGH_PROM = 0.08 # fderiv trough prominence
F_SPIKE_WIDTH = 5 # min x-dist between spike bins
F_FLAT_WIDTH = 20 # min x-dist of flat zones
MIN_AMPLITUDE = 0.8 # min flat zone amplitude

# chosen filenames
fList = []

# save filenames for good spike zone bins (& high flat zone amplitude)
def R1calc(mat, eps, wvlS, wvlTotal, t, ra, a, cnt, total, CHAR, RDATA):

    print(f"R1: Verifying shape of current parameter combination\tt ({cnt}/{total}): {t}\t ra: {ra}\t a: {a}", end = "\n"); sys.stdout.flush()
    
    root_filename = os.path.join(RDATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_t" + str(int(round(t, 4) * 1000)).zfill(4) + "_ra" + "{:.3f}".format(round(ra, 4)) + f"_a{int(round(a, 4) * 1000)}")
    vals_filename = root_filename + "_" + str(CHAR) + ".txt"

    vals = np.zeros(wvlTotal)
    compile(vals_filename, vals)
    fderiv = np.diff(vals); fderiv /= wvlS # first derivative

    # check for peaks and troughs
    fderivPeaks, _ = find_peaks(fderiv, prominence = F_PEAK_PROM)
    fderivTroughs, _ = find_peaks(-fderiv, prominence = F_TROUGH_PROM)
    if not (fderivPeaks.size and fderivTroughs.size): return # no peaks and/or no troughs

    fderivSpikes = np.concatenate((fderivPeaks, fderivTroughs))
    fderivSpikes = np.sort(fderivSpikes)
    fderivBins = []

    # bin spikes
    fderivBins.append(0)
    fderivBins.append(fderivSpikes[0])
    prev = 0
    for i in range(1, fderivSpikes.size):
        if fderivSpikes[i] - fderivSpikes[prev] > F_SPIKE_WIDTH:
            fderivBins.append(fderivSpikes[prev])
            fderivBins.append(fderivSpikes[i])
            prev = i
    fderivBins.append(fderivSpikes[fderivSpikes.size - 1])
    fderivBins.append(fderiv.size - 1)

    fderivBins = np.array(fderivBins)
    fderivBinDist = np.diff(fderivBins)
    fderivBinDist = fderivBinDist[::2] # only take even indices
    flat1_start, flat1_end, flat2_start, flat2_end = -1, -1, -1, -1

    # check flat zones width between spike bins
    for i in range(1, fderivBinDist.size):
        if fderivBinDist[i - 1] > F_FLAT_WIDTH and fderivBinDist[i] > F_FLAT_WIDTH:

            flat1_start = fderivBins[(i - 1) * 2]
            flat1_end = fderivBins[(i - 1) * 2 + 1]
            flat2_start = fderivBins[i * 2]
            flat2_end = fderivBins[i * 2 + 1]
            med1 = int(flat1_start + (flat1_end - flat1_start) / 2)
            med2 = int(flat2_start + (flat2_end - flat2_start) / 2)

            # high amplitude
            if vals[med1] > MIN_AMPLITUDE and vals[med2] > MIN_AMPLITUDE:
                
                # store filename
                global fList
                fList.append(vals_filename)

                return
    
    return

def R1(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raL, raH, raTotal, aL, aH, aTotal, DATA, RDATA):
    
    root_filename = os.path.join(DATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_")
    T_i_ra_filename = root_filename + "initial_ra_selected_T.txt"
    R_i_ra_filename = root_filename + "initial_ra_selected_R.txt"
    A_i_ra_filename = root_filename + "initial_ra_selected_A.txt"
    # if os.path.exists(A_i_ra_filename): return

    T_t_filename = root_filename + "t_selected_T.txt"
    TtList = np.zeros(line1(T_t_filename))
    compile(T_t_filename, TtList)

    R_t_filename = root_filename + "t_selected_R.txt"
    RtList = np.zeros(line1(R_t_filename))
    compile(R_t_filename, RtList)
    
    A_t_filename = root_filename + "t_selected_A.txt"
    AtList = np.zeros(line1(A_t_filename))
    compile(A_t_filename, AtList)

    cnt = 1
    for ra in np.linspace(raL, raH, raTotal):
        ra = round(ra, 4)

        for a in np.linspace(aL, aH, aTotal):
            a = round(a, 4)
            r = round(ra * a, 4)

            for t in TtList:
                t = round(t, 4)
                powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, r, a, cnt, int(len(TtList)), "R1", RDATA)
                R1calc(mat, eps, wvlS, wvlTotal, t, ra, a, cnt, int(len(TtList)), "T", RDATA)

            for t in RtList:
                t = round(t, 4)
                powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, r, a, cnt, int(len(RtList)), "R1", RDATA)
                R1calc(mat, eps, wvlS, wvlTotal, t, ra, a, cnt, int(len(RtList)), "R", RDATA)

            for t in AtList:
                t = round(t, 4)
                powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, r, a, cnt, int(len(AtList)), "R1", RDATA)
                R1calc(mat, eps, wvlS, wvlTotal, t, ra, a, cnt, int(len(AtList)), "A", RDATA)
            
        cnt += 1

    # save selected filenames
    T_i_ra_out = open(T_i_ra_filename, "w+")
    for filename in fList: T_i_ra_out.write(f"{filename}\n")
    T_i_ra_out.close()

    R_i_ra_out = open(R_i_ra_filename, "w+")
    for filename in fList: R_i_ra_out.write(f"{filename}\n")
    R_i_ra_out.close()

    A_i_ra_out = open(A_i_ra_filename, "w+")
    for filename in fList: A_i_ra_out.write(f"{filename}\n")
    A_i_ra_out.close()

    # # save graphs
    # selected_files = open(initial_ra_filename, "r")
    # fnames = selected_files.read().splitlines()
    # for fname in fnames:

    #     t = int(fname[-50:-46])
    #     ra = float(fname[-43:-38])
    #     a = int(fname[-13:-10])
        
    #     filename = os.path.join(RDATA, fname)

    #     vals = np.zeros(wvlTotal)
    #     compile(filename, vals)

    #     # graph
    #     ax = plt.subplot()
    #     ax.plot(np.linspace(wvlL, wvlH, wvlTotal), vals)
    #     # ax.plot(np.linspace(wvlL, wvlH - wvlS, wvlTotal - 1), fderiv)
    #     ax.set_title(f"{CHAR}: t = {round(t1, 4)}, {round(t2, 4)}, ra = {round(ra1, 4)}, {round(ra2, 4)}, a = {round(a, 4)}")
    #     ax.set_xlabel("Wavelength (nm)")
    #     ax.set_ylabel(f"{CHAR} (%)")
    #     ax.set_ylim([0, 1])

    #     # save plot...
    #     plotname = os.path.join(RPLOTS, fname[:-4] + ".jpg")
    #     if not os.path.exists(plotname): plt.savefig(plotname, bbox_inches = "tight") # ...if it doesn't already exist

    #     # clear plot
    #     plt.clf()
    
    return