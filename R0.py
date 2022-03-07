from aux_functions.sim import *
from aux_functions.calc import *
from scipy.signal import find_peaks

# selects t values
def R0calc(mat, eps, tL, tTotal, RDATA, DATA):

    root_filename = os.path.join(DATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_t_")
    T_f_filename = root_filename + "selected_T.txt"
    R_f_filename = root_filename + "selected_R.txt"
    A_f_filename = root_filename + "selected_A.txt"
    if os.path.exists(A_f_filename): return # don't run simulation if data already exists
    root_filename = os.path.join(RDATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_t_")


    print("Finding peak values for t...\n"); sys.stdout.flush()

    T_i_filename = root_filename + "T.txt"
    TtVals = np.zeros(tTotal)
    compile(T_i_filename, TtVals)
    
    R_i_filename = root_filename + "R.txt"
    RtVals = np.zeros(tTotal)
    compile(R_i_filename, RtVals)
    
    A_i_filename = root_filename + "A.txt"
    AtVals = np.zeros(tTotal)
    compile(A_i_filename, AtVals)


    TtPeaks, _ = find_peaks(TtVals, prominence = 0.1)
    TtList = []
    TtList.append([tL + TtPeaks[i] / 1000 for i in range(TtPeaks.size)])
    TtList = np.ndarray.tolist(np.array(np.squeeze(np.array(TtList)), ndmin = 1)) # remove extra (useless) dimension
    
    RtPeaks, _ = find_peaks(RtVals, prominence = 0.1)
    RtList = []
    RtList.append([tL + RtPeaks[i] / 1000 for i in range(RtPeaks.size)])
    RtList = np.ndarray.tolist(np.array(np.squeeze(np.array(RtList)), ndmin = 1)) # remove extra (useless) dimension    TtPeaks, _ = find_peaks(TtVals, prominence = 0.1)
    
    AtPeaks, _ = find_peaks(AtVals, prominence = 0.1)
    AtList = []
    AtList.append([tL + AtPeaks[i] / 1000 for i in range(AtPeaks.size)])
    AtList = np.ndarray.tolist(np.array(np.squeeze(np.array(AtList)), ndmin = 1)) # remove extra (useless) dimension


    out_Tt = open(T_f_filename, "w+")
    out_Tt.write(f"{len(TtList)}\n")
    for t in TtList: out_Tt.write(f"{t}\n")
    out_Tt.close()
    
    out_Rt = open(R_f_filename, "w+")
    out_Rt.write(f"{len(RtList)}\n")
    for t in RtList: out_Rt.write(f"{t}\n")
    out_Rt.close()

    out_At = open(A_f_filename, "w+")
    out_At.write(f"{len(AtList)}\n")
    for t in AtList: out_At.write(f"{t}\n")
    out_At.close()

def R0(mat, eps, NumBasis, wvlVal, wvlL, wvlH, wvlTotal, tL, tH, tTotal, DATA, RDATA):

    incSim(wvlL, wvlH, wvlTotal, NumBasis, RDATA)
    refSim(wvlL, wvlH, wvlTotal, NumBasis, RDATA)

    thicknessSim(mat, eps, NumBasis, wvlVal, tL, tH, tTotal, RDATA)
    R0calc(mat, eps, tL, tTotal, RDATA, DATA)

    return