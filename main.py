from R0 import *
from R1 import *
from R2 import *
from R3 import *

def main():

    # make directories
    PATH = os.getcwd()

    LOGS = os.path.join(PATH, "logs")
    if not os.path.isdir(LOGS): os.mkdir(LOGS)

    DATA = os.path.join(PATH, "data")
    if not os.path.isdir(DATA): os.mkdir(DATA)
    RDATA = os.path.join(DATA, "raw_data")
    if not os.path.isdir(RDATA): os.mkdir(RDATA)

    # if not os.path.isdir(RDATA): os.mkdir(RDATA)
    # PLOTS = os.path.join(PATH, "plots")
    # if not os.path.isdir(PLOTS): os.mkdir(PLOTS)
    # R1PLOTS = os.path.join(PLOTS, "R1_plots")
    # if not os.path.isdir(R1PLOTS): os.mkdir(R1PLOTS)
    # R2PLOTS = os.path.join(PLOTS, "R2_plots")
    # if not os.path.isdir(R2PLOTS): os.mkdir(R2PLOTS)
    # R3PLOTS = os.path.join(PLOTS, "R3_plots")
    # if not os.path.isdir(R3PLOTS): os.mkdir(R3PLOTS)

    WINDOW_RANGE = 10
    
    mat = "Silicon"
    eps = complex(12.011, 0)
    NumBasis = 32
    wvlVal = 940
    wvlL, wvlH = wvlVal - WINDOW_RANGE, wvlVal + WINDOW_RANGE
    wvlS = 0.100
    wvlTotal = int((wvlH - wvlL) / wvlS) + 1
    tL, tH, tS = 0.100, 1.000, 0.001
    tTotal = int((tH - tL) / tS) + 1
    raL, raH, raS = 0.380, 0.400, 0.010
    raTotal =int((raH - raL) / raS) + 1
    aL, aH, aS = 0.630, round((wvlVal - 300) / 1000, 3), 0.010 # 0.5 --> 0.64
    aTotal = int((aH - aL) / aS) + 1
    
    """
    R0: 
        - incident, reference
        - thickness sims
        - find/save good t values
    
    R1:
        - sim with coarse step size
        - find/save parameters with thin "spike zones" surrounded by two high-amplitude "flat zones"
    
    R2:
        - check for symmetric (to save), up-down (+ r/a increment), or down-up (- r/a increment) "spike" patterns
        - run corresponding fine step-size sim
        - save newly found symmetric "spikes"

    R3:
        - verify symmetric "spike" with more precise filters
        - determine and evaluate phase shift
        - save good final parameters
    """

    R0(mat, eps, NumBasis, wvlVal, wvlL, wvlH, wvlTotal, tL, tH, tTotal, DATA, RDATA)
    R1(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raL, raH, raTotal, aL, aH, aTotal, DATA, RDATA)
    R2(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, raH, DATA, RDATA)
    R3(mat, eps, NumBasis, wvlL, wvlH, wvlS, wvlTotal, DATA, RDATA)

    return
    

if __name__ == "__main__":
    main()