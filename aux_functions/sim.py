import S4 as S4
import sys, os
import numpy as np
import math

tBuff = 10.00
tPath = 50.00
initial = 1.00

# create S4 simulation object and initializes everything
def init(mat, eps, a, num):
    
    S = S4.New(Lattice = ((a, 0), (0, a)), NumBasis = num)

    # S.AddMaterial(Name = "SiO2", Epsilon = 2.31953 + 0j)
    S.AddMaterial(Name = "Glass", Epsilon = 2.1025 + 0j) # TODO
    S.AddMaterial(Name = "Graphene", Epsilon = 9 + 0j) # TODO
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)
    S.AddMaterial(Name = str(mat), Epsilon = eps)

    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    # S.AddLayer(Name = 'Graphene', Thickness = 0.00034, Material = 'Graphene')
    S.AddLayer(Name = str(mat), Thickness = initial, Material = str(mat))
    S.AddLayer(Name = 'Substrate', Thickness = 0.00, Material = 'Glass') # change to 100
    # S.AddLayer(Name = 'Reflector', Thickness = 0.00, Material = 'Gold')

    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    return S

# run incident simulation
def incSim(wvlL, wvlH, wvlTotal, num, RDATA):

    filename = os.path.join(RDATA, "incident.txt")
    if os.path.exists(filename): return # don't run simulation if data already exists

    S = S4.New(Lattice = ((1, 0), (0, 1)), NumBasis = num)
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)
    S.AddMaterial(Name = "SiO2", Epsilon = 2.31953 + 0j)

    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    S.AddLayer(Name = 'SOGBelow', Thickness = 0.00, Material = 'Vacuum')
    
    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    output = open(filename, "w+")
    for wvl in np.linspace(wvlL, wvlH, wvlTotal):
        freq = 1000 / wvl
        S.SetFrequency(freq)

        (E, _) = S.GetFields(0, 0, tBuff)
        output.write(f"{round(wvl, 5)},{E[1].real},{E[1].imag}\n")

        print("R0\tincident", round(wvl, 5), end = "\n"); sys.stdout.flush()

    output.close()

    return

# run reference simulation
def refSim(wvlL, wvlH, wvlTotal, num, RDATA):

    filename = os.path.join(RDATA, "reference.txt")
    if os.path.exists(filename): return # don't run simulation if data already exists

    S = S4.New(Lattice = ((1, 0), (0, 1)), NumBasis = num)
    S.AddMaterial(Name = "Vacuum", Epsilon = 1.0 + 0j)
    S.AddMaterial(Name = "SiO2", Epsilon = 2.31953 + 0j)

    S.AddLayer(Name = 'AirAbove', Thickness = 0.00, Material = 'Vacuum')
    S.AddLayer(Name = 'AirBuff', Thickness = tBuff, Material = 'Vacuum')
    S.AddLayer(Name = 'AirPath', Thickness = tPath, Material = 'Vacuum')
    S.AddLayer(Name = 'SOGBelow', Thickness = 0.00, Material = 'Vacuum')
    
    S.SetExcitationPlanewave(IncidenceAngles = (0, 0), sAmplitude = 1 + 0j, pAmplitude = 0 + 0j, Order = 0)

    output = open(filename, "w+")
    for wvl in np.linspace(wvlL, wvlH, wvlTotal):
        freq = 1000 / wvl
        S.SetFrequency(freq)

        (E, _) = S.GetFields(0, 0, tBuff + 2 * tPath)
        output.write(f"{round(wvl, 5)},{E[1].real},{E[1].imag}\n")

        print("R0\treference", round(wvl, 5), end = "\n"); sys.stdout.flush()

    output.close()

    return

# generates T, R, A at target wvl
def thicknessSim(mat, eps, NumBasis, wvlVal, tL, tH, tTotal, RDATA):

    root_filename = os.path.join(RDATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_t_")
    T_filename = root_filename + "T.txt"
    R_filename = root_filename + "R.txt"
    A_filename = root_filename + "A.txt"
    if os.path.exists(A_filename): return # don't run simulation if data already exists 
    
    S = init(mat, eps, initial, NumBasis)
    T_out = open(T_filename, "w+")
    R_out = open(R_filename, "w+")
    A_out = open(A_filename, "w+")
    T_out.write(f"{int(tTotal)}\n")
    R_out.write(f"{int(tTotal)}\n")
    A_out.write(f"{int(tTotal)}\n")

    for t in np.linspace(tL, tH, tTotal):
        t = round(t, 4)

        S.SetLayerThickness(Layer = str(mat), Thickness = t)

        freq = 1000 / wvlVal
        S.SetFrequency(freq)

        _, pyntbacktop = S.GetPowerFlux("AirBuff", 0)
        pyntfowbot, _ = S.GetPowerFlux("Substrate", 0)
        
        R_out.write(str(-pyntbacktop.real) + "\n")
        T_out.write(str(pyntfowbot.real) + "\n")
        A_out.write(str(1 - pyntfowbot.real + pyntbacktop.real) + "\n")

        print(f"R0\tthickness {t}", end = "\n"); sys.stdout.flush()
    
    T_out.close()
    R_out.close()
    A_out.close()
    
    return

# generates power (T, R, A) for wvl range
def powerSim(mat, eps, NumBasis, wvlL, wvlH, wvlTotal, t, ra, r, a, cnt, total, ROUND, RDATA):

    root_filename = os.path.join(RDATA, str(mat) + str(round(math.sqrt(eps.real), 3)) + "_t" + str(int(round(t, 4) * 1000)).zfill(4) + "_ra" + "{:.3f}".format(round(ra, 4)) + f"_a{int(round(a, 4) * 1000)}_")
    T_filename = root_filename + "T.txt"
    R_filename = root_filename + "R.txt"
    A_filename = root_filename + "A.txt"
    if os.path.exists(A_filename): return # don't run simulation if data already exists

    S = init(mat, eps, a, NumBasis)
    T_out = open(T_filename, "w+")
    T_out.write(f"{wvlTotal}\n")
    R_out = open(R_filename, "w+")
    R_out.write(f"{wvlTotal}\n")
    A_out = open(A_filename, "w+")
    A_out.write(f"{wvlTotal}\n")

    for wvl in np.linspace(wvlL, wvlH, wvlTotal):

        S.SetRegionCircle(Layer = str(mat), Material = "Vacuum", Center = (0, 0), Radius = r)
        S.SetLayerThickness(Layer = str(mat), Thickness = t)

        freq = 1000 / wvl
        S.SetFrequency(freq)

        _, pyntbacktop = S.GetPowerFlux("AirBuff", 0)
        pyntfowbot, _ = S.GetPowerFlux("Substrate", 0)
        
        R_out.write(str(-pyntbacktop.real) + "\n")
        T_out.write(str(pyntfowbot.real) + "\n")
        A_out.write(str(1 - pyntfowbot.real + pyntbacktop.real) + "\n")
        
        print(f"{ROUND} (thickness value {cnt}/{total})\tthickness: {round(t, 4)}\tra: {round(ra, 4)}\ta: {round(a, 4)}\twvl: {round(wvl, 5)}", end = "\n"); sys.stdout.flush()
    
    T_out.close()
    R_out.close()
    A_out.close()

    return

# custom field sims
def fieldSim(mat, eps, num, wvlL, wvlH, wvlTotal, filename, cnt, total, RDATA):

    t = int(filename[-18:-14])
    ra = float(filename[-11:-6])
    a = int(filename[-4:-1])
    r = round(ra * a, 4)
    
    field_filename = os.path.join(RDATA, str(int(round(t, 4) * 1000)).zfill(4) + "_ra" + "{:.3f}".format(round(ra, 4)) + f"_a{int(round(a, 4) * 1000)}_field.txt")
    if os.path.exists(field_filename): return # don"t run simulation if data already exists

    S = init(mat, eps, a, num)
    field_out = open(field_filename, "w+")
    field_out.write(f"{wvlTotal}\n")
    for wvl in np.linspace(wvlL, wvlH, wvlTotal):

        S.SetRegionCircle(Layer = str(mat), Material = "Vacuum", Center = (0, 0), Radius = r)
        S.SetLayerThickness(Layer = str(mat), Thickness = t)

        freq = 1000 / wvl
        S.SetFrequency(freq)

        (E, _) = S.GetFields(0, 0, tBuff + 2 * tPath)
        field_out.write(f"{round(wvl, 5)},{E[1].real},{E[1].imag}\n")

        print(f"R3 (file {cnt}/{total})\tthickness: {round(t, 4)}\tra: {round(ra, 4)}\ta: {round(a, 4)}\twvl: {round(wvl, 5)}", end = "\n"); sys.stdout.flush()
    field_out.close()

    return