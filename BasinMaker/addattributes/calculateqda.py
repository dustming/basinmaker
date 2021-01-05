from scipy.optimize import curve_fit


def func_Q_DA(A, k, c):
    return k * A ** c


def return_k_and_c_in_q_da_relationship(da_q):

    try:
        popt2, pcov2 = curve_fit(func_Q_DA, da_q[:, 0], da_q[:, 1])
    except RuntimeError:
        print("#######################################################")
        popt2 = np.full(2, -1)

    print(popt2)
    print(pcov2)
    return popt2[0], popt2[1]


def calculateChannaln(width, depth, Q, slope):
    zch = 2
    sidwd = zch * depth  ###river side width
    tab = "          "
    botwd = width - 2 * sidwd  ### river
    if botwd < 0:
        botwd = 0.5 * width
        sidwd = 0.5 * 0.5 * width
        zch = (width - botwd) / 2 / depth
    Ach = botwd * depth + 2 * zch * depth * depth / 2
    #    arcpy.AddMessage(depth)
    #    arcpy.AddMessage(zch)
    #    arcpy.AddMessage(botwd)
    #    arcpy.AddMessage(width)
    #    arcpy.AddMessage(slope)

    Pch = botwd + 2 * depth * (1 + zch ** 2) ** 0.5
    Rch = float(Ach) / float(Pch)  ### in meter
    V = float(Q) / float(Ach)
    if V > 0:
        n = (Rch ** (2.0 / 3.0)) * (slope ** (1.0 / 2.0)) / V
    else:
        n = -1.2345
    return n
