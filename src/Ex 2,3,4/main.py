# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import math
from math import sqrt
import pandas as pd
from matplotlib import pyplot as plt
import xlsxwriter

N_HEAT_CO2 = 0.057            # (mgCO2/J)
U_BLOW = 0.5                  # (none)
P_BLOW = 0.5 * pow(10, 6)     # (W)
A_FLR = 1.4 * 10000           # (m²) confirmed
U_EXT_CO2 = 0.5               # (none)
O_EXT_CO2 = 7.2 * pow(10, 4)  # (none) confirmed
U_PAD = 0.5                   # (none)
O_PAD = 0                     # (m³/s)
K_TH_SCR = 0.05 * pow(10, -3) #       confirmed
G = 9.81                      # (m/s²) confirmed
P_AIR = 1.215                 # (kg/m³) x
P_TOP = 1.2                   # (kg/m³) x
CAP_CO2_AIR = 4             # (m)
CAP_CO2_TOP = 0.8             # (m)
e = 2.718
CD = 0.75
CW = 0.09
A_ROOF = 0.14 * 10000         # (m²)
A_SIDE = 0                    # (m³)
U_SIDE = 0.5                  # (none)
H_SIDE_ROOF = 0               #
U_VENT_FORCED = 0.5
O_VENT_FORCED = 0             # (m³/s)
N_SIDE = 1.1
N_SIDE_THR = 0
S_INS_SCR = 1
C_LEAKAGE = pow(10, -4)
N_ROOF = 1.3                  # (none)
N_ROOF_THR = 0.9              # (none)
M_CH2O = 30 * pow(10, -3)     # (mgμ/mol)
H_ROOF = 1.6                  #
ALPHA = 0.385                 # (µmolµ/mol)
PAR_CAN = 100.0               # (µmol/m²s)
Q = 0.7                       # (J)
EJ = 37000.0                  # (J/mol)
T_CAN_K = 296.0               # (K)
T_25K = 298.15                # (K)
R_KLT = 8.314                 # (J/molK)
S = 710.0                     # (J/molK)
H = 220000.0                  # (J/mol)
LAI = 2.0                     # (leafs/m²)
J_MAX_25LEAF = 210.0          # (leafsJ)
N_CO2_AIR_STOM = 0.67         # (mol)
C_T = 1.7                     # (none)
T_CAN = 296.0                 # (K)


CO2_OUT = 707.96              # (mg/m³)
T_AIR = 293.25                 # (K)
T_TOP = 294.25                 # (K)
T_OUT = 290.85                 # (K)
V_WIND = 3.2                  # (m/s)
U_ROOF = 0.006                  # (none)
U_TH_SCR = 0                  # (none)
epxilon = pow(10, -8)
# formula 3
def MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr):
    return nHeatCO2 * UBlow * PBlow / AFlr

# formula 4
def MCExtAir(UExtCO2, phiExtCO2, AFlr):
    return UExtCO2 * phiExtCO2 / AFlr

# formula 5
def MCPadAir_1(fPad, CO2Out, CO2Air):
    return fPad * (CO2Out - CO2Air)

def MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air):
    fPad = UPad * phiPad / AFlr
    return fPad * (CO2Out - CO2Air)

# formula 6
def MCAirTop(f_ThScr, CO2Air, CO2Top):
    return f_ThScr * (CO2Air - CO2Top)

# formula 7
def fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop):
    a = UThScr * KThScr * pow(abs(TAir - TTop), 2 / 3)
    PMean_Air = (pAir + pTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(pAir - pTop) / (2 * PMean_Air), 1 / 2)
    return a + b

# formula 9
def MCAirOut(f_VentSide, f_VentForce, CO2Air, CO2Out):
    return (f_VentSide + f_VentForce) * (CO2Air - CO2Out)

# formula 10
def fVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind):
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2) + epxilon)
    TMean_Air = (TAir + TOut) / 2
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * Cw * pow(vWind, 2)
    return a * sqrt(b * c + d)

# formula 11
def nInsScr(sInsScr):
    return sInsScr * (2 - sInsScr)

# formula 12
def fleakage(cleakage, vWind):
    if vWind < 0.25:
        return 0.25 * cleakage
    else:
        return vWind * cleakage

# ppfVentSide with no stack eff
def ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw):
    return Cd * USide * ASide * vWind * sqrt(Cw) / (2 * AFlr)


# formula 13
def fVentSide(n_InsScr, ppf_VentSide, f_leakage, UThScr, ppf_VentRoofSide, nSide, nSide_Thr):
    if nSide >= nSide_Thr:
        return n_InsScr * ppf_VentSide + 0.5 * f_leakage
    else:
        return n_InsScr * (UThScr * ppf_VentSide + (1 - UThScr) * ppf_VentRoofSide * nSide) + 0.5 * f_leakage


# formula 14
def fVentForced(n_InsScr, UVentForced, phiVentForced, AFlr):
    return n_InsScr * UVentForced * phiVentForced / AFlr


# formula 15
def MCTopOut(f_VentRoof, CO2Top, CO2Out):
    return f_VentRoof * (CO2Top - CO2Out)


# formula 16
def fVentRoof(n_InsScr, f_leakage, UThScr, ppf_VentRoofSide, nRoof, nSide, nRoof_Thr, ppf_VentRoof):
    if nRoof >= nRoof_Thr:
        return n_InsScr * ppf_VentRoof + 0.5 * f_leakage
    else:
        return n_InsScr * (UThScr * ppf_VentRoof + (1 - UThScr) * ppf_VentRoofSide * nSide) + 0.5 * f_leakage


# formula 17
def ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind):
    TMeanAir = (TAir + TOut) / 2
    a = Cd * URoof * ARoof / (2 * AFlr)
    b = g * hVent * (TAir - TOut) / 2 / TMeanAir + Cw * pow(vWind, 2)
    return a * sqrt(b)


# formula 18, 19
def P(CO2Air, LAI):
    T_Can_K = 20 + 273
    J_POT = LAI * 210 * math.exp(37000 * (T_Can_K - 298.15) / (8.314 * T_Can_K * 298.15)) * (
                1 + math.exp((710 * 298.15 - 220000) / (8.314 * 298.15))) / (
                        1 + math.exp((710 * T_Can_K - 220000) / (8.314 * T_Can_K)))
    J = (J_POT + 38.5 - math.sqrt(math.pow(J_POT + 38.5, 2) - 2.8 * J_POT * 38.5)) / 1.4
    CO2Stom = 0.67 * CO2Air
    return (J * (CO2Stom - 498.1)) / (4 * (CO2Stom + 2 * 498.1))


def R(CO2Air, P):
    CO2Stom = 0.67 * CO2Air
    return P * 498.1 / CO2Stom


def MCAirCan(P, R, CBuf, CMaxBuf):
    MCH2O = 0.03
    hCBuf = 1
    if CBuf > CMaxBuf:
        hCBuf = 0
    return MCH2O * hCBuf * (P - R)


# formula 1
def dxCO2Air(CO2Air, CO2Top):
    # Calculate MCBlowAir
    nHeatCO2 = N_HEAT_CO2
    UBlow = U_BLOW
    PBlow = P_BLOW
    AFlr = A_FLR
    MC_BlowAir = MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr)

    # Calculate MCExtAir
    UExtCO2 = U_EXT_CO2
    phiExtCO2 = O_EXT_CO2
    MC_ExtAir = MCExtAir(UExtCO2, phiExtCO2, AFlr)

    # Calculate MCPadAir
    UPad = U_PAD
    phiPad = O_PAD
    CO2Out = CO2_OUT
    MCPadAir = MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air)

    # Calculate MCAirCan
    LaI = LAI
    p = P(CO2Air, LaI)
    R = 0
    CBuf = 0
    CMax_Buf = 4000
    MC_AirCan = MCAirCan(p, R, CBuf, CMax_Buf)

    # Calculate MCAirTop
    UThScr = U_TH_SCR
    KThScr = K_TH_SCR
    TAir = T_AIR
    TTop = T_TOP
    g = G
    pAir = P_AIR
    pTop = P_TOP
    f_ThScr = fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MC_AirTop = MCAirTop(f_ThScr, CO2Air, CO2Top)

    # Calculte MCAirOut
    cleakage = C_LEAKAGE
    vWind = V_WIND
    f_leakage = fleakage(cleakage, vWind)

    Cd = CD
    URoof = U_ROOF
    USide = U_SIDE
    ARoof = A_ROOF
    ASide = A_SIDE
    hSideRoof = H_SIDE_ROOF
    TOut = T_OUT
    Cw = CW
    f_VentRoofSide = fVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    ppf_VentSide = ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw)

    nSide = N_SIDE
    nSide_Thr = N_SIDE_THR
    sInsScr = S_INS_SCR
    n_InsScr = nInsScr(sInsScr)
    f_VentSide = fVentSide(n_InsScr, ppf_VentSide, f_leakage, UThScr, f_VentRoofSide, nSide, nSide_Thr)

    UVentForced = U_VENT_FORCED
    phiVentForced = O_VENT_FORCED
    f_VentForced = float(fVentForced(n_InsScr, UVentForced, phiVentForced, AFlr))

    MC_AirOut = MCAirOut(f_VentSide, f_VentForced, CO2Air, CO2Out)

    capCO2Air = CAP_CO2_AIR

    return (MC_BlowAir + MC_ExtAir + MCPadAir - MC_AirCan - MC_AirTop - MC_AirOut) / capCO2Air


# formula 2
def dxCO2Top(CO2Air, CO2Top):

    # Calculate MCAirTop
    UThScr = U_TH_SCR
    KThScr = K_TH_SCR
    TAir = T_AIR
    TTop = T_TOP
    g = G
    pAir = P_AIR
    pTop = P_TOP
    f_ThScr = fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MC_AirTop = MCAirTop(f_ThScr, CO2Air, CO2Top)

    # Calculate MCTopOut
    AFlr = A_FLR
    Cd = CD
    URoof = U_ROOF
    USide = U_SIDE
    ARoof = A_ROOF
    ASide = A_SIDE
    hSideRoof = H_SIDE_ROOF
    TOut = T_OUT
    Cw = CW
    vWind = V_WIND
    f_VentRoofSide = fVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)


    hVent = H_ROOF
    ppf_VentRoof = ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind)


    cleakage = C_LEAKAGE
    f_leakage = fleakage(cleakage, vWind)

    sInsScr = S_INS_SCR
    n_InsScr = nInsScr(sInsScr)
    nRoof = N_ROOF
    nSide = N_SIDE
    nRoof_Thr = N_ROOF_THR
    f_VentRoof = fVentRoof(n_InsScr, f_leakage, UThScr, f_VentRoofSide, nRoof, nSide, nRoof_Thr, ppf_VentRoof)

    CO2Out = CO2_OUT
    MC_TopOut = MCTopOut(f_VentRoof, CO2Top, CO2Out)

    capCO2Top = CAP_CO2_TOP
    return (MC_AirTop - MC_TopOut) / capCO2Top



def dx(CO2_AIR, CO2_TOP):
    return dxCO2Air(CO2_AIR, CO2_TOP), dxCO2Top(CO2_AIR, CO2_TOP)




N = 600000




def euler(dx, init_co2_air, init_co2_top, t0, h):
    number_of_steps = N
    t = np.zeros(number_of_steps + 1)
    f = np.zeros((number_of_steps + 1, 2))
    f[0, :] = init_co2_air, init_co2_top
    t[0] = t0
    data = pd.read_excel("data.xlsx")
    df = pd.DataFrame(data)
    i = 0
    for k in range(number_of_steps):
        global T_AIR, T_OUT, T_TOP, V_WIND, U_ROOF, U_TH_SCR, U_EXT_CO2
        if (k % 300) == 0:
            T_AIR = float(df.at[i, "Tair"])
            T_OUT = float(df.at[i, "Tout"])
            T_TOP = float(df.at[i, "Ttop"])
            V_WIND = float(df.at[i, "Vwind"])
            U_ROOF = U_TH_SCR = 0.5
            #U_ROOF = float(df.at[i, "Uroof"])
            #U_TH_SCR = float(df.at[i, "UThScr"])
            i += 1
        t[k + 1] = t[k] + h
        f[k + 1, 0] = f[k, 0] + h * dx(f[k, 0], f[k, 1])[0]
        f[k + 1, 1] = f[k, 1] + h * dx(f[k, 0], f[k, 1])[1]
    return f, t


def rk4(dx, init_co2_air, init_co2_top, t0, h):
    number_of_steps = N
    t = np.zeros(number_of_steps + 1)
    f = np.zeros((number_of_steps + 1, 2))
    t[0] = t0
    f[0, :] = init_co2_air, init_co2_top
    data = pd.read_excel("data.xlsx")
    df = pd.DataFrame(data)
    i = 0
    for k in range(number_of_steps):
        t[k + 1] = t[k] + h
        if (k % 300) == 0:
            global T_AIR, T_OUT, T_TOP, V_WIND, U_ROOF, U_TH_SCR
            T_AIR = float(df.at[i, "Tair"])
            T_OUT = float(df.at[i, "Tout"])
            T_TOP = float(df.at[i, "Ttop"])
            V_WIND = float(df.at[i, "Vwind"])
            U_ROOF = U_TH_SCR = 0.5
            #U_ROOF = float(df.at[i, "Uroof"])
            #U_TH_SCR = float(df.at[i, "UThScr"])
            i += 1
        k1_0 = dx(f[k, 0], f[k, 1])[0]
        k1_1 = dx(f[k, 0], f[k, 1])[1]
        k2_0 = dx(f[k, 0] + 0.5 * h * k1_0, f[k, 1] + 0.5 * h * k1_0)[0]
        k2_1 = dx(f[k, 0] + 0.5 * h * k1_1, f[k, 1] + 0.5 * h * k1_1)[1]
        k3_0 = dx(f[k, 0] + 0.5 * h * k2_0, f[k, 1] + 0.5 * h * k2_0)[0]
        k3_1 = dx(f[k, 0] + 0.5 * h * k2_1, f[k, 1] + 0.5 * h * k2_1)[1]
        k4_0 = dx(f[k, 0] + h * k3_0, f[k, 1] + h * k3_0)[0]
        k4_1 = dx(f[k, 0] + h * k3_1, f[k, 1] + h * k3_1)[1]
        f[k + 1, 0] = f[k, 0] + (1.0/6.0) * h * (k1_0 + 2 * (k2_0 + k3_0) + k4_0)
        f[k + 1, 1] = f[k, 1] + (1.0/6.0) * h * (k1_1 + 2 * (k2_1 + k3_1) + k4_1)
    return f, t




if __name__ == "__main__":
    init_co2_air = init_co2_top = 790.7675669
    t0 = 0
    h = 1
    f, t = euler(dx, init_co2_air, init_co2_top, t0, h)
    u, time = rk4(dx, init_co2_air, init_co2_top, t0, h)
    plt.plot(t, f[:, 0], label='euler')
    data = pd.read_excel("data.xlsx")
    df = pd.DataFrame(data)
    real = np.zeros(N + 1)
    i = -1
    for k in range(N + 1):
        if k % 300 == 0:
            i += 1
        real[k] = float(df.at[i, "CO2air(real data)"])
    plt.plot(t, real, label='actual')
    plt.xlabel('Time(s)')
    plt.ylabel('CO2 concentration')
    plt.show()
    workbook = xlsxwriter.Workbook('result_euler.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write('A1', 'Time')
    worksheet.write('B1', 'CO2_Air')
    worksheet.write('C1', 'CO2_Top')
    worksheet.write('D1', 'CO2_Air(actual)')
    i = 0
    for count in range(0, 7200, 300):
        worksheet.write(i + 1, 0, count)
        worksheet.write(i + 1, 1, f[count][0])
        worksheet.write(i + 1, 2, f[count][1])
        worksheet.write(i + 1, 3, real[count])
        i += 1
    workbook.close()
    workbook = xlsxwriter.Workbook('result_rk4.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write('A1', 'Time')
    worksheet.write('B1', 'CO2_Air')
    worksheet.write('C1', 'CO2_Top')
    worksheet.write('D1', 'CO2_Air(actual)')
    i = 0
    for count in range(0, 7200, 300):
        worksheet.write(i + 1, 0, count)
        worksheet.write(i + 1, 1, u[count][0])
        worksheet.write(i + 1, 2, u[count][1])
        worksheet.write(i + 1, 3, real[count])
        i += 1
    workbook.close()


