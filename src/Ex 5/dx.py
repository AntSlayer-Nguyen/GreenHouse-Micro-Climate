from math import sqrt

# --------Cap---------------------

def Cap_VP_Air(M_Water, H_Air, R, T_Air):
    return (M_Water * H_Air) / (R * T_Air)


def Cap_VP_Top(M_Water, H_Top, R, T_Top):
    return (M_Water * H_Top) / (R * T_Top)


# --------MV_BlowAir--------------

def MV_BlowAir(N_Heat_Vap, U_Blow, P_Blow, A_Flr):
    return (N_Heat_Vap * U_Blow * P_Blow) / A_Flr


# ---------MV_FogAir-------------

def MV_FogAir(U_Fog, O_Fog, A_Flr):
    return (U_Fog * O_Fog) / A_Flr


# ----------MV_PadAir------------

def MV_PadAir(U_Pad, O_Pad, P_Air, N_Pad, X_Pad, X_Out, A_Flr):
    return (U_Pad * O_Pad * P_Air * (N_Pad * (X_Pad - X_Out) + X_Out)) / A_Flr


# ------------MV_AirOut_Pad-------

def MV_AirOut_Pad(U_Pad, O_Pad, A_Flr, M_Water, R, VP_Air, T_Air):
    return (U_Pad * O_Pad * M_Water * VP_Air) / A_Flr * R * T_Air


# ------------MV_Air_ThScr--------------

def HEC_Air_ThScr(U_ThScr, T_ThScr, T_Air):
    return 1.7 * U_ThScr * pow(abs(T_Air - T_ThScr), 0.33)


def MV_Air_ThScr(S_MV_Air_ThScr, VP_Air, VP_ThScr, HEC_Air_ThScr):
    return (6.4 * pow(10, -9) * HEC_Air_ThScr * (VP_Air - VP_ThScr)) / (
                1 + pow(2.718, S_MV_Air_ThScr * (VP_Air - VP_ThScr)))


# ----------MV_Top_Cov_in-------------

def HEC_Top_Cov_in(c_Hec_in, T_Top, T_Cov_in, A_Cov, A_Flr):
    return c_Hec_in * pow(abs(T_Top - T_Cov_in), 0.33) * A_Cov / A_Flr


def MV_Top_Cov_in(S_MV_Top_Cov_in, VP_Top, VP_Cov_in, HEC_Top_Cov_in):
    a = 1 + pow(2.718, S_MV_Top_Cov_in * (VP_Top - VP_Cov_in))
    return (6.4 * pow(10, -9) * HEC_Top_Cov_in * (VP_Top - VP_Cov_in)) / a

# ----------MV_AirMech---------------

def HEC_AirMech(U_MechCool, COP_MechCool, P_MechCool, A_Flr, T_Air, T_Mech, Delta_H, VP_Air, VP_MechCool):
    a = U_MechCool * COP_MechCool * P_MechCool / A_Flr
    b = T_Air - T_Mech + 6.4 * pow(10, -9) * Delta_H * (VP_Air - VP_MechCool)
    return a / b


def MV_Air_Mech(S_MV_Air_Mech, VP_Mech, VP_Air, HEC_Air_Mech):
    return (6.4 * pow(10, -9) * HEC_Air_Mech * (VP_Air - VP_Mech)) / (
                1 + pow(2.718, S_MV_Air_Mech * (VP_Air - VP_Mech)))


# ----------f-Formula------------

def f_ThScr(U_ThScr, K_ThScr, T_Air, T_Top, g, p_Air, p_Top):
    a = U_ThScr * K_ThScr * pow(abs(T_Air - T_Top), 2 / 3)
    PMean_Air = (p_Air + p_Top) / 2
    b = (1 - U_ThScr) * pow(g * (1 - U_ThScr) * abs(p_Air - p_Top) / (2 * PMean_Air), 1 / 2)
    return a + b


def fleakage(cleakage, vWind):
    if vWind < 0.25:
        return 0.25 * cleakage
    else:
        return vWind * cleakage


def nInsScr(sInsScr):
    return sInsScr * (2 - sInsScr)


def f_Vent_Roof_Side(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind):
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2) + epxilon)
    TMean_Air = (TAir + TOut) / 2
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * Cw * pow(vWind, 2)
    return a * sqrt(b * c + d)


def ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw):
    return Cd * USide * ASide * vWind * sqrt(Cw) / (2 * AFlr)


def f_VentSide(n_InsScr, ppf_VentSide, f_leakage, UThScr, f_VentRoofSide, nSide, nSide_Thr):
    if nSide >= nSide_Thr:
        return n_InsScr * ppf_VentSide + 0.5 * f_leakage
    else:
        return n_InsScr * (UThScr * ppf_VentSide + (1 - UThScr) * f_VentRoofSide * nSide) + 0.5 * f_leakage


def f_VentForced(n_InsScr, U_VentForced, phi_VentForced, A_Flr):
    return n_InsScr * U_VentForced * phi_VentForced / A_Flr


def ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind):
    TMeanAir = (TAir + TOut) / 2
    part1 = Cd * URoof * ARoof / (2 * AFlr)
    part2 = g * hVent * (TAir - TOut) / 2 / TMeanAir + Cw * pow(vWind, 2)
    return part1 * sqrt(part2)


def fVentRoof(n_InsScr, f_leakage, UThScr, ppf_VentRoofSide, nRoof, nSide, nRoof_Thr, ppf_VentRoof):
    if nRoof >= nRoof_Thr:
        return n_InsScr * ppf_VentRoof + 0.5 * f_leakage
    else:
        return n_InsScr * (UThScr * ppf_VentRoof + (1 - UThScr) * ppf_VentRoofSide * nSide) + 0.5 * f_leakage


# -----------MV_AirTop-----------

def MV_AirTop(M_Water, R, f_Th_Scr, VP_Air, VP_Top, T_Air, T_Top):
    return M_Water * f_Th_Scr * (VP_Air / T_Air - VP_Top / T_Top) / R


# -----------MV_AirOut----------

def MV_AirOut(M_Water, R, VP_Air, VP_Out, T_Air, T_Out, f_vent_Side, f_vent_Forced):
    return M_Water * (f_vent_Forced + f_vent_Side) * (VP_Air / T_Air - VP_Out / T_Out) / R


# -----------MV_TopOut----------

def MV_TopOut(M_Water, R, VP_Top, VP_Out, T_Top, T_Out, f_vent_roof):
    return M_Water * f_vent_roof * (VP_Top / T_Top - VP_Out / T_Out) / R


# -----------MV_CanAir----------

def VEC_CanAir(p_Air, c_p_Air, LAI, Delta_H, y, r_b, r_s):
    return (2 * p_Air * c_p_Air * LAI) / (Delta_H * y * (r_s + r_b))


def MV_CanAir(VEC_CanAir, VP_Can, VP_Air):
    return VEC_CanAir * (VP_Can - VP_Air)


# -----------result_dx----------

def dx(VP_Air, VP_Top):
    # MV_BlowAir
    mv_blow_air = MV_BlowAir(N_HEAT_VAP, U_BLOW, P_BLOW, A_FLR)

    # MV_FogAir
    mv_fog_air = MV_FogAir(U_FOG, O_FOG, A_FLR)

    # MV_PadAir
    mv_pad_air = MV_PadAir(U_PAD, O_PAD, P_AIR, N_PAD, X_PAD, X_OUT, A_FLR)

    # MV_AirOut_Pad
    mv_airout_pad = MV_AirOut_Pad(U_PAD, O_PAD, A_FLR, M_WATER, R, VP_Air, T_AIR)

    # MV_Air_ThScr
    hec_air_thscr = HEC_Air_ThScr(U_THSCR, T_THSCR, T_AIR)
    mv_air_thscr = MV_Air_ThScr(S_MV_12, VP_Air, VP_THSCR, hec_air_thscr)

    # MV_Top_Cov_in
    hec_top_cov_in = HEC_Top_Cov_in(C_HEC_IN, T_TOP, T_COV_IN, A_COV, A_FLR)
    mv_top_cov_in = MV_Top_Cov_in(S_MV_12, VP_Top, VP_COV_IN, hec_top_cov_in)

    # MV_AirMech
    hec_air_mech = HEC_AirMech(U_MECH_COOL, COP_MECHCOOL, P_MECHCOOL, A_FLR, T_AIR, T_MECH, DELTA_H, VP_Air, VP_MECH)
    mv_air_mech = MV_Air_Mech(S_MV_12, VP_MECH, VP_Air, hec_air_mech)

    # MV_CanAir
    vec_can_air = VEC_CanAir(P_AIR, CP_AIR, LAI, DELTA_H, Y, R_B, R_S)
    mv_can_air = MV_CanAir(vec_can_air, VP_CAN, VP_Air)

    # f-Formula
    n_ins_scr = nInsScr(S_INS_SCR)
    f_th_scr = f_ThScr(U_THSCR, K_THSCR, T_AIR, T_TOP, G, P_AIR, P_TOP)
    f_leakage = fleakage(C_LEAKAGE, V_WIND)
    f_vent_roof_side = f_Vent_Roof_Side(CD, A_FLR, U_ROOF, U_SIDE, A_ROOF, A_SIDE, G, H_SIDE_ROOF, T_AIR, T_OUT, CW,
                                        V_WIND)
    ppf_vent_side = ppfVentSide(CD, U_SIDE, A_SIDE, V_WIND, A_FLR, CW)
    f_vent_side = f_VentSide(n_ins_scr, ppf_vent_side, f_leakage, U_THSCR, f_vent_roof_side, N_SIDE, N_SIDE_TH)
    f_vent_forced = f_VentForced(n_ins_scr, U_VENTFORCED, O_VENTFORCED, A_FLR)
    ppf_vent_roof = ppfVentRoof(CD, U_ROOF, A_ROOF, A_FLR, G, H_VENT, T_AIR, T_OUT, CW, V_WIND)
    f_vent_roof = fVentRoof(n_ins_scr, f_leakage, U_THSCR, f_vent_roof_side, N_ROOF, N_SIDE, N_ROOF_TH, ppf_vent_roof)

    # MV_AirTop
    mv_air_top = MV_AirTop(M_WATER, R, f_th_scr, VP_Air, VP_Top, T_AIR, T_TOP)

    # MV_AirOut
    mv_air_out = MV_AirOut(M_WATER, R, VP_Air, VP_OUT, T_AIR, T_OUT, f_vent_side, f_vent_forced)

    # MV_TopOut
    mv_top_out = MV_TopOut(M_WATER, R, VP_Top, VP_OUT, T_TOP, T_OUT, f_vent_roof)

    # dx
    cap_VP_air = Cap_VP_Air(M_WATER, H_AIR, R, T_AIR)
    cap_VP_top = Cap_VP_Top(M_WATER, H_TOP, R, T_TOP)
    dx_VP_air = (mv_can_air + mv_pad_air + mv_fog_air + mv_blow_air - mv_air_mech - mv_air_thscr
                 - mv_air_out - mv_air_top - mv_airout_pad) / cap_VP_air
    dx_VP_top = (mv_air_top - mv_top_cov_in - mv_top_out) / cap_VP_top
    return dx_VP_air, dx_VP_top