# AHRI 1600: 11.2.1.3
# functions: Variable capacity System
from typing import List, Union, Tuple

def get_q_tj(
    N_j_c: float,
    BL_tj: float,
    q_dot_Low_tj: float,
    q_dot_Full_tj: float,
) -> Tuple[float, float, float, float]:
    """
    get cooling energy (Btu) at each bin temperature
    In R version, ratioTotalCooling
    :param N_j_c: cooling conditioning hours at each bin temperature
    :param BL_tj: building load at each bin temperature
    :param q_dot_Low_tj: Low stage capacity at each bin temperature
    :param q_dot_Full_tj: Full stage capacity at each bin temperature
    :return: cooling energy at each bin temperature
    """
    ratio_total_cool = 0.0
    X_Tj_k1 = 0.0
    X_Tj_k2 = 0.0
    X_j = 0.0

    if q_dot_Low_tj >= BL_tj:
        X_Tj_k1 = BL_tj / q_dot_Low_tj
        X_j = X_Tj_k1
        ratio_total_cool = X_Tj_k1 * q_dot_Low_tj * N_j_c

    elif q_dot_Low_tj < BL_tj < q_dot_Full_tj:
        X_j = 1 # variable load always on
        ratio_total_cool = BL_tj * N_j_c

    if BL_tj >= q_dot_Full_tj:
        X_j = 1 # full load always on
        ratio_total_cool = q_dot_Full_tj * N_j_c

    return X_Tj_k1, X_Tj_k2, X_j, ratio_total_cool