# AHRI 1600
# functons: Two-stage System
from typing import List, Union, Tuple, Optional

def get_delta_Low_cutout_factor(
        tj_i: float,
        q_dot_Low_tj: float,
        P_Low_tj: float,
        is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: bool = False,
        t_on: Union[float, None] = None,
        t_off: Union[float, None] = None,
) -> float:
    """
    A function to calculate a cutout factor (delta) for Low stage
    :param tj_i: bin temperature
    :param q_dot_Low_tj: calculated Low stage heating capacity at each bin temperature
    :param P_Low_tj: calculated Low stage power consumption at each bin temperature
    :param is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: the controls of the unit prohibit compressor operation based on outdoor temperature?
    :param t_on: the outdoor temperature at which the compressor reinitiates operation
    :param t_off: The outdoor temperature below which the compressor ceases to operate
    :return: a cutout factor (delta) for Low stage
    """
    delta = 0.0
    # validation: if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature is True, both t_on and t_off should exist
    if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature:
        if not t_on and not t_off:
            raise ValueError(
                "If is_unit_prohibit_compressor_operation_based_on_outdoor_temperature is True, both t_on and t_off should exist.")

    if t_off is None:
        t_off = float('-inf')
    if t_on is None:
        t_on = float('-inf')

    # the latter is different from Appendix M1
    if (tj_i <= t_off) | (q_dot_Low_tj / (3.412 * P_Low_tj) < 1.0):
        delta = 0
    elif t_off < tj_i <= t_on:
        delta = 1/2
    elif tj_i > t_on:
        delta = 1

    return delta


def get_delta_Full_cutout_factor(
    tj_i: float,
    q_dot_Full_tj: float,  # calculated heating capacity at each bin temperature
    P_Full_tj: float,
    # the controls of the unit prohibit compressor operation based on outdoor temperature?
    is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: bool = False,
    # the outdoor temperature at which the compressor reinitiates operation
    t_on: Optional[float] = None,
    # The outdoor temperature below which the compressor ceases to operate
    t_off: Optional[float] = None,
) -> float:
    """
    get Heat pump low-temperature cutout factor eqn. 11.129, 11.130, 11.131
    :param tj_i: bin temperature
    :param q_dot_Full_tj: calculated heating capacity at each bin temperature
    :param P_Full_tj: calculated heating power at each bin temperature
    :param t_on: the outdoor temperature at which the compressor reinitiates operation
    :param t_off: The outdoor temperature below which the compressor ceases to operate
    :return: Heat pump low-temperature cutout factor
    """
    # validation: if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature is True, both t_on and t_off should exist
    if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature:
        if not t_on and not t_off:
            raise ValueError(
                "If is_unit_prohibit_compressor_operation_based_on_outdoor_temperature is True, both t_on and t_off should exist.")

    if t_off is None:
        t_off = float('-inf')
    if t_on is None:
        t_on = float('-inf')

    if tj_i <= t_off or q_dot_Full_tj / (3.412 * P_Full_tj) < 1:
        return 0
    elif t_off < tj_i <= t_on and q_dot_Full_tj / (3.412 * P_Full_tj) >= 1:
        return 0.5
    elif tj_i > t_on and q_dot_Full_tj / (3.412 * P_Full_tj) >= 1:
        return 1

def get_E_tj(
    tj_i: float,
    N_j_h: float,
    BL_tj: float,
    q_dot_Full_tj: float,
    q_dot_Low_tj: float,
    P_Low_tj: float,
    P_Full_tj: float,
    delta_Low_cutout: float,
    delta_Full_cutout: float,
    lock_out_low_capacity_ops: bool,
    od_temp_when_lock_out: Union[float, None],
    C_h_D_Full: float,
    C_h_D_Low: float,
) -> Tuple[float, float, float, float, float]:
    """
    a function to calculate heating energy (Wh) consumption at each bin temperature
    In R version, ratioTotalPower
    :param BL_tj: building load at each bin temperature
    :param N_j_h: fraction of bin hours
    :param q_dot_Full_tj: Full stage capacity at each bin temperature
    :param q_dot_Low_tj: Low stage capacity at each bin temperature
    :param P_Low_tj: Low stage power at each bin temperature
    :param P_Full_tj: Full stage power at each bin temperature
    :param delta_Low_cutout: a cutout factor (delta) for Low stage
    :param delta_Full_cutout: a cutout factor (delta) for Full stage
    :param lock_out_low_capacity_ops: the controls of the unit prohibit compressor operation based on outdoor temperature?
    :param od_temp_when_lock_out: The outdoor temperature below which the compressor ceases to operate
    :param C_h_D_Full: the degradation coefficient for Full stage
    :param C_h_D_Low: the degradation coefficient for Low stage
    :return: heating energy (Wh) consumption at each bin temperature
    """
    X_Tj_k1 = 0.0
    X_Tj_k2 = 0.0
    PLF = 0.0
    X_j = 0.0  # on-time fractions for the lowest compressor stage for bin j.
    E_tj = 0.0

    if BL_tj <= q_dot_Low_tj:
        X_Tj_k1 = BL_tj / q_dot_Low_tj
        X_j = X_Tj_k1
        PLF = 1 - C_h_D_Low * (1 - X_Tj_k1)
        E_tj = X_Tj_k1 * P_Low_tj * delta_Low_cutout / PLF * N_j_h

    elif q_dot_Low_tj < BL_tj < q_dot_Full_tj:  # HLF_Low = X_Tj_k1; HLF_Full = X_Tj_k2
        X_Tj_k1 = (q_dot_Full_tj - BL_tj) / (q_dot_Full_tj - q_dot_Low_tj)
        X_Tj_k2 = 1 - X_Tj_k1
        X_j = X_Tj_k1
        E_tj = (X_Tj_k1 * P_Low_tj + X_Tj_k2 * P_Full_tj) * \
            delta_Low_cutout * N_j_h

    if (q_dot_Low_tj < BL_tj < q_dot_Full_tj) and lock_out_low_capacity_ops:
        if od_temp_when_lock_out:  # when od_temp_when_lock_out is not None
            if tj_i < od_temp_when_lock_out:  # lock out temperature is reached or exceeded, alternating full and off
                X_Tj_k2 = BL_tj / q_dot_Full_tj
                X_Tj_k1 = 1 - X_Tj_k2
                X_j = X_Tj_k2
                PLF = 1 - C_h_D_Full * (1 - X_Tj_k2)
                E_tj = (X_Tj_k2 * P_Full_tj * delta_Full_cutout) / \
                    PLF * N_j_h
            else:  # lock out temperature is not reached, just as the above  q_dot_Low_tj < BL_tj < q_dot_Full_tj
                X_Tj_k1 = (q_dot_Full_tj - BL_tj) / \
                    (q_dot_Full_tj - q_dot_Low_tj)
                X_Tj_k2 = 1 - X_Tj_k1
                X_j = X_Tj_k1
                E_tj = (X_Tj_k1 * P_Low_tj + X_Tj_k2 * P_Full_tj) * \
                    delta_Low_cutout * N_j_h
        else:  # when od_temp_when_lock_out is None, just lock out low capacity operations
            X_Tj_k2 = BL_tj / q_dot_Full_tj
            X_Tj_k1 = 1 - X_Tj_k2
            X_j = X_Tj_k2
            PLF = 1 - C_h_D_Full * (1 - X_Tj_k2)
            E_tj = (X_Tj_k2 * P_Full_tj * delta_Full_cutout) / \
                PLF * N_j_h

    if BL_tj >= q_dot_Full_tj:
        X_Tj_k2 = 1.0
        X_j = 1.0
        E_tj = P_Full_tj * delta_Full_cutout * N_j_h

    return PLF, X_Tj_k1, X_Tj_k2, X_j, E_tj
