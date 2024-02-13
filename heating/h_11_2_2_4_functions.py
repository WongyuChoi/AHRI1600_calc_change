# AHRI 1600
# functions: Variable capacity system
from typing import Union, Optional
from typing import List, Tuple

def get_delta_Low_cutout_factor(
        tj_i: float,
        q_dot_Low_tj: float,
        P_Low_tj: float,
        is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: bool = False,
        t_on: Union[float, None] = None,
        t_off: Union[float, None] = None,
        # use_COP: bool = False,
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

    if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature:
        if tj_i <= t_off:
            return 0
        elif t_off < tj_i <= t_on:
            return 0.5
        elif tj_i > t_on:
            return 1
    else:
        return 1

def get_delta_Int_cutout_factor(
        tj_i: float,
        q_dot_Int_tj: float,
        P_Int_tj: float,
        is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: bool = False,
        t_on: Union[float, None] = None,
        t_off: Union[float, None] = None,
        # use_COP: bool = False,
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

    if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature:
        if tj_i <= t_off:
            return 0
        elif t_off < tj_i <= t_on:
            return 0.5
        elif tj_i > t_on:
            return 1
    else:
        return 1


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
    use_COP: bool = False,
) -> float:
    """
    get Heat pump low-temperature cutout factor eqn. 11.129, 11.130, 11.131
    :param tj_i: bin temperature
    :param q_dot_Full_tj: calculated heating capacity at each bin temperature
    :param P_Full_tj: calculated heating power at each bin temperature
    :param t_on: the outdoor temperature at which the compressor reinitiates operation
    :param t_off: The outdoor temperature below which the compressor ceases to operate
    :param is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: the controls of the unit prohibit compressor operation based on outdoor temperature?
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

    COP = q_dot_Full_tj / (3.412 * P_Full_tj)
    if use_COP: # goes regardless cut off temperature in Appendix M1
        if COP < 1:
            return 0
    if is_unit_prohibit_compressor_operation_based_on_outdoor_temperature:
        if use_COP:
            if COP < 1:
                return 0
            else:
                if tj_i <= t_off:
                    return 0
                elif t_off < tj_i <= t_on:
                    return 0.5
                elif tj_i > t_on:
                    return 1
        else:
            if tj_i <= t_off:
                return 0
            elif t_off < tj_i <= t_on:
                return 0.5
            elif tj_i > t_on:
                return 1
    else:
        return 1

def get_E_tj(
    N_j_h: float,
    BL_tj: float,
    q_dot_Full_tj: float,
    q_dot_Int_tj: float,
    q_dot_Low_tj: float,
    P_Full_tj: float,
    P_Int_tj: float,
    P_Low_tj: float,
    delta_Full_cutout: float,
    delta_Int_cutout: float,
    delta_Low_cutout: float,
    is_variable_capacity_certified_two_stage_system: bool,
    C_h_D: float,
) -> Tuple[float, float, float, float]:
    """
    a function to get compressor energy consumption E_tj (Wh)
    :param N_j_h: a capacity range ratio (N_j_h) between Full and Low
    :param BL_tj: calculated heating capacity at each bin temperature
    :param q_dot_Low_tj: calculated Low stage heating capacity at each bin temperature
    :param q_dot_Full_tj: calculated Full stage heating capacity at each bin temperature
    :param q_dot_Int_tj: calculated Int stage heating capacity at each bin temperature
    :param P_Full_tj: calculated Full stage power consumption at each bin temperature
    :param P_Int_tj: calculated Int stage power consumption at each bin temperature
    :param P_Low_tj: calculated Low stage power consumption at each bin temperature
    :param delta_Full_cutout: a cutout factor (delta) for Full stage
    :param delta_Int_cutout: a cutout factor (delta) for Int stage
    :param delta_Low_cutout: a cutout factor (delta) for Low stage
    :param is_variable_capacity_certified_two_stage_system: True if the system is variable capacity certified two stage system
    :param C_h_D: a coefficient for defrost energy consumption
    :return: E_tj
    """
    X_k1 = 0.0
    X_k2 = 0.0
    PLF = 0.0
    ratioPower = 0.0
    COP_inter = 0.0
    if is_variable_capacity_certified_two_stage_system:
        if q_dot_Low_tj >= BL_tj:
            X_k1 = BL_tj / q_dot_Low_tj 
            X_j = X_k1 * delta_Low_cutout
            PLF = 1 - C_h_D * (1 - X_k1)
            ratioPower = X_k1 * P_Low_tj * delta_Low_cutout/ PLF * N_j_h
        elif q_dot_Low_tj < BL_tj < q_dot_Full_tj:
            X_k1 = (q_dot_Full_tj - BL_tj) / (q_dot_Full_tj - q_dot_Low_tj) 
            X_j = delta_Low_cutout
            X_k2 = 1 - X_k1
            ratioPower = (X_k1 * P_Low_tj + X_k2 *
                          P_Full_tj) * delta_Low_cutout * N_j_h

        if BL_tj >= q_dot_Full_tj:
            X_j = delta_Full_cutout
            ratioPower = P_Full_tj * delta_Full_cutout * N_j_h

    else:
        COP_k1 = q_dot_Low_tj / P_Low_tj / 3.412
        COP_k2 = q_dot_Full_tj / P_Full_tj / 3.412
        COP_kv = q_dot_Int_tj / P_Int_tj / 3.412

        if q_dot_Low_tj >= BL_tj:
            X_k1 = BL_tj / q_dot_Low_tj 
            X_j = X_k1 * delta_Low_cutout
            PLF = 1 - C_h_D * (1 - X_k1)
            ratioPower = X_k1 * P_Low_tj * delta_Low_cutout / PLF * N_j_h
        elif q_dot_Low_tj < BL_tj < q_dot_Full_tj:
            if q_dot_Low_tj < BL_tj <= q_dot_Int_tj:
                COP_inter = COP_k1 + (COP_kv - COP_k1) / (q_dot_Int_tj -
                                                          q_dot_Low_tj) * (BL_tj - q_dot_Low_tj)
            else:  # q_dot_Int_tj < BL_tj and BL_tj <= q_dot_Full_tj:
                COP_inter = COP_kv + (COP_k2 - COP_kv) / (q_dot_Full_tj -
                                                          q_dot_Int_tj) * (BL_tj - q_dot_Int_tj)
            X_j = delta_Int_cutout  # variable mode is always on
            ratioPower = BL_tj / (3.412 * COP_inter) * delta_Int_cutout * N_j_h

        elif BL_tj >= q_dot_Full_tj:
            X_j = delta_Full_cutout
            X_k2 = 1
            ratioPower = P_Full_tj * delta_Full_cutout * N_j_h

    return PLF, X_k1, X_k2, X_j, ratioPower
