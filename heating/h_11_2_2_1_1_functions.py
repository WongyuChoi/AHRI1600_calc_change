# AHRI 1600: 11.2.2.1.1
# functions: Single Stage System with Either a Fixed-Speed Indoor Blower or a Constant-Air-Volume-Rate Indoor Blower, or a Single-Speed Coil-Only System Heat Pump
from typing import List, Union, Tuple, Optional

def get_delta_Full_cutout_factor(
    tj_i: float,
    q_dot_Full_tj: float,  # calculated heating capacity at each bin temperature
    P_Full_tj: float,
    # the controls of the unit prohibit compressor operation based on outdoor temperature?
    is_unit_prohibit_compressor_operation_based_on_outdoor_temperature: bool,
    # the outdoor temperature at which the compressor reinitiates operation
    t_on: Optional[float],
    # The outdoor temperature below which the compressor ceases to operate
    t_off: Optional[float],
    use_COP: bool = False,
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
    P_Full_tj: float,
    delta_Full_cutout: float,
    deg_coeff_heat: float,
) -> Tuple[float, float, float]:
    """
    get total compressor energy (Wh) at each bin temperature
    In R version, ratioTotalPower
    :param N_j_h: heating conditioning hours at each bin temperature
    :param BL_tj: building heating load at bin temperature
    :param q_dot_Full_tj: calculated heating capacity at each bin temperature
    :param P_Full_tj: calculated heating power at each bin temperature
    :param delta_Full_cutout: Heat pump low-temperature cutout factor
    :param deg_coeff_heat: degradation coefficient after validation
    :return: PLF, X_j, E_tj total energy at each bin temperature
    """
    E_tj = 0.0
    HLF_Full =  min(BL_tj / q_dot_Full_tj, 1)
    X_j = min(BL_tj / q_dot_Full_tj, 1) * delta_Full_cutout
    # considering the following revision for Feb.

    PLF_Full_tj = 1 - deg_coeff_heat * (1 - X_j)
    
    E_tj = (HLF_Full * P_Full_tj * delta_Full_cutout) / PLF_Full_tj * N_j_h
    return PLF_Full_tj, X_j, E_tj
