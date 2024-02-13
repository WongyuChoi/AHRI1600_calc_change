# AHRI 1600: 11.2.1.2
# functions: Two-stage System and Variable capacity Certified, Two-Stage Capacity System
from typing import List, Union, Tuple

def get_q_tj(
    tj_i: float,
    N_j_c: float,
    BL_tj: float,
    q_dot_Low_tj: float,
    q_dot_Full_tj: float,
    lock_out_low_capacity_ops: bool,
    od_temp_when_lock_out: Union[float, None],
) -> Tuple[float, float, float, float]:
  """
  get cooling energy (Btu) at each bin temperature
  In R version, ratioTotalCooling
  :param tj_i: bin temperature
  :param BL_tj: building load at each bin temperature
  :param q_dot_Low_tj: Low stage capacity at each bin temperature
  :param q_dot_Full_tj: Full stage capacity at each bin temperature
  :param lock_out_low_capacity_ops: lock out low capacity operations
  :param od_temp_when_lock_out: outdoor temperature when lock out low capacity operations if None, just lock out low capacity operations for that range
  :param N_j_c: fraction of bin hours
  :return: CLF (or X_j) cooling energy at each bin temperature
  """
  X_Tj_k1 = 0.0
  X_Tj_k2 = 0.0
  X_j = 0.0 # on-time fractions for the lowest compressor stage for bin j.
  ratio_total_cool = 0.0

  if q_dot_Low_tj >= BL_tj:
    X_Tj_k1 = BL_tj / q_dot_Low_tj
    X_j = X_Tj_k1
    ratio_total_cool = X_Tj_k1 * q_dot_Low_tj * N_j_c

  elif q_dot_Low_tj < BL_tj < q_dot_Full_tj:
    X_Tj_k1 = (q_dot_Full_tj - BL_tj) / (q_dot_Full_tj - q_dot_Low_tj)
    X_Tj_k2 = 1 - X_Tj_k1
    # X_j = X_Tj_k1
    X_j = 1
    ratio_total_cool = (X_Tj_k1 * q_dot_Low_tj + X_Tj_k2 * q_dot_Full_tj) * N_j_c

  if (q_dot_Low_tj < BL_tj < q_dot_Full_tj) and lock_out_low_capacity_ops:
    if od_temp_when_lock_out: # when od_temp_when_lock_out is not None
      if tj_i >= od_temp_when_lock_out: # lock out temperature is reached or exceeded, alternating full and off
        X_Tj_k2 = BL_tj / q_dot_Full_tj
        X_j = X_Tj_k2
        ratio_total_cool = X_Tj_k2 * q_dot_Full_tj * N_j_c
      else: # lock out temperature is not reached, just as the above  q_dot_Low_tj < BL_tj < q_dot_Full_tj
        X_Tj_k1 = (q_dot_Full_tj - BL_tj) / (q_dot_Full_tj - q_dot_Low_tj)
        X_Tj_k2 = 1 - X_Tj_k1
        # X_j = X_Tj_k1
        X_j = 1
        ratio_total_cool = (X_Tj_k1 * q_dot_Low_tj + X_Tj_k2 * q_dot_Full_tj) * N_j_c
    else: # when od_temp_when_lock_out is None, just lock out low capacity operations
      X_Tj_k2 = BL_tj / q_dot_Full_tj
      X_j = X_Tj_k2
      ratio_total_cool = X_Tj_k2 * q_dot_Full_tj * N_j_c

  if BL_tj >= q_dot_Full_tj:
    X_Tj_k2 = 1.0
    X_j = 1.0 # full load always on
    ratio_total_cool = q_dot_Full_tj * N_j_c

  return X_Tj_k1, X_Tj_k2, X_j, ratio_total_cool


