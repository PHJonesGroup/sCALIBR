import pandas as pd
import numpy as np
from .p_control_any_implement import p_control_any_implement
from .make_tables_Z_two import make_tables_Z_two
from .computeZ import computeZ

def z_p_CTR_any(
    st: float,
    en: float,
    step: float,
    T_any: pd.DataFrame,
    binn: np.ndarray,
    p_cont: np.ndarray,
    me_sd: np.ndarray,
    cond1: str,
    cond2: str,
):
    """
    Calibrate any control set (e.g. intergenic or NT) against zGE controls,
    assign p/q values, compute Z‑scores, and return per‑replicate tables + histos.
    """
    T_ij_q_verti_any = p_control_any_implement(T_any, binn, p_cont)

    LFC_t = T_ij_q_verti_any.iloc[:, 2].to_numpy()

    Z_t1, binn, perc_t1, perc_zt1 = computeZ(
        st, en, step, LFC_t, me_sd
    )

    T_z_q_any1 = make_tables_Z_two(
        T_ij_q_verti_any, Z_t1
    )

    return T_z_q_any1, binn, perc_t1, perc_zt1
