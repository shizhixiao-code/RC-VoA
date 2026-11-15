# Unmodified functions captured from archived dual_verification.py.
import numpy as np
from exact_completion import RATES as RATE, B as BASE, R as REF, A as UNION, UP, DOWN, MAND as MANDATORY, MISSING

def floor_bounds(seen, unresolved, decision):
    """Bounds for deficits relative to any public binary baseline decision."""
    n_j = seen[..., :5, :].sum(axis=-1)
    tp_j = (seen[..., :5, :]*decision).sum(axis=-1)
    u_pos = (unresolved*decision).sum(axis=-1)
    u_neg = (unresolved*~decision).sum(axis=-1)
    j_lo = np.maximum(0, (RATE*(n_j+u_pos[...,None])+99)//100-tp_j-u_pos[...,None])
    j_hi = np.maximum(0, (RATE*(n_j+u_neg[...,None])+99)//100-tp_j)
    n_c = seen[..., :6, :].sum(axis=(-1,-2))
    tp_c = (seen[..., :6, :]*decision).sum(axis=(-1,-2))
    p_lo = np.maximum(0, (90*(n_c+u_pos)+99)//100-tp_c-u_pos)
    p_hi = np.maximum(0, (90*(n_c+u_neg)+99)//100-tp_c)
    t_lo = np.maximum(j_lo.sum(axis=-1), p_lo)
    t_hi = np.minimum(np.maximum(j_hi.sum(axis=-1), p_hi), n_c-tp_c+u_neg)
    return j_lo, j_hi, t_lo, t_hi

def certificate(seen, unresolved, excluded_down, base_flags, flag_cap, K):
    seen = np.asarray(seen, dtype=np.int64)
    unresolved = np.asarray(unresolved, dtype=np.int64)
    assert seen.shape[-2:] == (7,8) and unresolved.shape[-1] == 8
    assert np.all(seen >= 0) and np.all(unresolved >= 0)
    ca = seen[..., :6, :].sum(axis=-2)
    j_lo, j_hi, _, _ = floor_bounds(seen, unresolved, BASE)
    _, _, union_lo, union_hi = floor_bounds(seen, unresolved, UNION)
    w_lo = (ca*MANDATORY).sum(axis=-1)
    w_hi = w_lo+(unresolved*MANDATORY).sum(axis=-1)
    t_lo = w_lo+union_lo
    t_hi = np.minimum(w_hi+union_hi, ((ca+unresolved)*~BASE).sum(axis=-1))
    missing_lo = (ca*MISSING).sum(axis=-1)
    missing_hi = missing_lo+(unresolved*MISSING).sum(axis=-1)
    bf_lo = (seen[..., 6, :]*BASE).sum(axis=-1)
    bf_hi = bf_lo+(unresolved*BASE).sum(axis=-1)
    sf_lo = (seen[..., 6, :]*REF).sum(axis=-1)
    sf_hi = sf_lo+(unresolved*REF).sum(axis=-1)
    m_lo = np.maximum((bf_lo+9)//10, bf_lo-9*sf_hi//10)
    m_hi = np.maximum((bf_hi+9)//10, bf_hi-9*sf_lo//10)
    delta = int(base_flags)-int(flag_cap)
    r_lo, r_hi = np.maximum(m_lo,t_lo+delta), np.maximum(m_hi,t_hi+delta)
    cmj_lo = (seen[..., :5, :]*UP).sum(axis=-1)
    cmj_hi = cmj_lo+(unresolved*UP).sum(axis=-1)[...,None]
    cm_lo = (ca*UP).sum(axis=-1)
    cm_hi = cm_lo+(unresolved*UP).sum(axis=-1)
    bp_lo = (seen[...,6,:]*DOWN).sum(axis=-1)
    bp_hi = bp_lo+(unresolved*DOWN).sum(axis=-1)
    lower = np.concatenate([-missing_hi[...,None],cmj_lo-j_hi,
        np.stack([cm_lo-t_hi,bp_lo-m_hi,bp_lo+excluded_down-r_hi,K-t_hi-r_hi],axis=-1)],axis=-1)
    upper = np.concatenate([-missing_lo[...,None],cmj_hi-j_lo,
        np.stack([cm_hi-t_lo,bp_hi-m_lo,bp_hi+excluded_down-r_lo,K-t_lo-r_lo],axis=-1)],axis=-1)
    assert np.all(lower <= upper)
    passed, failed = (lower >= 0).all(axis=-1), (upper < 0).any(axis=-1)
    assert not np.any(passed & failed)
    verdict = np.where(passed,1,np.where(failed,-1,0)).astype(np.int8)
    return verdict, lower, upper
