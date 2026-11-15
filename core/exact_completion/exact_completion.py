"""Exact compatible-completion verifier for the Figure 4 cap/free-quota audit.

Public inputs: seen[7,8], unresolved[8], excluded downgrade count, fixed flags,
flag cap and imaging cap. No hidden endpoint totals, truth or order enter here.
Unknown eligible endpoints can take exactly one of seven diagnostic categories.
Integer cell counts are equivalent to individual labels (no within-group label
equality is assumed). Group atomicity belongs to the caller's ledger updates.
"""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_array

RATES = np.array([95, 95, 95, 90, 90], dtype=np.int64)
CELLS = np.arange(8)
B, P, R = (CELLS & 4) != 0, (CELLS & 2) != 0, (CELLS & 1) != 0
A = B | R
UP, DOWN = ~B & P, B & ~P
MAND, MISSING = ~B & R, ~B & R & ~P
NAMES = ['mandatory_unrescuable', 'Melanoma', 'Invasive SCC', 'SCC NOS',
         'BCC', 'SCCIS', 'rescue_total', 'FP_reduction', 'flag_reduction', 'capacity']


def full_margins(counts, excluded_down, base_flags, flag_cap, K):
    """Integer-only full-record oracle, with no optimization auxiliaries."""
    x = np.asarray(counts, dtype=np.int64)
    ca = x[:6].sum(axis=0)
    nj = x[:5].sum(axis=1)
    req = (RATES * nj + 99) // 100
    tjb = np.maximum(0, req - (x[:5] * B).sum(axis=1))
    tja = np.maximum(0, req - (x[:5] * A).sum(axis=1))
    pooled = max(0, (90 * int(ca.sum()) + 99) // 100 - int((ca * A).sum()))
    t = int((ca * MAND).sum()) + max(int(tja.sum()), pooled)
    bf, rf = int((x[6] * B).sum()), int((x[6] * R).sum())
    m = max((bf + 9) // 10, bf - 9 * rf // 10)
    r = max(m, t + base_flags - flag_cap)
    cm, bp = int((ca * UP).sum()), int((x[6] * DOWN).sum())
    return np.r_[-int((ca * MISSING).sum()), (x[:5] * UP).sum(axis=1) - tjb,
                 cm - t, bp - m, bp + excluded_down - r, K - t - r].astype(np.int64)


class E:
    """Sparse integer affine expression."""
    def __init__(self, co=None, c=0):
        self.co, self.c = co or {}, int(c)

    @staticmethod
    def cast(v):
        return v if isinstance(v, E) else E(c=v)

    def __add__(self, other):
        other = E.cast(other)
        d = self.co.copy()
        for k, v in other.co.items():
            d[k] = d.get(k, 0) + v
            if not d[k]: del d[k]
        return E(d, self.c + other.c)

    __radd__ = __add__
    def __neg__(self): return self * -1
    def __sub__(self, other): return self + -E.cast(other)
    def __rsub__(self, other): return E.cast(other) + -self
    def __mul__(self, v): return E({k: a * int(v) for k, a in self.co.items()}, self.c * int(v))
    __rmul__ = __mul__


class Model:
    def __init__(self):
        self.lb, self.ub, self.rows, self.rl, self.ru = [], [], [], [], []

    def var(self, lo, hi):
        assert int(lo) <= int(hi)
        i = len(self.lb)
        self.lb.append(int(lo)); self.ub.append(int(hi))
        return E({i: 1})

    def bounds(self, expr):
        expr = E.cast(expr)
        lo = hi = expr.c
        for i, v in expr.co.items():
            lo += v * (self.lb[i] if v >= 0 else self.ub[i])
            hi += v * (self.ub[i] if v >= 0 else self.lb[i])
        return lo, hi

    def add(self, expr, lo=-np.inf, hi=np.inf):
        expr = E.cast(expr)
        self.rows.append(expr.co); self.rl.append(lo - expr.c); self.ru.append(hi - expr.c)

    def ceil(self, expr, den):
        lo, hi = self.bounds(expr)
        y = self.var(-(-lo // den), -(-hi // den))
        self.add(den * y - expr, 0, den - 1)
        return y

    def maximum(self, a, b):
        a, b = E.cast(a), E.cast(b)
        la, ua = self.bounds(a); lb, ub = self.bounds(b)
        if la >= ub: return a
        if lb >= ua: return b
        y = self.var(max(la, lb), max(ua, ub))
        bit = self.var(0, 1)
        self.add(y - a, 0); self.add(y - b, 0)
        self.add(y - a - (ub - la) * bit, hi=0)
        self.add(y - b + (ua - lb) * bit, hi=ua - lb)
        return y

    def solve(self, time_limit=10.0):
        n = len(self.lb)
        mat = np.zeros((len(self.rows), n), dtype=np.int64)
        for i, row in enumerate(self.rows):
            for j, v in row.items(): mat[i, j] = v
        res = milp(np.zeros(n), integrality=np.ones(n),
                   bounds=Bounds(self.lb, self.ub),
                   constraints=LinearConstraint(csc_array(mat.astype(float)), self.rl, self.ru),
                   options={'time_limit': time_limit, 'mip_rel_gap': 0.0, 'presolve': True})
        # A time/node limit, numerical error or unknown status never proves UNSAT.
        if res.status == 2: return 'UNSAT', None, int(res.status)
        if res.status != 0 or res.x is None: return 'UNKNOWN', None, int(res.status)
        z = np.rint(res.x).astype(np.int64)
        # Recheck every primal constraint and bound in integer arithmetic.
        if not (np.max(np.abs(z - res.x)) < 1e-5 and
                np.all(z >= self.lb) and np.all(z <= self.ub) and
                np.all(mat @ z >= self.rl) and np.all(mat @ z <= self.ru)):
            return 'UNKNOWN', None, 4
        return 'SAT', z, int(res.status)


def variables(seen, unresolved):
    seen = np.asarray(seen, dtype=np.int64)
    unresolved = np.asarray(unresolved, dtype=np.int64)
    assert seen.shape == (7, 8) and unresolved.shape == (8,)
    assert np.all(seen >= 0) and np.all(unresolved >= 0)
    model = Model()
    x = [[model.var(seen[j, c], seen[j, c] + unresolved[c]) for c in range(8)] for j in range(7)]
    for c in range(8):
        total = int(seen[:, c].sum() + unresolved[c])
        model.add(sum(x[j][c] for j in range(7)), total, total)
    return model, x


def count(x, labels=range(7), mask=np.ones(8, dtype=bool)):
    return sum((x[j][c] for j in labels for c in range(8) if mask[c]), E())


def pass_model(seen, unresolved, excluded_down, base_flags, flag_cap, K):
    """Existential feasibility via integer epigraphs; no big-M needed here.

    Auxiliary deficits may exceed their exact values, but every margin is
    monotone decreasing in these auxiliaries. Thus a satisfying epigraph model
    exists iff the exact minima/maxima satisfy all margins simultaneously.
    """
    model, x = variables(seen, unresolved)
    n = int(np.sum(seen) + np.sum(unresolved))
    model.add(count(x, range(6), MISSING), 0, 0)
    da = []
    for j, rate in enumerate(RATES):
        nj = count(x, [j])
        # C_-j >= [ceil(a_j*n_j)-TP_Bj]_+ is a direct integer inequality.
        model.add(100 * (count(x, [j], B) + count(x, [j], UP)) - int(rate) * nj, 0)
        d = model.var(0, n)
        model.add(100 * (d + count(x, [j], A)) - int(rate) * nj, 0)
        da.append(d)
    w = count(x, range(6), MAND)
    t, m, r = model.var(0, n), model.var(0, n), model.var(0, 2 * n + abs(base_flags - flag_cap))
    model.add(t - w - sum(da), 0)
    model.add(100 * (t - w + count(x, range(6), A)) - 90 * count(x, range(6)), 0)
    bf, rf = count(x, [6], B), count(x, [6], R)
    model.add(10 * m - bf, 0)
    model.add(10 * (m - bf) + 9 * rf, 0)
    model.add(r - m, 0); model.add(r - t, base_flags - flag_cap)
    model.add(count(x, range(6), UP) - t, 0)
    model.add(count(x, [6], DOWN) - m, 0)
    model.add(count(x, [6], DOWN) + excluded_down - r, 0)
    model.add(t + r, hi=K)
    return model


def adverse_model(seen, unresolved, excluded_down, base_flags, flag_cap, K, margin):
    """Exact integer ceilings and binary max graphs for one adverse margin."""
    model, x = variables(seen, unresolved)
    if margin == 0:
        model.add(count(x, range(6), MISSING), 1)
        return model
    if 1 <= margin <= 5:
        j = margin - 1
        # For integer support s >= 0, s < [ceil(a*n)-tp]_+ iff
        # 100*(s+tp) <= rate*n-1. No floating strict inequality.
        model.add(100 * (count(x, [j], UP) + count(x, [j], B)) - int(RATES[j]) * count(x, [j]), hi=-1)
        return model
    if margin in (6, 8, 9):
        deficits = []
        for j, rate in enumerate(RATES):
            req = model.ceil(int(rate) * count(x, [j]), 100)
            deficits.append(model.maximum(0, req - count(x, [j], A)))
        pooled = model.maximum(0, model.ceil(90 * count(x, range(6)), 100) - count(x, range(6), A))
        t = count(x, range(6), MAND) + model.maximum(sum(deficits), pooled)
    if margin in (7, 8, 9):
        bf, rf = count(x, [6], B), count(x, [6], R)
        m = model.maximum(model.ceil(bf, 10), bf + model.ceil(-9 * rf, 10))
    if margin in (8, 9):
        r = model.maximum(m, t + base_flags - flag_cap)
    if margin == 6: expr = count(x, range(6), UP) - t
    elif margin == 7: expr = count(x, [6], DOWN) - m
    elif margin == 8: expr = count(x, [6], DOWN) + excluded_down - r
    elif margin == 9: expr = K - t - r
    else: raise ValueError(margin)
    model.add(expr, hi=-1)
    return model


@dataclass
class Result:
    verdict: int
    possible_pass: str
    adverse_margin: int
    calls: int
    unknowns: int
    seconds: float
    pass_witness: object = None
    adverse_witness: object = None


def solve_query(seen, unresolved, excluded_down, base_flags, flag_cap, K, margin=None, time_limit=10.0):
    args = (seen, unresolved, excluded_down, base_flags, flag_cap, K)
    model = pass_model(*args) if margin is None else adverse_model(*args, margin)
    status, z, raw_status = model.solve(time_limit)
    witness = None if z is None else z[:56].reshape(7, 8)
    if witness is not None:
        assert np.all(witness >= seen)
        assert np.array_equal(witness.sum(0), np.asarray(seen).sum(0) + unresolved)
        actual = full_margins(witness, excluded_down, base_flags, flag_cap, K)
        assert (actual >= 0).all() if margin is None else actual[margin] < 0
    return status, witness, raw_status


def verify(seen, unresolved, excluded_down, base_flags, flag_cap, K, time_limit=10.0):
    start = perf_counter()
    args = (seen, unresolved, excluded_down, base_flags, flag_cap, K)
    passed, good, _ = solve_query(*args, time_limit=time_limit)
    calls, unknowns = 1, int(passed == 'UNKNOWN')
    if passed == 'UNSAT':
        return Result(-1, passed, -1, calls, unknowns, perf_counter() - start)
    # Conservatively defer immediately if feasible completion is not established.
    if passed == 'UNKNOWN':
        return Result(0, passed, -1, calls, unknowns, perf_counter() - start)
    for j in range(10):
        adverse, bad, _ = solve_query(*args, margin=j, time_limit=time_limit)
        calls += 1; unknowns += int(adverse == 'UNKNOWN')
        if adverse != 'UNSAT':
            return Result(0, passed, j, calls, unknowns, perf_counter() - start, good, bad)
    return Result(1, passed, -1, calls, unknowns, perf_counter() - start, good)
