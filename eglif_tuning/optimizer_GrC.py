from cerebellar_models.optimization.optimizer import Optimizer, Constraint, _apply_constraints, ArchiveND
import random, os, warnings, logging
import numpy as np
import pandas as pd
from cerebellar_models.optimization.features import *

# === WARNING SETTING ===
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
logging.getLogger().setLevel(logging.ERROR)
for name in ("deap", "matplotlib", "numexpr"):
    logging.getLogger(name).setLevel(logging.ERROR)
np.seterr(all="ignore")
os.environ.setdefault("PYTHONWARNINGS", "ignore")


def constrain_bounds(ind, ctx):
    bounds = ctx["bounds"]
    for pos, (a, b) in bounds.items():
        x = float(ind[pos]) if np.isfinite(ind[pos]) else (a + b) / 2.0
        if a > b:
            a, b = b, a
        ind[pos] = min(max(x, a), b)
    return ind


def kadap_constrain(ind, ctx):
    """
    k_adap >= f(k2)
    f(k2) = (C_m/4) * (k2 + 1/tau_m)^2
    """
    idx_kadap = ctx["idx_kadap"]
    C_m = ctx["C_m"]
    tau_m = ctx["tau_m"]
    k2 = 1.0 / tau_m
    kad_low = (C_m / 4.0) * (k2 + 1.0 / tau_m) ** 2
    if ind[idx_kadap] < kad_low:
        ind[idx_kadap] = kad_low
    return ind


def random_init(optimizer):
    ind = [random.uniform(*optimizer.bounds[p]) for p in optimizer.opt_params]
    ind = _apply_constraints(ind, optimizer.constraints)
    return ind


def _scalar(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return x[0] if len(x) else np.nan
    return x


def _smape(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.abs(a) + np.abs(b)
    num = 2.0 * np.abs(a - b)
    return np.where(denom > 0, num / denom, 0.0)


def _weighted_median(values, weights):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if v.size == 0:
        return np.nan
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    cutoff = 0.5 * cw[-1]
    return float(v[np.searchsorted(cw, cutoff)])


def fit_fi_slope(df, threshold, selected_currents):
    sub = df[df['current'].isin(selected_currents)].sort_values('current')
    x = (sub['current'].values - threshold).astype(float)
    y = sub['mean_frequency'].apply(_scalar).astype(float).values
    m = (x > 0) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size == 0:
        return np.nan
    r = y / x
    w = np.sqrt(x)
    return _weighted_median(r, w)


def _cvar90(x):
    """Conditional Value at Risk al 90° percentile (più robusto della media)."""
    arr = np.sort(np.asarray(x, dtype=float))
    if arr.size == 0:
        return np.nan
    k = int(0.1 * arr.size)
    return float(np.mean(arr[-k:])) if k > 0 else float(np.mean(arr))


def _seg_slopes(x, y):
    """Slopes tra punti consecutivi di una f–I curve."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx = np.diff(x)
    dy = np.diff(y)
    m = dx > 0
    return np.divide(dy[m], dx[m], out=np.full_like(dy[m], np.nan), where=m)


def evaluate_GrC(optimizer, ind, cell_params=None):
    params = dict(cell_params if cell_params is not None else optimizer.model_params)
    for i, p in enumerate(optimizer.opt_params):
        params[p] = float(ind[i])
    optimizer.model_params = params

    # Feature extraction
    target_metrics = optimizer._extract_multicomp_features()
    nest_metrics = optimizer._extract_nest_features()

    # 1. Rheobase error
    valid_thr = nest_metrics.loc[nest_metrics['spike_count_stimint'] > 0, 'current']
    nest_thr = valid_thr.min() if not valid_thr.empty else None
    valid_thr_t = target_metrics.loc[target_metrics['spike_count_stimint'] > 0, 'current']
    targ_thr = valid_thr_t.min() if not valid_thr_t.empty else None

    if (nest_thr is None) or (targ_thr is None) or (not np.isfinite(targ_thr)):
        thr_loss = 1.0
    else:
        rel = abs(nest_thr - targ_thr) / max(abs(targ_thr), 1e-9)
        thr_loss = float(np.clip(rel / (1.0 + rel), 0.0, 1.0))

    # 2. Slope core
    uniq_currents = np.sort(target_metrics['current'].unique())
    dI = float(np.median(np.diff(uniq_currents))) if len(uniq_currents) > 1 else 0.0
    cut = (targ_thr if np.isfinite(targ_thr) else -np.inf) + dI
    selected_currents = target_metrics.loc[target_metrics['current'] >= cut, 'current'].values

    slope_targ = fit_fi_slope(target_metrics, targ_thr, selected_currents)
    slope_nest = fit_fi_slope(nest_metrics, targ_thr, selected_currents)
    if (not np.isfinite(slope_targ)) or (slope_targ == 0) or (not np.isfinite(slope_nest)):
        slope_core = 1.0
    else:
        r = abs(slope_nest - slope_targ) / max(abs(slope_targ), 1e-9)
        slope_core = float(np.clip(r / (1.0 + r), 0.0, 1.0))

    # 3. Shape guardrail (+ tail slope pen)
    nm_sel = nest_metrics[nest_metrics['current'].isin(selected_currents)].sort_values('current')
    cur = nm_sel['current'].values
    nf = nm_sel['mean_frequency'].apply(_scalar).astype(float).values
    thr_ref = nest_thr if (nest_thr is not None and np.isfinite(nest_thr)) else targ_thr
    slope_own = fit_fi_slope(nest_metrics, thr_ref, selected_currents)

    if not np.isfinite(slope_own) or nf.size == 0:
        shape_cvar = spiky_last = cv_pen = early_loss = tail_slope_pen = 1.0
    else:
        x = np.maximum(0.0, cur - thr_ref)
        y_hat = slope_own * x
        sm = _smape(nf, y_hat)

        shape_cvar = _cvar90(sm)
        shape_cvar = 1.0 if not np.isfinite(shape_cvar) else float(np.clip(shape_cvar, 0, 1))

        mpos = x > 0
        xs, ys = x[mpos], nf[mpos]
        seg = _seg_slopes(xs, ys)
        if seg.size >= 2 and np.isfinite(seg).all():
            med = float(np.median(seg))
            last = float(seg[-1])
            spiky_last = 1.0 if med <= 0 else float(np.clip(max(0.0, (last - med) / med) / (1.0 + max(0.0, (last - med) / med)), 0, 1))
            mu = float(np.mean(seg))
            sd = float(np.std(seg))
            cv_pen = 1.0 if mu <= 0 else float(np.clip(abs(sd / mu) / (1.0 + abs(sd / mu)), 0, 1))
        else:
            spiky_last = cv_pen = 1.0

        early_loss = 1.0
        try:
            i1 = float(np.min(selected_currents))
            tm1 = target_metrics.loc[target_metrics['current'] == i1].iloc[0]
            nm1 = nest_metrics.loc[nest_metrics['current'] == i1].iloc[0]
            lat_t = float(_scalar(tm1['time_to_first_spike']))
            lat_n = float(_scalar(nm1['time_to_first_spike']))
            if np.isfinite(lat_t) and lat_t > 0 and np.isfinite(lat_n):
                rlat = abs(lat_n - lat_t) / max(lat_t, 1e-3)
                lat_loss = rlat / (1.0 + rlat)
            else:
                lat_loss = 1.0
            inv1_t = float(_scalar(tm1['inv_first_ISI']))
            inv1_n = float(_scalar(nm1['inv_first_ISI']))
            sm1 = _smape([inv1_n], [inv1_t])[0] if np.isfinite(inv1_t) and np.isfinite(inv1_n) else 1.0
            early_loss = float(np.clip(0.5 * lat_loss + 0.5 * sm1, 0.0, 1.0))
        except Exception:
            early_loss = 1.0

        try:
            tail_mask = cur >= (np.nanmax(cur) - 3 * max(dI, 1e-9))
            xt, yt = x[tail_mask], nf[tail_mask]
            seg_tail = _seg_slopes(xt, yt)
            if seg_tail.size >= 1 and np.isfinite(slope_targ):
                slope_tail = float(np.median(seg_tail))
                if slope_tail > slope_targ:
                    rtail = (slope_tail - slope_targ) / max(abs(slope_targ), 1e-9)
                    tail_slope_pen = float(np.clip(rtail / (1.0 + rtail), 0.0, 1.0))
                else:
                    tail_slope_pen = 0.0
            else:
                tail_slope_pen = 1.0
        except Exception:
            tail_slope_pen = 1.0

    shape_loss = float(np.clip(max(shape_cvar, spiky_last, cv_pen, early_loss, tail_slope_pen), 0.0, 1.0))

    # 4. Gap loss
    tm_sel = target_metrics[target_metrics['current'].isin(selected_currents)][['current', 'mean_frequency']].sort_values('current').copy()
    tm_sel['tf'] = tm_sel['mean_frequency'].apply(_scalar).astype(float)

    nm_join = nest_metrics[nest_metrics['current'].isin(selected_currents)][['current', 'mean_frequency']].sort_values('current').copy()
    nm_join['nf'] = nm_join['mean_frequency'].apply(_scalar).astype(float)

    comp = pd.merge(nm_join[['current', 'nf']], tm_sel[['current', 'tf']], on='current', how='inner')
    if comp.shape[0] == 0:
        gap_loss = 1.0
    else:
        I = comp['current'].values.astype(float)
        tf = comp['tf'].values.astype(float)
        nf = comp['nf'].values.astype(float)

        tau_thr = max(2.0 * dI, 1e-9)
        tau_tail = max(2.0 * dI, 1e-9)
        Imax = float(np.nanmax(I))

        w_thr = np.exp(-(I - targ_thr) / tau_thr)
        w_thr[I <= targ_thr] = 0.0
        w_tail = np.exp(-(Imax - I) / tau_tail)
        w_tail[I < (Imax - 3 * tau_tail)] *= 0.25

        alpha = 0.55
        w_mix = alpha * w_thr + (1.0 - alpha) * w_tail

        sm = _smape(nf, tf)
        under = (nf < tf).astype(float)
        over = (nf > tf).astype(float)

        asym_thr = 1.5 * under * np.exp(-(I - targ_thr) / tau_thr) + 1.0
        asym_tail = 1.5 * over * np.exp(-(Imax - I) / tau_tail) + 1.0
        asym = alpha * asym_thr + (1.0 - alpha) * asym_tail

        sm_weighted = sm * asym
        if (not np.isfinite(w_mix).any()) or np.nansum(w_mix) == 0:
            gap_loss = 1.0
        else:
            gap_loss = float(np.clip(np.nansum((w_mix / np.nansum(w_mix)) * sm_weighted), 0, 1))

    return (thr_loss, slope_core, shape_loss, gap_loss)


if __name__ == "__main__":
    # Example for GrC
    data_folder = os.path.join(os.path.dirname(__file__), 'tofit_eglif/results_tofitEglif/GrC/')
    threshold_GrC = -41.0
    protocol = {
        'start_stim': 100.0,
        'end_stim': 600.0,
        'duration': 700.0,
        'threshold': threshold_GrC,
    }
    cell_params = {
        "t_ref": 1.5, "V_min": -150, "C_m": 7, "V_th": -41, "V_reset": -70,
        "E_L": -62, "I_e": -0.888, "V_m": -62.0, "tau_m": 24.15,
        "k_adap": 0.022, "k_1": 0.311, "k_2": 0.041407868,
        "A1": 0.01, "A2": -0.94,
    }
    cell_params['k_2'] = 1.0 / cell_params['tau_m']  # oscillatory regime

    bounds = {
        "I_e": (-10.0, 10.0),
        "k_adap": (-2, 1),
        "k_1": (0.01, 10.0),
        "A1": (0.01, 50.0),
        "A2": (-10.0, 30.0),
    }

    fitness = {
        'rheobase_error': -1,
        'slope_error': -1,
        'shape_error': -1,
        'gap_error': -1,
    }

    archive = ArchiveND(
        eps_vec=[0.02, 0.05, 0.05, 0.05],
        tau_vec=[0.08, 0.15, 0.20, 0.12],
        cap=1000
    )

    optimizer = Optimizer(
        multicomp_data=data_folder,
        protocol=protocol,
        nest_model='eglif_multirec_opt',
        opt_params=list(bounds.keys()),
        model_params=cell_params,
        fitness=fitness,
        bounds=bounds,
        archive=archive,
        knee_weights=[2.0, 1.0, 1.0, 1.0]
    )

    param_positions = {name: i for i, name in enumerate(optimizer.opt_params)}
    constraint1 = Constraint(
        func=constrain_bounds,
        name="bounds",
        ctx={"bounds": {param_positions[k]: v for k, v in bounds.items()}}
    )
    constraint2 = Constraint(
        func=kadap_constrain,
        name="kadap>=f(k2)",
        ctx={
            "idx_kadap": param_positions["k_adap"],
            "C_m": cell_params["C_m"],
            "tau_m": cell_params["tau_m"]
        }
    )
    optimizer.add_constraint(constraint1)
    optimizer.add_constraint(constraint2)

    optimizer.init_params_fn = random_init
    optimizer.evaluate_fn = evaluate_GrC
    optimizer.print_evolution = True
    toolbox = optimizer.set_optimizer()

    best_individual, best_fitness = optimizer.optimize()
    print(best_individual)
