"""
Set of evaluation function for the optimization of NEST single cell parameters vs NEURON current injection experiments.
"""

import numpy as np
import pandas as pd

from cerebellar_models.optimization.features import inv_first_isi, inv_second_isi


# -------   UTILS   -------
def _scalar(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return x[0] if len(x) else np.nan
    return x


def currents_above_thr(df, thr):
    """
    Select only experiments of current injection above a given threshold.
    """
    uniq = np.sort(df["current"].unique())
    dI = float(np.median(np.diff(uniq))) if len(uniq) > 1 else 0.0
    return df.loc[df["current"] >= (thr + dI), "current"].values


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


def _fit_fi_slope(df, threshold, selected_currents):
    sub = df[df["current"].isin(selected_currents)].sort_values("current")
    x = (sub["current"].values - threshold).astype(float)
    y = sub["mean_frequency"].apply(_scalar).astype(float).values
    m = (x > 0) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size == 0:
        return np.nan
    r = y / x
    w = np.sqrt(x)
    return _weighted_median(r, w)


def spike_penalty(diff_spikes, scale=10):
    return min(1.0, float(diff_spikes) / float(scale))


def _fi_curve(df, sel=None):

    sub = df if sel is None else df[df["current"].isin(sel)]
    sub = sub.sort_values("current")

    I = sub["current"].values.astype(float)
    f = sub["mean_frequency"].apply(_scalar).astype(float).values

    m = np.isfinite(I) & np.isfinite(f)
    return I[m], f[m]


def extract_info(df, start, end):
    out = []
    for curr in df["current"].unique():
        spikes = np.array(df.loc[df["current"] == curr, "peak_time"].values[0])
        spikes = spikes[(spikes >= start) & (spikes <= end)]

        if spikes.size > 0:
            first_spike = spikes.min()
            n_spikes = len(spikes)
            if spikes.size > 1:
                second_spike = spikes[1]
            else:
                second_spike = np.nan
        else:
            first_spike = np.nan
            second_spike = np.nan
            n_spikes = 0

        out.append(
            {
                "current": curr,
                "first_spike": first_spike,
                "second_spike": second_spike,
                "n_spikes": n_spikes,
            }
        )
    return pd.DataFrame(out)


# -------   LOSSES   -------
def rheobase_loss(nest_df, targ_df):
    try:
        nest_thr = nest_df.loc[nest_df["spike_count_stimint"] > 0, "current"].min()
        targ_thr = targ_df.loc[targ_df["spike_count_stimint"] > 0, "current"].min()
        rel_thr = abs(nest_thr - targ_thr) / max(abs(targ_thr), 1e-9)
        return np.clip(rel_thr / (1 + rel_thr), 0, 1), targ_thr
    except Exception:
        return 1.0, np.nan


def pacemaking_loss(targ, nest):
    try:
        tf = float(_scalar(targ.loc[targ["current"] == 0, "mean_frequency"]))
        nf = float(_scalar(nest.loc[nest["current"] == 0, "mean_frequency"]))
        rel = abs(nf - tf) / max(tf, 1e-9)
        return np.clip(rel / (1 + rel), 0, 1)
    except Exception:
        return 1.0


def cv_loss(targ, nest, delta=0.1, alpha=0.3, regularity=False):
    try:
        tcv = float(_scalar(targ.loc[targ["current"] == 0, "ISI_CV"]))
        ncv = float(_scalar(nest.loc[nest["current"] == 0, "ISI_CV"]))
        diff = abs(ncv - tcv)

        if diff <= delta:
            return (diff / delta) * alpha

        loss_base = diff / (diff + delta)
        if regularity:
            if ncv < tcv:
                return alpha * loss_base
            else:
                return loss_base
        return loss_base

    except Exception:
        return 1.0


def slope_loss(targ, nest, thr, sel):
    try:
        st, sn = _fit_fi_slope(targ, thr, sel), _fit_fi_slope(nest, thr, sel)
        r = abs(sn - st) / max(abs(st), 1e-9)
        return np.clip(r / (1 + r), 0, 1)
    except Exception:
        return 1.0


def gap_loss(targ, nest, sel, weighted=None, missing_penalty=1.0):
    try:
        tm = targ[targ["current"].isin(sel)][["current", "mean_frequency"]]
        nm = nest[nest["current"].isin(sel)][["current", "mean_frequency"]]

        comp = pd.merge(
            nm.assign(nf=nm["mean_frequency"].apply(_scalar)),
            tm.assign(tf=tm["mean_frequency"].apply(_scalar)),
            on="current",
        )

        if comp.empty:
            return float(missing_penalty)

        I = comp["current"].values.astype(float)
        tf = comp["tf"].values.astype(float)
        nf = comp["nf"].values.astype(float)

        m = np.isfinite(I) & np.isfinite(tf) & np.isfinite(nf)
        I, tf, nf = I[m], tf[m], nf[m]

        if I.size == 0:
            return float(missing_penalty)

        diff = np.abs(nf - tf)

        max_distance = np.max(tf) - np.min(tf)
        if max_distance > 0:
            err_vec = diff / max_distance
        else:
            err_vec = np.zeros_like(diff)

        err_vec = np.clip(err_vec, 0.0, 1.0)

        if weighted is not None:
            if weighted == "gaussian":
                Im = 0.5 * (I.min() + I.max())
                s = 0.35 * (I.max() - I.min() + 1e-9)
                w = np.exp(-((I - Im) ** 2) / (2.0 * s**2))

            elif weighted == "inverse":
                w = 1.0 / (I + 1e-9)

            elif weighted == "high":
                I_shift = I - I.min()
                if np.allclose(I_shift, 0.0):
                    w = np.ones_like(I_shift)
                else:
                    w = I_shift
            else:
                raise ValueError(f"Weight strategy {weighted} not implemented")

            w = np.asarray(w, dtype=float)
            w[m]
            if np.all(w <= 0):
                mean_err = float(np.mean(err_vec))
            else:
                mean_err = float(np.sum(w * err_vec) / np.sum(w))
        else:
            mean_err = float(np.mean(err_vec))

        return float(np.clip(mean_err, 0.0, 1.0))

    except Exception:
        return float(missing_penalty)


def curvature_loss(
    targ,
    nest,
    sel=None,
    missing_penalty=0.5,
    use_max_slope_penalty=False,
    max_slope_weight=0.5,
):
    try:
        It, Ft = _fi_curve(targ, sel)
        In, Fn = _fi_curve(nest, sel)

        common = np.intersect1d(It, In)
        if common.size < 3:
            return float(missing_penalty)

        mt = np.isin(It, common)
        mn = np.isin(In, common)
        It, Ft = It[mt], Ft[mt]
        In, Fn = In[mn], Fn[mn]

        if not np.allclose(It, In):
            return float(missing_penalty)

        I = It
        dI = np.diff(I)
        st = np.diff(Ft) / dI
        sn = np.diff(Fn) / dI

        if st.size < 2:
            return float(missing_penalty)

        ct = np.diff(st)
        cn = np.diff(sn)

        scale = np.max(np.abs(ct))
        if scale <= 0:
            scale = 1.0

        diff_curv = np.abs(cn - ct)
        mean_rel_curv = float(np.mean(diff_curv) / scale)
        curv_loss = mean_rel_curv / (1.0 + mean_rel_curv)
        curv_loss = float(np.clip(curv_loss, 0.0, 1.0))

        if use_max_slope_penalty and max_slope_weight > 0.0:
            smax_t = float(np.max(np.abs(st)))
            smax_n = float(np.max(np.abs(sn)))
            denom = max(smax_t, 1e-9)
            rel_slope_err = abs(smax_n - smax_t) / denom
            slope_pen = rel_slope_err / (1.0 + rel_slope_err)
            slope_pen = float(np.clip(slope_pen, 0.0, 1.0))
            total_loss = curv_loss + max_slope_weight * slope_pen
            total_loss = float(np.clip(total_loss, 0.0, 1.0))
            return total_loss

        return curv_loss

    except Exception:
        return float(missing_penalty)


def post_first_spike_loss(
    targ, nest, protocol, thr=0.0, sign="pos", missing_penalty=0.5, mode="max"
):
    start, end = protocol["end_stim"], protocol["duration"]

    if getattr(targ, "empty", True) or getattr(nest, "empty", True):
        return float(missing_penalty)

    common = np.intersect1d(targ["current"].unique(), nest["current"].unique())
    if common.size == 0:
        return float(missing_penalty)

    tm = extract_info(targ[targ["current"].isin(common)], start, end)
    nm = extract_info(nest[nest["current"].isin(common)], start, end)

    if sign == "pos":
        tm_sub = tm[tm["current"] >= thr]
        nm_sub = nm[nm["current"] >= thr]
    elif sign == "neg":
        tm_sub = tm[tm["current"] < thr]
        nm_sub = nm[nm["current"] < thr]
    else:
        raise ValueError("sign must be 'pos' or 'neg'")

    losses = []
    for c in tm_sub["current"]:
        t = tm_sub[tm_sub["current"] == c].iloc[0]
        n = nm_sub[nm_sub["current"] == c].iloc[0] if c in nm_sub["current"].values else None

        if t["n_spikes"] == 0:
            if n is None or n["n_spikes"] == 0:
                losses.append(0.0)
            else:
                losses.append(spike_penalty(n["n_spikes"]))
        else:
            if n is None or n["n_spikes"] == 0:
                losses.append(1.0)
            else:
                err = abs(n["first_spike"] - t["first_spike"]) / max(t["first_spike"], 1e-9)
                losses.append(np.clip(err / (1 + err), 0.0, 1.0))
    if mode == "max":
        return float(np.max(losses)) if losses else float(missing_penalty)
    return float(np.mean(losses)) if losses else float(missing_penalty)


def poststim_autorhythm_loss(
    targ_post, nest_post, protocol, thr, transient=300.0, missing_penalty=0.5
):
    start = protocol["end_stim"] + transient
    end = protocol["duration"]
    window = max(end - start, 1e-9)

    if getattr(targ_post, "empty", True) or getattr(nest_post, "empty", True):
        return float(missing_penalty)

    common = np.intersect1d(targ_post["current"].unique(), nest_post["current"].unique())
    if common.size == 0:
        return float(missing_penalty)
    common = common[common < thr]
    if common.size == 0:
        return float(missing_penalty)

    tm = extract_info(targ_post[targ_post["current"].isin(common)], start, end)
    nm = extract_info(nest_post[nest_post["current"].isin(common)], start, end)

    losses = []
    for c in tm["current"]:
        t = tm[tm["current"] == c].iloc[0]
        n = nm[nm["current"] == c].iloc[0] if c in nm["current"].values else None

        ft = t["n_spikes"] / window
        fn = (n["n_spikes"] / window) if n is not None else 0.0

        if t["n_spikes"] == 0:
            continue

        if n is None or n["n_spikes"] == 0:
            losses.append(1.0)
        else:
            rel = abs(fn - ft) / max(ft, 1e-9)
            losses.append(float(np.clip(rel / (1.0 + rel), 0.0, 1.0)))

    if not losses:
        return float(missing_penalty)

    return float(np.nanmean(losses))


def post_rebound_loss(
    targ, nest, protocol, window=50.0, thr=0.0, sign="pos", missing_penalty=0.5, scale=10
):
    win_start = protocol["end_stim"]
    win_end = win_start + window
    if getattr(targ, "empty", True) or getattr(nest, "empty", True):
        return float(missing_penalty)

    common = np.intersect1d(targ["current"].unique(), nest["current"].unique())
    if common.size == 0:
        return float(missing_penalty)

    tm = extract_info(targ[targ["current"].isin(common)], win_start, win_end)
    nm = extract_info(nest[nest["current"].isin(common)], win_start, win_end)

    if sign == "pos":
        tm = tm[tm["current"] >= thr]
        nm = nm[nm["current"] >= thr]
    elif sign == "neg":
        tm = tm[tm["current"] < thr]
        nm = nm[nm["current"] < thr]
    else:
        raise ValueError("sign must be 'pos' or 'neg'")

    losses = []

    for c in tm["current"]:
        t = tm[tm["current"] == c].iloc[0]
        n = nm[nm["current"] == c].iloc[0] if c in nm["current"].values else None

        if t["n_spikes"] == 0:
            if n is None or n["n_spikes"] == 0:
                losses.append(0.0)
            else:
                losses.append(spike_penalty(n["n_spikes"], scale))
            continue

        if n is None or n["n_spikes"] == 0:
            losses.append(1.0)
            continue

        diff = abs(n["n_spikes"] - t["n_spikes"])
        losses.append(spike_penalty(diff, scale))

    return float(np.nanmean(losses)) if losses else float(missing_penalty)
