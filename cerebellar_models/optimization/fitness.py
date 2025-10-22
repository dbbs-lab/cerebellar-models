"""
Set of evaluation function for the optimization of NEST single cell parameters vs NEURON current injection experiments.
"""

import numpy as np
import pandas as pd


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


def _smape(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.abs(a) + np.abs(b)
    num = 2.0 * np.abs(a - b)
    return np.where(denom > 0, num / denom, 0.0)


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


def cv_loss(targ, nest):
    try:
        tcv = float(_scalar(targ.loc[targ["current"] == 0, "ISI_CV"]))
        ncv = float(_scalar(nest.loc[nest["current"] == 0, "ISI_CV"]))
        rel = abs(ncv - tcv) / max(tcv, 1e-9)
        return np.clip(rel / (1 + rel), 0, 1)
    except Exception:
        return 1.0


def slope_loss(targ, nest, thr, sel):
    try:
        st, sn = _fit_fi_slope(targ, thr, sel), _fit_fi_slope(nest, thr, sel)
        r = abs(sn - st) / max(abs(st), 1e-9)
        return np.clip(r / (1 + r), 0, 1)
    except Exception:
        return 1.0


def gap_loss(targ, nest, sel, weighted=None):
    try:
        tm = targ[targ["current"].isin(sel)][["current", "mean_frequency"]]
        nm = nest[nest["current"].isin(sel)][["current", "mean_frequency"]]

        comp = pd.merge(
            nm.assign(nf=nm["mean_frequency"].apply(_scalar)),
            tm.assign(tf=tm["mean_frequency"].apply(_scalar)),
            on="current",
        )

        if comp.empty:
            return 1.0

        if weighted:
            I = comp["current"].values
            if weighted == "gaussian":
                Im, s = 0.5 * (I.min() + I.max()), 0.35 * (I.max() - I.min())
                w = np.exp(-((I - Im) ** 2) / (2 * s**2))
            elif weighted == "inverse":
                w = 1 / (I + 1e-10)
            else:
                raise ValueError(f"Weight strategy {weighted} not implemented")

            smape_loss = _smape(comp["nf"], comp["tf"])
            weighted_smape = np.sum(w * smape_loss) / np.sum(w)
            return np.clip(weighted_smape, 0, 1)

        smape_loss = _smape(comp["nf"], comp["tf"])
        return np.clip(np.nanmean(smape_loss), 0, 1)
    except Exception:
        return 1.0


def poststim_latency_loss(targ_post, nest_post, protocol, missing_penalty=0.5):
    start, end = protocol["end_stim"], protocol["duration"]
    if getattr(targ_post, "empty", True) or getattr(nest_post, "empty", True):
        return float(missing_penalty), float(missing_penalty)

    common = np.intersect1d(targ_post["current"].unique(), nest_post["current"].unique())
    if common.size == 0:
        return float(missing_penalty), float(missing_penalty)

    def extract_info(df):
        out = []
        for curr in df["current"].unique():
            spikes = np.array(df.loc[df["current"] == curr, "peak_time"].values[0])
            spikes = spikes[(spikes >= start) & (spikes <= end)]
            first_spike = spikes.min() if spikes.size > 0 else np.nan
            n_spikes = len(spikes)
            out.append({"current": curr, "first_spike": first_spike, "n_spikes": n_spikes})
        return pd.DataFrame(out)

    tm = extract_info(targ_post[targ_post["current"].isin(common)])
    nm = extract_info(nest_post[nest_post["current"].isin(common)])

    pos_tm, pos_nm = tm[tm["current"] >= 0], nm[nm["current"] >= 0]
    neg_tm, neg_nm = tm[tm["current"] < 0], nm[nm["current"] < 0]

    def spike_penalty(n_spikes, max_spikes=5):
        """Valore tra 0 e 1: cresce con n_spikes."""
        return min(1.0, n_spikes / max_spikes)

    pos_losses = []
    for c in pos_tm["current"]:
        t = pos_tm[pos_tm["current"] == c].iloc[0]
        n = pos_nm[pos_nm["current"] == c].iloc[0] if c in pos_nm["current"].values else None

        if t["n_spikes"] == 0:
            if n is None or n["n_spikes"] == 0:
                pos_losses.append(0.0)
            else:
                pos_losses.append(spike_penalty(n["n_spikes"]))
        else:
            if n is None or n["n_spikes"] == 0:
                pos_losses.append(1.0)
            else:
                err = abs(n["first_spike"] - t["first_spike"]) / max(t["first_spike"], 1e-9)
                pos_losses.append(np.clip(err / (1 + err), 0.0, 1.0))

    pos_loss = float(np.nanmean(pos_losses)) if pos_losses else float(missing_penalty)

    # --- NEGATIVE ---
    neg_losses = []
    for c in neg_tm["current"]:
        t = neg_tm[neg_tm["current"] == c].iloc[0]
        n = neg_nm[neg_nm["current"] == c].iloc[0] if c in neg_nm["current"].values else None

        if t["n_spikes"] == 0:
            if n is None or n["n_spikes"] == 0:
                neg_losses.append(0.0)
            else:
                neg_losses.append(spike_penalty(n["n_spikes"]))
        else:
            if n is None or n["n_spikes"] == 0:
                neg_losses.append(1.0)
            else:
                err = abs(n["first_spike"] - t["first_spike"]) / max(t["first_spike"], 1e-9)
                neg_losses.append(np.clip(err / (1 + err), 0.0, 1.0))

    neg_loss = float(np.nanmean(neg_losses)) if neg_losses else float(missing_penalty)

    return pos_loss, neg_loss
