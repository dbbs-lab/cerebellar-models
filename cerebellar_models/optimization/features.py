import efel
import os
import numpy as np
from typing import Optional
import re
import pandas as pd
import inspect
from scipy.stats import linregress


def time_to_first_spike(spike_train, start_stim):
    if len(spike_train) > 0:
        time_to_first_spike = spike_train[0] - start_stim
    else:
        time_to_first_spike = None
    return time_to_first_spike


def time_to_second_spike(spike_train, start_stim):
    if len(spike_train) > 1:
        time_to_second_spike = spike_train[1] - start_stim
    else:
        time_to_second_spike = None
    return time_to_second_spike


def time_to_last_spike(spike_train, start_stim):
    if len(spike_train) > 0:
        time_to_last_spike = spike_train[-1] - start_stim
    else:
        time_to_last_spike = None
    return time_to_last_spike


def inv_first_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 0:
        inv_first_ISI = 1000.0 / all_isi_values_vec[0]
    else:
        inv_first_ISI = None
    return inv_first_ISI


def inv_second_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 1:
        inv_second_ISI = 1000.0 / all_isi_values_vec[1]
    else:
        inv_second_ISI = None
    return inv_second_ISI


def inv_third_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 2:
        inv_third_ISI = 1000.0 / all_isi_values_vec[2]
    else:
        inv_third_ISI = None
    return inv_third_ISI


def inv_fourth_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 3:
        inv_fourth_ISI = 1000.0 / all_isi_values_vec[3]
    else:
        inv_fourth_ISI = None
    return inv_fourth_ISI


def inv_fifth_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 4:
        inv_fifth_ISI = 1000.0 / all_isi_values_vec[4]
    else:
        inv_fifth_ISI = None
    return inv_fifth_ISI


def inv_last_isi(spike_train):
    all_isi_values_vec = np.diff(spike_train)
    if len(all_isi_values_vec) > 0:
        inv_last_ISI = 1000.0 / all_isi_values_vec[-1]
    else:
        inv_last_ISI = None
    return inv_last_ISI


def ISI_log_slope(spike_train):
    ISI = np.diff(spike_train)
    if len(ISI) < 2:
        return np.nan

    log_x = np.log(np.arange(1, len(ISI) + 1))
    log_ISI_values = np.log(ISI)

    mask = np.isfinite(log_x) & np.isfinite(log_ISI_values)
    log_x = log_x[mask]
    log_ISI_values = log_ISI_values[mask]

    if len(log_x) < 2:
        return np.nan

    try:
        slope, _ = np.polyfit(log_x, log_ISI_values, 1)
    except np.linalg.LinAlgError:
        slope = np.nan

    return slope


def ISI_CV(spike_train):
    ISI_values = np.diff(spike_train)
    ISI_mean = np.mean(ISI_values)
    ISI_CV = np.std(ISI_values, ddof=1) / ISI_mean
    return ISI_CV


def adaptation_index(
    spike_train, stim_start, stim_end, offset=0, spike_skipf=0.1, max_spike_skip=2
):
    if spike_skipf < 0 or spike_skipf >= 1:
        raise ValueError("spike_skipf should be in [0, 1).")

    spike_time = spike_train[
        (spike_train >= stim_start - offset) & (spike_train <= stim_end - offset)
    ]

    start_idx = min(max_spike_skip, round(len(spike_time) * spike_skipf))
    spike_time = spike_time[start_idx:]

    if len(spike_time) < 3:
        return np.nan

    ISI = np.diff(spike_time)
    ISI_prev = ISI[:-1]
    ISI_next = ISI[1:]

    valid_idx = (ISI_prev + ISI_next) != 0
    if not np.any(valid_idx):
        return np.nan

    AI_values = (ISI_next[valid_idx] - ISI_prev[valid_idx]) / (
        ISI_next[valid_idx] + ISI_prev[valid_idx]
    )

    return np.mean(AI_values)


def spike_count_stimint(spike_train, stim_start, stim_end):
    peaktimes_stimint = np.where((spike_train >= stim_start) & (spike_train <= stim_end))[0]
    spike_count_stimint = len(peaktimes_stimint)
    return spike_count_stimint


def mean_frequency(spike_train, stim_start, stim_end):
    condition = (stim_start < spike_train) & (spike_train < stim_end)
    selected_spikes = spike_train[condition]
    spikecount = len(selected_spikes)

    if spikecount < 2:
        return np.nan

    last_spike_time = selected_spikes[-1]
    duration = last_spike_time - stim_start
    if duration <= 0:
        return np.nan

    mean_frequency = 1000.0 * spikecount / duration  # Hz
    return mean_frequency


FEATURE_FUNCS = {
    "time_to_first_spike": time_to_first_spike,
    "time_to_second_spike": time_to_second_spike,
    "time_to_last_spike": time_to_last_spike,
    "inv_first_ISI": inv_first_isi,
    "inv_second_ISI": inv_second_isi,
    "inv_third_ISI": inv_third_isi,
    "inv_fourth_ISI": inv_fourth_isi,
    "inv_fifth_ISI": inv_fifth_isi,
    "inv_last_ISI": inv_last_isi,
    "ISI_log_slope": ISI_log_slope,
    "ISI_CV": ISI_CV,
    "adaptation_index": adaptation_index,
    "spike_count_stimint": spike_count_stimint,
    "mean_frequency": mean_frequency,
}


def _check_features(features):
    features_ = []
    for feature in features:
        if feature not in FEATURE_FUNCS.keys():
            print("Feature {} not available.".format(feature))
        else:
            features_.append(feature)
    return features_


def _call_feature(func, spike_train, start_stim, end_stim):
    sig = inspect.signature(func)
    n_params = len(sig.parameters)
    if n_params == 1:
        return func(spike_train)
    elif n_params == 2:
        return func(spike_train, start_stim)
    elif n_params >= 3:
        return func(spike_train, start_stim, end_stim)
    return None


def multicomp_features(
    data_folder: str,
    threshold: float,
    features: Optional[list[str]] = FEATURE_FUNCS.keys(),
    start_stim: Optional[float] = None,
    end_stim: Optional[float] = None,
) -> pd.DataFrame:
    if threshold is not None:
        efel.api.set_setting("Threshold", threshold)
    else:
        print("Threshold not set.")

    features_ = _check_features(features)
    data_files = [f for f in os.listdir(data_folder) if f.endswith(".txt")]

    traces = []
    currents = np.empty(len(data_files), dtype=float)

    for idx, f in enumerate(data_files):
        filepath = os.path.join(data_folder, f)
        data = np.loadtxt(filepath)

        time = data[:, 0]
        voltage = data[:, 1]

        if not start_stim:
            start_stim = time[0]
        if not end_stim:
            end_stim = time[-1]

        trace = {
            "T": time,
            "V": voltage,
            "stim_start": [start_stim],
            "stim_end": [end_stim],
        }
        traces.append(trace)

        match = re.search(r"inj(-?\d+(?:\.\d+)?)", f)
        if match:
            currents[idx] = float(match.group(1))
        else:
            currents[idx] = np.nan

    results = efel.get_feature_values(
        traces,
        features_,
    )

    df_list = []
    for curr, res in zip(currents, results):
        row = {"current": curr}
        row.update(res)
        df_list.append(row)

    df = pd.DataFrame(df_list)
    df = df.sort_values("current").reset_index(drop=True)

    return df


def point_neuron_features(
    spike_trains: list[np.ndarray],
    currents: list[np.ndarray],
    features: Optional[list[str]] = FEATURE_FUNCS.keys(),
    start_stim: Optional[float] = None,
    end_stim: Optional[float] = None,
) -> pd.DataFrame:

    features_ = _check_features(features)
    all_results = []
    for spike_train, curr in zip(spike_trains, currents):
        r = {}
        for feat in features_:
            func = FEATURE_FUNCS[feat]
            r[feat] = _call_feature(func, spike_train, start_stim, end_stim)
        r["current"] = curr
        all_results.append(r)

    df = pd.DataFrame(all_results)
    df = df.sort_values("current").reset_index(drop=True)

    return df
