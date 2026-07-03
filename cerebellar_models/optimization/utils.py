import os
import sys
from contextlib import contextmanager

import matplotlib.pyplot as plt
import numpy as np


@contextmanager
def _suppress_output():
    devnull = open(os.devnull, "w")
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = devnull, devnull
    try:
        yield
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        devnull.close()


def plot_fi_curve(
    spikes_df, start_stim_ms, end_stim_ms, duration_ms, ax=None, label=None, color="0.2"
):
    """Plot firing rate (Hz) vs injected current (pA).

    Parameters
    ----------
    spikes_df : DataFrame with columns 'current' and 'peak_time'
    start_stim_ms : stimulus onset in ms
    end_stim_ms : stimulus offset in ms
    duration_ms : total simulation duration in ms (used for current=0)
    ax : matplotlib Axes (created if None)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))

    currents = sorted(spikes_df["current"].astype(float).unique())
    rates = []
    for curr in currents:
        mask = np.isclose(spikes_df["current"].astype(float).values, curr)
        peak_vals = spikes_df.loc[mask, "peak_time"].values
        spike_times = []
        for x in peak_vals:
            if x is None:
                continue
            if isinstance(x, (list, tuple, np.ndarray)):
                spike_times.extend(np.asarray(x, dtype=float).ravel().tolist())
            else:
                try:
                    spike_times.append(float(x))
                except Exception:
                    continue
        spike_times = np.asarray(spike_times, dtype=float)
        spike_times = spike_times[np.isfinite(spike_times)]

        if curr == 0.0:
            window_duration_s = duration_ms / 1000.0
        else:
            spike_times = spike_times[(spike_times >= start_stim_ms) & (spike_times <= end_stim_ms)]
            window_duration_s = (end_stim_ms - start_stim_ms) / 1000.0

        rates.append(spike_times.size / window_duration_s)

    print("FI curve:")
    for c, r in zip(currents, rates):
        print(f"  {c:>8.1f} pA  ->  {r:.2f} Hz")

    ax.plot(currents, rates, marker="o", linewidth=1.5, markersize=5, color=color, label=label)
    ax.set_xlabel("I [pA]", fontsize=11)
    ax.set_ylabel("Firing rate [Hz]", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=10)

    return ax
