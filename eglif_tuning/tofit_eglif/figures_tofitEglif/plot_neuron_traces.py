import os
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from cerebellar_models.optimization.features import multicomp_features


def parse_current_from_filename(fname: str) -> float:
    m = re.search(r"inj(-?\d+(?:\.\d+)?)", fname)
    return float(m.group(1)) if m else np.nan


def flatten_peak_times(peak_values):
    out = []
    for x in peak_values:
        if x is None:
            continue
        if isinstance(x, (list, tuple, np.ndarray)):
            out.extend(np.asarray(x, dtype=float).ravel().tolist())
        else:
            try:
                out.append(float(x))
            except Exception:
                continue
    out = np.asarray(out, dtype=float)
    return out[np.isfinite(out)]


def style_vm_axis(ax):
    ax.tick_params(axis="both", which="major", labelsize=10, length=4, width=1)
    ax.tick_params(axis="both", which="minor", labelsize=10, length=2, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)


def style_I_axis(axI, show_xlabel=False):
    axI.spines["top"].set_visible(False)
    axI.spines["left"].set_visible(False)
    axI.spines["right"].set_visible(True)
    axI.spines["bottom"].set_visible(True)

    axI.spines["right"].set_linewidth(1)
    axI.spines["right"].set_color("#7a1f2b")
    axI.spines["bottom"].set_linewidth(1)

    axI.yaxis.set_ticks_position("right")
    axI.yaxis.set_label_position("right")
    axI.tick_params(
        axis="y", which="major", labelsize=5, colors="#7a1f2b", length=3, width=1, pad=2
    )
    axI.set_ylabel("I [pA]", color="#7a1f2b", fontsize=6, labelpad=3)
    axI.tick_params(axis="x", which="major", labelsize=9, length=3, width=1, pad=12)
    axI.tick_params(axis="x", which="minor", length=2, width=0.8, pad=12)

    if show_xlabel:
        axI.set_xlabel("Time [ms]", fontsize=10, labelpad=8)
    else:
        axI.set_xlabel("")


if __name__ == "__main__":
    protocol = {"start_stim": 200.0, "end_stim": 700.0, "duration": 2000.0}
    cell_name = "PC"
    data_folder = f"../results_tofitEglif/{cell_name}/"
    threshold = -43.0

    features = multicomp_features(data_folder, threshold=threshold)
    spikes_df = features[["peak_time", "current"]]

    data_files = [f for f in os.listdir(data_folder) if f.endswith(".txt")]
    files_curr = [(parse_current_from_filename(f), f) for f in data_files]
    files_curr = [(c, f) for c, f in files_curr if np.isfinite(c)]
    files_curr.sort(key=lambda x: x[0])

    unique_currents = sorted({c for c, _ in files_curr})
    if not unique_currents:
        raise RuntimeError(f"No valid .txt traces found in: {data_folder}")

    traces = []
    global_vmax = -np.inf
    global_tmin, global_tmax = np.inf, -np.inf
    for curr in unique_currents:
        fname = [f for c, f in files_curr if c == curr][0]
        data = np.loadtxt(os.path.join(data_folder, fname))
        t = data[:, 0].astype(float)
        v = data[:, 1].astype(float)

        traces.append((float(curr), t, v))
        global_vmax = max(global_vmax, float(np.max(v)))
        global_tmin = min(global_tmin, float(np.min(t)))
        global_tmax = max(global_tmax, float(np.max(t)))

    print(f"Global max Vm across all traces: {global_vmax:.6g} mV")

    Imax = max(1.0, max(abs(float(c)) for c in unique_currents))
    Iylim = 1.10 * Imax
    Iticks = [-Imax, 0.0, Imax]


    major_locator = mticker.MaxNLocator(nbins=5)
    minor_locator = mticker.AutoMinorLocator(2)

    n = len(traces)
    fig_h = max(3.0, 2.65 * n)
    fig = plt.figure(figsize=(7.6, fig_h), constrained_layout=False)

    height_ratios = []
    for k in range(n):
        height_ratios += [4.0, 1.6]
        height_ratios += [0.55 if k < n - 1 else 0.15]

    gs = fig.add_gridspec(
        nrows=len(height_ratios),
        ncols=1,
        height_ratios=height_ratios,
        hspace=0.10,
    )

    row = 0
    for i, (curr, t, v) in enumerate(traces):
        axV = fig.add_subplot(gs[row, 0])
        axI = fig.add_subplot(gs[row + 1, 0], sharex=axV)
        row += 3

        axV.set_facecolor("white")
        axI.set_facecolor("white")

        if curr != 0:
            axV.axvspan(
                protocol["start_stim"],
                protocol["end_stim"],
                facecolor="#ffe680",
                alpha=0.18,
                edgecolor="none",
                zorder=0,
            )

        # Vm trace + threshold
        axV.plot(t, v, color="0.35", linewidth=1.1, zorder=2)
        axV.axhline(y=threshold, color="0.6", linestyle="--", linewidth=0.9, zorder=1)

        axV.set_ylabel(r"$V_m$ [mV]", fontsize=11)
        axV.set_ylim(-100, 25)
        axV.set_xlim(global_tmin, global_tmax)
        style_vm_axis(axV)
        axV.tick_params(axis="x", labelbottom=False)

        # Current step below (0 outside stim)
        I = np.zeros_like(t, dtype=float)
        I[(t >= protocol["start_stim"]) & (t <= protocol["end_stim"])] = float(curr)

        axI.plot(t, I, drawstyle="steps-post", linewidth=1.6, color="#7a1f2b")
        axI.set_ylim(-Iylim, Iylim)
        axI.set_yticks(Iticks)

        axI.xaxis.set_major_locator(major_locator)
        axI.xaxis.set_minor_locator(minor_locator)

        show_xlabel = i == n - 1
        style_I_axis(axI, show_xlabel=show_xlabel)

        for lab in axI.get_xticklabels(which="both"):
            lab.set_clip_on(False)

        # Current label in the I panel (clipped to axI)
        txt = axI.text(
            0.98,
            0.30 if curr >= 0 else 0.70,
            f"I = {curr:g} pA",
            transform=axI.transAxes,
            ha="right",
            va="center",
            fontsize=8,
            color="#7a1f2b",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.80),
            zorder=10,
        )
        txt.set_clip_on(True)
        txt.set_clip_path(axI.patch)

        # Spike markers (optional)
        mask = np.isclose(spikes_df["current"].astype(float).values, float(curr))
        peak_vals = spikes_df.loc[mask, "peak_time"].values
        spike_times = flatten_peak_times(peak_vals)
        if spike_times.size > 0:
            pass
            # y0, y1 = global_vmax + 2.0, global_vmax + 6.0
            # axV.vlines(spike_times, y0, y1, color="darkred", linewidth=0.9)

    fig.subplots_adjust(top=0.985, bottom=0.085, left=0.12, right=0.93)
    fig.savefig(f"{cell_name}_exp.png", dpi=700, bbox_inches="tight")
    fig.savefig(f"{cell_name}_exp.pdf", dpi=700, bbox_inches="tight")
    plt.show()
