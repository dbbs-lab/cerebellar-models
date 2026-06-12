import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from bsb import from_storage
from cerebellar_models.analysis import StructureReport
from cerebellar_models.MFM.transfer_functions import TransferFunctionSimulator

SCAFFOLD_PATH = 'mouse_cerebellum.hdf5'
OUTPUT_DIR    = 'output_TF'


def load_tf(cell_name):
    path = os.path.join(OUTPUT_DIR, f'TF_{cell_name}.npz')
    data = np.load(path, allow_pickle=True)
    return data['tf_mean'], data['tf_std'], data['tags'].tolist(), [
        data[f'freq_axis_{tag}'] for tag in data['tags']
    ]


if __name__ == '__main__':
    scaffold = from_storage(SCAFFOLD_PATH)
    report = StructureReport(scaffold)
    report.plots['density_table'].update()
    counts = report.plots['density_table'].get_counts()

    n_bc = counts['BC']
    n_sc = counts['SC']
    total = n_bc + n_sc
    w_bc = n_bc / total
    w_sc = n_sc / total
    print(f"basket_cell: {n_bc}  stellate_cell: {n_sc}  (w_BC={w_bc:.3f}, w_SC={w_sc:.3f})")

    mean_bc, std_bc, tags_bc, axes_bc = load_tf('basket_cell')
    mean_sc, std_sc, tags_sc, axes_sc = load_tf('stellate_cell')

    if tags_bc != tags_sc:
        raise ValueError(f"TF axes mismatch: {tags_bc} vs {tags_sc}")
    for ax_b, ax_s, tag in zip(axes_bc, axes_sc, tags_bc):
        if not np.allclose(ax_b, ax_s):
            raise ValueError(f"freq_axis mismatch for label '{tag}'")

    tf_mean = w_bc * mean_bc + w_sc * mean_sc
    tf_std = np.sqrt(w_bc**2 * std_bc**2 + w_sc**2 * std_sc**2)

    tags = tags_bc
    freq_axes = axes_bc

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_path = os.path.join(OUTPUT_DIR, 'TF_MLI.npz')
    np.savez(
        out_path,
        tf_mean=tf_mean,
        tf_std=tf_std,
        tags=np.array(tags),
        **{f'freq_axis_{t}': ax for t, ax in zip(tags, freq_axes)},
    )
    print(f"Saved → {out_path}")

    colors = TransferFunctionSimulator._CELL_COLORS
    mli_color = [
        (a + b) / 2
        for a, b in zip(colors['basket_cell'], colors['stellate_cell'])
    ]
    cmap = LinearSegmentedColormap.from_list('MLI', [(1, 1, 1, 1), mli_color])

    extent = [freq_axes[0][0], freq_axes[0][-1], freq_axes[1][0], freq_axes[1][-1]]

    fig1, ax1 = plt.subplots(figsize=(5, 4))
    im1 = ax1.imshow(tf_mean.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
    ax1.set_title("MLI — Mean firing rate (pop. weighted)")
    ax1.set_xlabel(tags[0])
    ax1.set_ylabel(tags[1])
    fig1.colorbar(im1, ax=ax1, label='Hz')
    fig1.tight_layout()
    fig1.savefig('TF_MLI.png', dpi=150)

    fig2, ax2 = plt.subplots(figsize=(5, 4))
    im2 = ax2.imshow(tf_std.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
    ax2.set_title("MLI — STD (pop. weighted)")
    ax2.set_xlabel(tags[0])
    ax2.set_ylabel(tags[1])
    fig2.colorbar(im2, ax=ax2, label='Hz')
    fig2.tight_layout()
    fig2.savefig('TF_std_MLI.png', dpi=150)

    plt.show()