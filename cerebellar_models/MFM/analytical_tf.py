"""
Semi-analytical MFM transfer functions module.
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
import yaml
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.optimize import minimize


# order expected by compute_statistics_eglif
_PSEQ_KEYS = ("Qe", "Te", "Ee", "Qi", "Ti", "Ei", "Gl", "Cm", "El", "Ke", "Ki")


def pseq_params_eglif(p):
    return tuple(p[k] for k in _PSEQ_KEYS)


def erfc_func(muV, sV, TvN, Vthre, Gl, Cm, alpha):
    """TF expression: nu_out."""
    return 0.5 / TvN * Gl / Cm * special.erfc((Vthre - muV) / np.sqrt(2) / sV) * alpha


def effective_Vthre(Y, muV, sV, TvN, Gl, Cm, alpha):
    """Inverse of erfc_func: threshold that reproduces the numerical rate Y."""
    return muV + np.sqrt(2) * sV * special.erfcinv((1 / alpha) * (Y * 2 * TvN * Cm / Gl))


def threshold_func(muV, sV, TvN, muGn, *P):
    """Phenomenological threshold in the fluctuation moments.

    1st order: linear in (muV, sV, TvN) plus a log(muGn) term (EITN notebook
    form), i.e. ``len(P) == 5``.
    """
    muV0, DmuV0 = -60e-3, 10e-3
    sV0, DsV0 = 4e-3, 6e-3
    TvN0, DTvN0 = 0.5, 1.0
    dmu = (muV - muV0) / DmuV0
    dsig = (sV - sV0) / DsV0
    dtau = (TvN - TvN0) / DTvN0
    P0, P1, P2, P3, P4 = P
    return P0 + P1 * dmu + P2 * dsig + P3 * dtau + P4 * np.log(muGn)


def compute_statistics_eglif(Fe, Fi, XX, Qe, Te, Ee, Qi, Ti, Ei, Gl, Cm, El, Ke, Ki):
    """
    Fluctuation moments (muV, sigmaV, tauV) from shot-noise theory.
    XX is the adaptation variable (disabled, XX = 0).
    """
    fe, fi = Fe, Fi
    muGe, muGi = Qe * Ke * Te * fe, Qi * Ki * Ti * fi  # population conductances
    muG = Gl + muGe + muGi  # total input conductance
    muV = (np.e * (muGe * Ee + muGi * Ei + Gl * El) - XX) / muG
    muGn, Tm = muG / Gl, Cm / muG
    Ue, Ui = Qe / muG * (Ee - muV), Qi / muG * (Ei - muV)
    sVe = (2 * Tm + Te) * ((np.e * Ue * Te) / (2 * (Te + Tm))) ** 2 * Ke * fe
    sVi = (2 * Tm + Ti) * ((np.e * Ui * Ti) / (2 * (Ti + Tm))) ** 2 * Ki * fi
    sV = np.sqrt(sVe + sVi)
    fe, fi = fe + 1e-9, fi + 1e-9  # avoid division by zero
    Tv_num = Ke * fe * Ue**2 * Te**2 * np.e**2 + Ki * fi * Ui**2 * Ti**2 * np.e**2
    Tv = 0.5 * Tv_num / ((sV + 1e-20) ** 2)
    TvN = Tv * Gl / Cm
    return muGe, muGi, muG, muV, sV + 1e-20, muGn, TvN, Tv


def analytical_tf_rate(fe, fi, XX, alpha, params, P):
    """Analytical TF output rate for excitatory rate fe and inhibitory rate fi."""
    _, _, _, muV, sV, muGn, TvN, _ = compute_statistics_eglif(
        fe, fi, XX, *pseq_params_eglif(params)
    )
    return erfc_func(
        muV, sV, TvN, threshold_func(muV, sV, TvN, muGn, *P), params["Gl"], params["Cm"], alpha
    )


def fit_two_step(Fout, Fe_eff, fiSim, w, params, alpha, maxiter=50000, verbose=False):
    """Step 1: fit in threshold space; step 2: fit in rate space. Returns (P, MSE)."""
    _, _, _, muV, sV, muGn, TvN, _ = compute_statistics_eglif(
        Fe_eff, fiSim, w, *pseq_params_eglif(params)
    )

    i = np.where((Fout > 0.0) & (Fout < Fout.max()))
    Vth = effective_Vthre(Fout[i], muV[i], sV[i], TvN[i], params["Gl"], params["Cm"], alpha)
    ok = np.isfinite(Vth)
    if ok.sum() < 6:  # alpha too small: threshold not defined
        return None, np.inf

    res1 = lambda p: np.mean(
        (Vth[ok] - threshold_func(muV[i][ok], sV[i][ok], TvN[i][ok], muGn[i][ok], *p)) ** 2
    )
    p0 = np.zeros(5)
    p0[0] = -50e-3
    P = minimize(res1, p0).x

    res2 = lambda p: np.mean((Fout - analytical_tf_rate(Fe_eff, fiSim, w, alpha, params, p)) ** 2)
    r = minimize(
        res2,
        P,
        method="nelder-mead",
        options={"xatol": 1e-7, "fatol": 1e-6, "maxiter": maxiter, "maxfev": maxiter},
    )
    if verbose:
        print(f"  step 1 P [mV]: {np.round(P * 1e3, 3)}")
        print(f"  step 2 P [mV]: {np.round(r.x * 1e3, 3)}   MSE = {r.fun:.3f} Hz^2")
    return r.x, r.fun


def fit_alpha_and_P(Fout, Fe_eff, fiSim, w, params, alpha0, P0):
    """Joint refinement of (alpha, P) in rate space, from a starting point."""

    def obj(x):
        if x[0] <= 0:
            return np.inf
        return np.mean((Fout - analytical_tf_rate(Fe_eff, fiSim, w, x[0], params, x[1:])) ** 2)

    r = minimize(
        obj,
        np.r_[alpha0, P0],
        method="nelder-mead",
        options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 100000, "maxfev": 100000},
    )
    return r.x[0], r.x[1:], r.fun


def nest_to_mf(constants, syn, exc_receptor=1, inh_receptor=2, name="cell"):
    """
    Convert NEST-unit E-GLIF ``constants`` (ms, mV, pF, pA, nS/ms) to the SI
    dictionary used by the MF formulas. ``syn`` carries the network quantities
    (Qe, Qi in S; Ke, Ki in-degrees); ``exc_receptor`` / ``inh_receptor`` select
    which ``tau_syn{r}`` / ``E_rev{r}`` are the excitatory / inhibitory inputs.
    """
    ms, mV, pF, pA = 1e-3, 1e-3, 1e-12, 1e-12
    out = {
        "name": name,
        "N": 1,
        "Cm": constants["C_m"] * pF,
        "tau_m": constants["tau_m"] * ms,
        "El": constants["E_L"] * mV,
        "V_reset": constants["V_reset"] * mV,
        "V_th": constants["V_th"] * mV,
        "t_ref": constants["t_ref"] * ms,
        "A1": constants["A1"] * pA,
        "A2": constants["A2"] * pA,
        "I_e": constants["I_e"] * pA,
        "k_1": constants["k_1"] / ms,  # 1/ms -> 1/s
        "k_2": constants["k_2"] / ms,
        "k_adap": constants["k_adap"] * 1e-9 / ms,  # nS/ms -> S/s
        "Te": constants[f"tau_syn{exc_receptor}"] * ms,
        "Ee": constants[f"E_rev{exc_receptor}"] * mV,
        "Ti": constants[f"tau_syn{inh_receptor}"] * ms,
        "Ei": constants[f"E_rev{inh_receptor}"] * mV,
    }
    out["Gl"] = out["Cm"] / out["tau_m"]  # S
    out.update(syn)
    return out


_CELL_COLORS = {
    "granule_cell": [0.7, 0.15, 0.15, 1.0],
    "golgi_cell": [0, 0.45, 0.7, 1.0],
    "purkinje_cell": [0.275, 0.800, 0.275, 1.0],
    "purkinje_cell_Z+": [0.275, 0.800, 0.275, 1.0],
    "purkinje_cell_Z-": [0.275, 0.800, 0.275, 1.0],
    "basket_cell": [1, 0.647, 0, 1.0],
    "stellate_cell": [1, 0.84, 0, 1.0],
    "MLI": [1, 0.7435, 0, 1.0],
}


def _cell_cmap(cell_name):
    color = _CELL_COLORS.get(cell_name, [0.275, 0.800, 0.275, 1.0])
    return LinearSegmentedColormap.from_list(cell_name, [(1, 1, 1, 1), color])




class AnalyticalTransferFunction:
    """
    Holds a fitted (or loaded) analytical TF: the gain ``alpha``, the threshold
    coefficients ``P`` and the SI parameter set needed to evaluate it. Two
    inputs only: an excitatory and an inhibitory presynaptic rate.
    """

    def __init__(self, cell_name, params, alpha, P, tags, freq_axes=None, mse=None):
        self.cell_name = cell_name
        self.params = params
        self.alpha = float(alpha)
        self.P = np.asarray(P, dtype=float)
        self.tags = list(tags)  # [exc_label, inh_label]
        self.freq_axes = None if freq_axes is None else [np.asarray(a) for a in freq_axes]
        self.mse = None if mse is None else float(mse)

    def evaluate(self, fe, fi, w=0.0):
        """Analytical output rate at excitatory rate ``fe`` and inhibitory rate ``fi``."""
        return analytical_tf_rate(fe, fi, w, self.alpha, self.params, self.P)

    __call__ = evaluate

    def evaluate_on_grid(self, freq_axes=None):
        """Evaluate on the (exc, inh) grid; returns an [inh, exc] array."""
        axes = freq_axes if freq_axes is not None else self.freq_axes
        if axes is None:
            raise RuntimeError("No freq_axes available; pass one explicitly.")
        fe_ax, fi_ax = axes
        FI, FE = np.meshgrid(np.asarray(fi_ax), np.asarray(fe_ax), indexing="ij")
        return self.evaluate(FE, FI)

    @classmethod
    def load(cls, cell_name, output_dir="output_analytical_TFs"):
        """Load a saved analytical TF from ``output_dir/analyticalTF_{cell_name}.json``."""
        path = os.path.join(output_dir, f"analyticalTF_{cell_name}.json")
        with open(path) as f:
            d = json.load(f)
        return cls(
            cell_name,
            d["params"],
            d["alpha"],
            d["P"],
            d["tags"],
            freq_axes=d.get("freq_axes"),
            mse=d.get("mse"),
        )

    def save(self, output_dir="output_analytical_TFs"):
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f"analyticalTF_{self.cell_name}.json")
        params = {
            k: (v.item() if isinstance(v, np.generic) else v) for k, v in self.params.items()
        }
        d = {
            "cell_name": self.cell_name,
            "synaptic_inputs": {
                self.tags[0]: {"role": "exc", "Q_nS": params["Qe"] * 1e9, "K": params["Ke"]},
                self.tags[1]: {"role": "inh", "Q_nS": params["Qi"] * 1e9, "K": params["Ki"]},
            },
            "params": params,
            "alpha": self.alpha,
            "P": self.P.tolist(),
            "tags": self.tags,
            "freq_axes": None if self.freq_axes is None else [a.tolist() for a in self.freq_axes],
            "mse": self.mse,
        }
        with open(path, "w") as f:
            json.dump(d, f, indent=2)
        print(f"Analytical TF saved → {path}")
        return path

    def metrics(self, numerical_mean):
        """RMSE / MAE / max|err| / R² of the analytical TF against a numerical grid."""
        an = self.evaluate_on_grid()
        r = np.asarray(numerical_mean) - an
        return {
            "RMSE": float(np.sqrt((r**2).mean())),
            "MAE": float(np.abs(r).mean()),
            "maxAE": float(np.abs(r).max()),
            "R2": float(1 - (r**2).sum() / ((numerical_mean - numerical_mean.mean()) ** 2).sum()),
        }

    def summary(self):
        print(f"Cell    : {self.cell_name.replace('_', ' ')}")
        print(f"Inputs  : exc={self.tags[0]}, inh={self.tags[1]}")
        print(f"alpha   : {self.alpha:.3f}")
        print("P [mV]  : " + ", ".join(f"P{i}={p * 1e3:.2f}" for i, p in enumerate(self.P)))
        if self.freq_axes is not None:
            fe, fi = self.freq_axes
            print(
                f"Domain  : {self.tags[0]} [{fe[0]:.1f}–{fe[-1]:.1f} Hz], "
                f"{self.tags[1]} [{fi[0]:.1f}–{fi[-1]:.1f} Hz]"
            )
        if self.mse is not None:
            print(f"Fit MSE : {self.mse:.3f} Hz^2  (RMSE {np.sqrt(self.mse):.2f} Hz)")

    def plot_tf(self, numerical_tf=None, save=True, n=60, output_dir="."):
        """3D surface of the analytical TF, optionally with the numerical dots."""
        if self.freq_axes is None:
            raise RuntimeError("No freq_axes stored; cannot pick a plotting domain.")
        fe_ax, fi_ax = self.freq_axes
        FE, FI = np.meshgrid(
            np.linspace(fe_ax[0], fe_ax[-1], n), np.linspace(fi_ax[0], fi_ax[-1], n)
        )
        Z = self.evaluate(FE, FI)

        num_mean = num_FE = num_FI = None
        zmax = Z.max()
        if numerical_tf is not None:
            ie = numerical_tf.tags.index(self.tags[0])
            ii = numerical_tf.tags.index(self.tags[1])
            exc_axis = numerical_tf.freq_axes[ie]
            inh_axis = numerical_tf.freq_axes[ii]
            num_mean = np.moveaxis(numerical_tf.tf_mean, (ie, ii), (0, 1)).T  # [inh, exc]
            num_FE, num_FI = np.meshgrid(exc_axis, inh_axis)
            zmax = max(zmax, num_mean.max())

        cmap = _cell_cmap(self.cell_name)
        norm = Normalize(0, zmax)
        cell_label = self.cell_name.replace("_", " ")

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(
            FE, FI, Z, facecolors=cmap(norm(Z)), rstride=1, cstride=1, lw=0, alpha=0.85, shade=False
        )
        ax.plot_wireframe(FE, FI, Z, rstride=6, cstride=6, color="k", lw=0.3, alpha=0.35)
        if num_mean is not None:
            ax.scatter(num_FE, num_FI, num_mean, c="k", s=9, depthshade=False)
        ax.set_xlabel(f"{self.tags[0]} [Hz]")
        ax.set_ylabel(f"{self.tags[1]} [Hz]")
        ax.set_zlabel(f"{cell_label} [Hz]")
        ax.set_zlim(0, zmax)
        ax.view_init(elev=18, azim=-150)
        title = f"{cell_label}   α = {self.alpha:.2f}"
        if self.mse is not None:
            title += f"   RMSE = {np.sqrt(self.mse):.1f} Hz"
        ax.set_title(title)
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5, pad=0.1)
        cb.set_label(f"{cell_label} [Hz]")
        fig.tight_layout()

        if save:
            os.makedirs(output_dir, exist_ok=True)
            fig.savefig(os.path.join(output_dir, f"analyticalTF_{self.cell_name}.png"), dpi=150)
        return fig


# --------------------------------------------------------------------------- #
# Fitting driver                                                              #
# --------------------------------------------------------------------------- #


class AnalyticalTransferFunctionFitting:
    """
    Fits the analytical TF on a numerical template. Single-cell and network
    parameters come from the BSB scaffold + circuit YAML (same sources as
    ``NumericalTransferFunctionSimulator``) and are converted to SI units.

    The two TF inputs (one excitatory, one inhibitory) are inferred from the
    connections in ``tf_dict``: connections are grouped by label, each label's
    receptor fixes its reversal potential, and the label with the higher
    reversal is taken as excitatory. Multiple connections sharing a label are
    aggregated into one (Q, K): Q is the arithmetic mean of the per-connection
    quantal conductances (each weight scaled by its synapses-per-pair, since
    synaptic time constants are shared within a receptor); K sums the in-degree
    within a presynaptic population and averages across the distinct
    populations fused into the label (which collapse to one mean cell in the
    MFM), i.e. K = sum(convergences) / n_presynaptic_populations.
    """

    def __init__(
        self, scaffold_path, yaml_path, tf_dict_path, cell_name, numerical_tf, cut=2,
    ):
        self.cell_name = cell_name
        self.cut = cut
        self.numerical_tf = numerical_tf

        self._convergences, self._synapses_per_pair = self._scaffold_convergences(scaffold_path)
        self.params, self.exc_label, self.inh_label = self._build_params(
            yaml_path, tf_dict_path, cell_name
        )
        self._prepare_template()
        self._fit = None

    # -- SNN info ----------------------------------------------------------- #

    @staticmethod
    def _scaffold_convergences(scaffold_path):
        from bsb import from_storage

        from cerebellar_models.analysis import StructureReport

        scaffold = from_storage(scaffold_path)
        report = StructureReport(scaffold)
        conn_table = report.plots["connectivity_table"]
        conn_table.update()
        conv = {
            ps.tag: float(np.mean(arr))
            for ps, arr in zip(
                scaffold.get_connectivity_sets(), conn_table.get_convergences().values()
            )
        }
        spp = {
            ps.tag: float(np.mean(arr))
            for ps, arr in zip(
                scaffold.get_connectivity_sets(), conn_table.get_nb_synapse_per_pair().values()
            )
        }
        return conv, spp

    @staticmethod
    def _load_yaml(yaml_path):
        with open(yaml_path) as f:
            data = yaml.safe_load(f)
        if "$import" not in data:
            return data
        import_info = data.pop("$import")
        ref = import_info["ref"].split("#")[0]
        base_path = (
            ref
            if os.path.isabs(ref)
            else os.path.join(os.path.dirname(os.path.abspath(yaml_path)), ref)
        )
        base = AnalyticalTransferFunctionFitting._load_yaml(base_path)
        return AnalyticalTransferFunctionFitting._deep_merge(base, data)

    @staticmethod
    def _deep_merge(base, override):
        result = dict(base)
        for k, v in override.items():
            if k in result and isinstance(result[k], dict) and isinstance(v, dict):
                result[k] = AnalyticalTransferFunctionFitting._deep_merge(result[k], v)
            else:
                result[k] = v
        return result

    def _build_params(self, yaml_path, tf_dict_path, cell_name):
        circuit = self._load_yaml(yaml_path)
        with open(tf_dict_path) as f:
            tf_dict = json.load(f)

        nest_cfg = circuit["simulations"]["nest_basal_activity"]
        constants = nest_cfg["cell_models"][cell_name]["constants"]
        conn_models = nest_cfg["connection_models"]
        connectivity = circuit.get("connectivity", {})
        cell_conns = tf_dict[cell_name]["connections"]

        # group connections by label. K sums the in-degree within each
        # presynaptic population but averages across the distinct populations
        # fused into the label (in the MFM they collapse to one mean cell),
        # i.e. K = sum(convergences) / n_presynaptic_populations. Q is the
        # arithmetic mean of theanc per-connection quantal conductances.
        groups = {}
        for tag, info in cell_conns.items():
            label = info["label"]
            syn = conn_models[tag]["synapses"][0]
            receptor = int(syn.get("receptor_type", 1))
            K = self._convergences[tag]
            spp = self._synapses_per_pair.get(tag, 1.0)
            Q = float(syn["weight"]) * spp * 1e-9  # nS -> S
            pre = connectivity.get(tag, {}).get("presynaptic", {}).get("cell_types", [tag])
            g = groups.setdefault(label, {"receptor": receptor, "K": 0.0, "Qs": [], "pops": set()})
            if g["receptor"] != receptor:
                raise ValueError(f"label '{label}' mixes receptor types")
            g["K"] += K
            g["Qs"].append(Q)
            g["pops"].update(pre)

        if len(groups) != 2:
            raise ValueError(
                f"analytical TF needs exactly 2 input labels (exc, inh); "
                f"{cell_name} has {sorted(groups)}"
            )

        for g in groups.values():
            g["E_rev"] = constants[f"E_rev{g['receptor']}"]
            g["Q"] = float(np.mean(g["Qs"])) if g["Qs"] else 0.0
            g["K"] = g["K"] / max(len(g["pops"]), 1)
        exc_label, inh_label = sorted(groups, key=lambda l: groups[l]["E_rev"], reverse=True)

        syn = {
            "Qe": groups[exc_label]["Q"],
            "Ke": groups[exc_label]["K"],
            "Qi": groups[inh_label]["Q"],
            "Ki": groups[inh_label]["K"],
        }
        params = nest_to_mf(
            constants,
            syn,
            exc_receptor=groups[exc_label]["receptor"],
            inh_receptor=groups[inh_label]["receptor"],
            name=cell_name,
        )
        return params, exc_label, inh_label


    def _prepare_template(self):
        ntf = self.numerical_tf
        ie = ntf.tags.index(self.exc_label)
        ii = ntf.tags.index(self.inh_label)
        c = self.cut
        self.Fe = ntf.freq_axes[ie][c:]
        self.Fi = ntf.freq_axes[ii][c:]
        # orient to [exc, inh], drop the cut bins, store as [inh, exc]
        mean = np.moveaxis(ntf.tf_mean, (ie, ii), (0, 1))[c:, c:]
        std = np.moveaxis(ntf.tf_std, (ie, ii), (0, 1))[c:, c:]
        self.MEAN = mean.T
        self.SD = std.T
        self.fiSim, self.Fe_eff = np.meshgrid(self.Fi, self.Fe, indexing="ij")
        self.w = np.zeros_like(self.MEAN)

    def _make_result(self, alpha, P, mse):
        return AnalyticalTransferFunction(
            self.numerical_tf.cell_name,
            self.params,
            alpha,
            P,
            [self.exc_label, self.inh_label],
            freq_axes=[self.Fe, self.Fi],
            mse=mse,
        )

    # -- fitting ------------------------------------------------------------ #

    def fit(self, alpha=4.5, save=True, output_dir="output_analytical_TFs", verbose=True):
        """Two-step fit at a fixed gain ``alpha``. Returns an AnalyticalTransferFunction."""
        P, mse = fit_two_step(
            self.MEAN, self.Fe_eff, self.fiSim, self.w, self.params, alpha,
            verbose=verbose,
        )
        self._fit = self._make_result(alpha, P, mse)
        if save:
            self._fit.save(output_dir)
        return self._fit

    def _scan_alpha(self, alphas):
        out = [
            fit_two_step(
                self.MEAN, self.Fe_eff, self.fiSim, self.w, self.params, a,
            )
            for a in alphas
        ]
        return np.array([o[1] for o in out]), [o[0] for o in out]

    def optimize_alpha(
        self,
        alpha_coarse=None,
        alpha_fine=None,
        save=True,
        output_dir="output_analytical_TFs",
        verbose=True,
    ):
        """
        Coarse + fine scan of alpha followed by a joint (alpha, P) refinement.
        Returns the optimised AnalyticalTransferFunction.
        """
        if alpha_coarse is None:
            alpha_coarse = np.round(np.arange(1.0, 15.001, 0.25), 2)
        if alpha_fine is None:
            alpha_fine = np.round(np.arange(2.0, 6.001, 0.05), 2)

        self.scan_coarse = (alpha_coarse, self._scan_alpha(alpha_coarse)[0])
        mse_fine, Ps = self._scan_alpha(alpha_fine)
        self.scan_fine = (alpha_fine, mse_fine)

        b = int(np.nanargmin(mse_fine))
        a_opt, P_opt, mse_opt = fit_alpha_and_P(
            self.MEAN, self.Fe_eff, self.fiSim, self.w, self.params, alpha_fine[b], Ps[b]
        )
        self._fit = self._make_result(a_opt, P_opt, mse_opt)
        if verbose:
            print(
                f"{self.numerical_tf.cell_name}: best fine-scan α = {alpha_fine[b]:.2f} "
                f"→ joint α* = {a_opt:.3f}   RMSE = {np.sqrt(mse_opt):.2f} Hz"
            )
        if save:
            self._fit.save(output_dir)
        return self._fit

    # -- plotting ----------------------------------------------------------- #

    def plot_tf(self, save=True, output_dir="."):
        """3D surface of the current fit on top of the numerical dots."""
        if self._fit is None:
            raise RuntimeError("Run fit() or optimize_alpha() first.")
        return self._fit.plot_tf(numerical_tf=self.numerical_tf, save=save, output_dir=output_dir)

    def plot_fit(self, save=True, output_dir="."):
        """Sections at fixed inhibition: numerical (diamonds, ±SD) vs analytical (lines)."""
        if self._fit is None:
            raise RuntimeError("Run fit() or optimize_alpha() first.")
        fe = np.linspace(self.Fe[0], self.Fe[-1], 200)
        cmap = plt.cm.viridis
        n = len(self.Fi)

        fig, ax = plt.subplots(figsize=(6, 4.5))
        for j, fi in enumerate(self.Fi):
            col = cmap(j / max(n - 1, 1))
            ax.errorbar(
                self.Fe, self.MEAN[j], yerr=self.SD[j], color=col, marker="D", ms=4,
                capsize=2, elinewidth=1, lw=0,
            )
            ax.plot(fe, self._fit.evaluate(fe, fi), color=col)
        m = self._fit.metrics(self.MEAN)
        ax.set_title(
            f"{self.numerical_tf.cell_name}   α = {self._fit.alpha:.2f}   RMSE = {m['RMSE']:.1f} Hz"
        )
        ax.set_xlabel(f"{self.exc_label} [Hz]")
        ax.set_ylabel(f"{self.numerical_tf.cell_name.replace('_', ' ')} [Hz]")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(self.Fi[0], self.Fi[-1]))
        fig.colorbar(sm, ax=ax, label=f"{self.inh_label} [Hz]")
        fig.tight_layout()

        if save:
            os.makedirs(output_dir, exist_ok=True)
            fig.savefig(
                os.path.join(output_dir, f"analyticalTF_fit_{self.numerical_tf.cell_name}.png"),
                dpi=150,
            )
        return fig