import json
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
import yaml
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.optimize import minimize


def erfc_func(muV, sV, TvN, Vthre, Gl, Cm, alpha):
    """TF expression: nu_out."""
    return 0.5 / TvN * Gl / Cm * special.erfc((Vthre - muV) / np.sqrt(2) / sV) * alpha


def effective_Vthre(Y, muV, sV, TvN, Gl, Cm, alpha):
    """Inverse of erfc_func: threshold that reproduces the numerical rate Y."""
    return muV + np.sqrt(2) * sV * special.erfcinv((1 / alpha) * (Y * 2 * TvN * Cm / Gl))


def threshold_func(muV, sV, TvN, muGn, *P):
    """
    Phenomenological threshold in the fluctuation moments (Zerlaut/EITN form).
    The number of coefficients sets the order: 5 -> first order (muV, sV, TvN,
    log(muGn)); up to 11 -> second order with the quadratic and cross terms.
    """
    muV0, DmuV0 = -60e-3, 10e-3
    sV0, DsV0 = 4e-3, 6e-3
    TvN0, DTvN0 = 0.5, 1.0
    v = (muV - muV0) / DmuV0
    s = (sV - sV0) / DsV0
    t = (TvN - TvN0) / DTvN0
    # canonical basis: [1, v, s, t, log(muGn), v², s², t², v·s, v·t, s·t]
    basis = [1.0, v, s, t, np.log(muGn), v * v, s * s, t * t, v * s, v * t, s * t]
    if len(P) > len(basis):
        raise ValueError(f"threshold_func takes at most {len(basis)} coefficients, got {len(P)}")
    return sum(P[k] * basis[k] for k in range(len(P)))


def compute_statistics_eglif(fe, fi, XX, params):
    """
    Fluctuation moments (muV, sigmaV, tauV) from shot-noise theory for a TF with
    ``len(params["exc"])`` excitatory and ``len(params["inh"])`` inhibitory
    channels. ``fe`` / ``fi`` are lists of rates aligned with the respective
    channel lists (one entry per channel); ``XX`` is the adaptation variable.
    """
    exc, inh = params["exc"], params["inh"]
    Gl, Cm, El = params["Gl"], params["Cm"], params["El"]
    if len(fe) != len(exc):
        raise ValueError(f"expected {len(exc)} excitatory rates, got {len(fe)}")
    if len(fi) != len(inh):
        raise ValueError(f"expected {len(inh)} inhibitory rates, got {len(fi)}")

    muGe = [c["Q"] * c["K"] * c["Te"] * f for c, f in zip(exc, fe)]
    muGi = [c["Q"] * c["K"] * c["Ti"] * f for c, f in zip(inh, fi)]
    muG = Gl + sum(muGe) + sum(muGi)  # total input conductance
    muV = (
        np.e
        * (
            sum(g * c["Ee"] for g, c in zip(muGe, exc))
            + sum(g * c["Ei"] for g, c in zip(muGi, inh))
            + Gl * El
        )
        - XX
    ) / muG
    muGn, Tm = muG / Gl, Cm / muG

    sV2, Tv_num = 0.0, 0.0
    for c, f in zip(exc, fe):
        U, Te = c["Q"] / muG * (c["Ee"] - muV), c["Te"]
        sV2 = sV2 + (2 * Tm + Te) * ((np.e * U * Te) / (2 * (Te + Tm))) ** 2 * c["K"] * f
        Tv_num = Tv_num + c["K"] * (f + 1e-9) * U**2 * Te**2 * np.e**2
    for c, f in zip(inh, fi):
        U, Ti = c["Q"] / muG * (c["Ei"] - muV), c["Ti"]
        sV2 = sV2 + (2 * Tm + Ti) * ((np.e * U * Ti) / (2 * (Ti + Tm))) ** 2 * c["K"] * f
        Tv_num = Tv_num + c["K"] * (f + 1e-9) * U**2 * Ti**2 * np.e**2
    sV = np.sqrt(sV2)
    Tv = 0.5 * Tv_num / ((sV + 1e-20) ** 2)
    TvN = Tv * Gl / Cm
    return muV, sV + 1e-20, muGn, TvN


def analytical_tf_rate(fe, fi, XX, alpha, params, P):
    """Analytical TF output rate for excitatory rate(s) ``fe`` and inhibitory rate(s) ``fi``."""
    muV, sV, muGn, TvN = compute_statistics_eglif(fe, fi, XX, params)
    return erfc_func(
        muV, sV, TvN, threshold_func(muV, sV, TvN, muGn, *P), params["Gl"], params["Cm"], alpha
    )


def fit_two_step(Fout, fe, fi, w, params, alpha, maxiter=50000, verbose=False):
    """Step 1: fit in threshold space; step 2: fit in rate space. Returns (P, MSE)."""
    muV, sV, muGn, TvN = compute_statistics_eglif(fe, fi, w, params)

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

    res2 = lambda p: np.mean((Fout - analytical_tf_rate(fe, fi, w, alpha, params, p)) ** 2)
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


def fit_alpha_and_P(Fout, fe, fi, w, params, alpha0, P0):
    """Joint refinement of (alpha, P) in rate space, from a starting point."""

    def obj(x):
        if x[0] <= 0:
            return np.inf
        return np.mean((Fout - analytical_tf_rate(fe, fi, w, x[0], params, x[1:])) ** 2)

    r = minimize(
        obj,
        np.r_[alpha0, P0],
        method="nelder-mead",
        options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 100000, "maxfev": 100000},
    )
    return r.x[0], r.x[1:], r.fun


def nest_to_mf(constants, exc_specs, inh_specs, name="cell"):
    """
    Convert NEST-unit E-GLIF ``constants`` (ms, mV, pF, pA, nS/ms) to the SI
    parameter dictionary used by the MF formulas. ``exc_specs`` / ``inh_specs``
    are lists of channels, each a dict with ``Q`` (S), ``K`` (in-degree) and
    ``receptor`` (the NEST ``tau_syn{r}`` / ``E_rev{r}`` that channel uses).
    """
    ms, mV, pF, pA = 1e-3, 1e-3, 1e-12, 1e-12

    def channel(spec, t_key, e_key):
        r = spec["receptor"]
        return {
            "Q": spec["Q"],
            "K": spec["K"],
            t_key: constants[f"tau_syn{r}"] * ms,
            e_key: constants[f"E_rev{r}"] * mV,
        }

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
        "k_1": constants["k_1"] / ms,
        "k_2": constants["k_2"] / ms,
        "k_adap": constants["k_adap"] * 1e-9 / ms,
        "exc": [channel(s, "Te", "Ee") for s in exc_specs],
        "inh": [channel(s, "Ti", "Ei") for s in inh_specs],
    }
    out["Gl"] = out["Cm"] / out["tau_m"]  # S
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
    coefficients ``P`` and the SI parameter set needed to evaluate it. The TF
    has ``n_exc`` excitatory and ``n_inh`` inhibitory inputs; ``tags`` lists the
    excitatory labels first, the inhibitory labels last.
    """

    def __init__(self, cell_name, params, alpha, P, tags, freq_axes=None, mse=None):
        self.cell_name = cell_name
        self.params = params
        self.alpha = float(alpha)
        self.P = np.asarray(P, dtype=float)
        self.tags = list(tags)  # [exc labels..., inh labels...]
        self.n_exc = len(self.params["exc"])
        self.n_inh = len(self.params["inh"])
        if len(self.tags) != self.n_exc + self.n_inh:
            raise ValueError(
                f"tags {self.tags} inconsistent with {self.n_exc} exc + {self.n_inh} inh channels"
            )
        self.freq_axes = None if freq_axes is None else [np.asarray(a) for a in freq_axes]
        self.mse = None if mse is None else float(mse)

    @property
    def exc_tags(self):
        return self.tags[: self.n_exc]

    @property
    def inh_tags(self):
        return self.tags[self.n_exc :]

    def evaluate(self, *rates):
        """
        Analytical output rate. Call with the excitatory rate(s) first, the
        inhibitory rate(s) next and an optional adaptation ``w`` last, e.g.
        ``tf(fe, fi)`` / ``tf(fe, fi, w)`` for a 1-exc/1-inh TF, or
        ``tf(fe1, fe2, fi, w)`` for the two-excitatory Golgi TF.
        """
        n = self.n_exc + self.n_inh
        if len(rates) == n:
            w = 0.0
        elif len(rates) == n + 1:
            w = rates[n]
        else:
            raise TypeError(
                f"{self.cell_name} TF takes {n} rates (+optional w), got {len(rates)}"
            )
        fe = list(rates[: self.n_exc])
        fi = list(rates[self.n_exc : n])
        return analytical_tf_rate(fe, fi, w, self.alpha, self.params, self.P)

    __call__ = evaluate

    def evaluate_on_grid(self, freq_axes=None):
        """Evaluate on the stored (exc..., inh...) grid; returns an array in tags order."""
        axes = freq_axes if freq_axes is not None else self.freq_axes
        if axes is None:
            raise RuntimeError("No freq_axes available; pass one explicitly.")
        grids = np.meshgrid(*[np.asarray(a) for a in axes], indexing="ij")
        return self.evaluate(*grids)

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

        def _clean(v):
            return v.item() if isinstance(v, np.generic) else v

        def _clean_channels(channels):
            return [{k: _clean(v) for k, v in c.items()} for c in channels]

        params = {}
        for k, v in self.params.items():
            if k in ("exc", "inh"):
                params[k] = _clean_channels(v)
            else:
                params[k] = _clean(v)

        syn = {
            tag: {"role": "exc", "Q_nS": c["Q"] * 1e9, "K": c["K"]}
            for tag, c in zip(self.exc_tags, self.params["exc"])
        }
        syn.update(
            {
                tag: {"role": "inh", "Q_nS": c["Q"] * 1e9, "K": c["K"]}
                for tag, c in zip(self.inh_tags, self.params["inh"])
            }
        )
        d = {
            "cell_name": self.cell_name,
            "synaptic_inputs": syn,
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
        """RMSE / MAE / max|err| / R² of the analytical TF against a numerical grid (tags order)."""
        an = self.evaluate_on_grid()
        numerical_mean = np.asarray(numerical_mean)
        r = numerical_mean - an
        return {
            "RMSE": float(np.sqrt((r**2).mean())),
            "MAE": float(np.abs(r).mean()),
            "maxAE": float(np.abs(r).max()),
            "R2": float(1 - (r**2).sum() / ((numerical_mean - numerical_mean.mean()) ** 2).sum()),
        }

    def summary(self):
        print(f"Cell    : {self.cell_name.replace('_', ' ')}")
        print(f"Inputs  : exc={self.exc_tags}, inh={self.inh_tags}")
        print(f"alpha   : {self.alpha:.3f}")
        print("P [mV]  : " + ", ".join(f"P{i}={p * 1e3:.2f}" for i, p in enumerate(self.P)))
        if self.freq_axes is not None:
            for tag, ax in zip(self.tags, self.freq_axes):
                print(f"Domain  : {tag} [{ax[0]:.1f}–{ax[-1]:.1f} Hz]")
        if self.mse is not None:
            print(f"Fit MSE : {self.mse:.3f} Hz^2  (RMSE {np.sqrt(self.mse):.2f} Hz)")

    def plot_tf(self, numerical_tf=None, save=True, n=60, output_dir="."):
        """3D surface of the analytical TF (single excitatory + single inhibitory only)."""
        if self.n_exc != 1 or self.n_inh != 1:
            raise NotImplementedError(
                "plot_tf surface is only defined for 1-exc/1-inh TFs; use the "
                "numerical 3D scatter for TFs with more inputs (e.g. the Golgi cell)."
            )
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

    def __init__(
        self,
        scaffold_path,
        yaml_path,
        tf_dict_path,
        cell_name,
        numerical_tf,
        cut=2,
    ):
        self.cell_name = cell_name
        self.cut = cut
        self.numerical_tf = numerical_tf

        self._convergences, self._synapses_per_pair = self._scaffold_convergences(scaffold_path)
        self.params, self.exc_labels, self.inh_labels = self._build_params(
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

        # group connections by label (see class docstring for the Q/K rule)
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

        for g in groups.values():
            g["E_rev"] = constants[f"E_rev{g['receptor']}"]
            g["Q"] = float(np.mean(g["Qs"])) if g["Qs"] else 0.0
            g["K"] = g["K"] / max(len(g["pops"]), 1)

        # excitatory = reversal above leak (depolarizing), inhibitory otherwise;
        # order by reversal (descending) for a deterministic channel ordering
        El = constants["E_L"]
        labels = sorted(groups, key=lambda l: groups[l]["E_rev"], reverse=True)
        exc_labels = [l for l in labels if groups[l]["E_rev"] > El]
        inh_labels = [l for l in labels if groups[l]["E_rev"] <= El]
        if not exc_labels or not inh_labels:
            raise ValueError(
                f"{cell_name}: expected both exc and inh inputs, got "
                f"exc={exc_labels}, inh={inh_labels}"
            )

        def specs(lbls):
            return [
                {"Q": groups[l]["Q"], "K": groups[l]["K"], "receptor": groups[l]["receptor"]}
                for l in lbls
            ]

        params = nest_to_mf(constants, specs(exc_labels), specs(inh_labels), name=cell_name)
        return params, exc_labels, inh_labels

    def _prepare_template(self):
        ntf = self.numerical_tf
        order = self.exc_labels + self.inh_labels
        idx = [ntf.tags.index(l) for l in order]
        c = self.cut
        axes = [ntf.freq_axes[k][c:] for k in idx]
        n_exc = len(self.exc_labels)
        self.Fe_axes = axes[:n_exc]
        self.Fi_axes = axes[n_exc:]

        sl = tuple(slice(c, None) for _ in idx)
        # reorder numerical axes to [exc..., inh...], drop the cut bins
        self.MEAN = np.moveaxis(ntf.tf_mean, idx, range(len(idx)))[sl]
        self.SD = np.moveaxis(ntf.tf_std, idx, range(len(idx)))[sl]

        grids = np.meshgrid(*axes, indexing="ij")
        self.fe_grids = grids[:n_exc]  # one per excitatory channel
        self.fi_grids = grids[n_exc:]  # one per inhibitory channel
        self.w = np.zeros_like(self.MEAN)

    def _make_result(self, alpha, P, mse):
        return AnalyticalTransferFunction(
            self.numerical_tf.cell_name,
            self.params,
            alpha,
            P,
            self.exc_labels + self.inh_labels,
            freq_axes=self.Fe_axes + self.Fi_axes,
            mse=mse,
        )

    # -- fitting ------------------------------------------------------------ #

    def fit(self, alpha=4.5, save=True, output_dir="output_analytical_TFs", verbose=True):
        """Two-step fit at a fixed gain ``alpha``. Returns an AnalyticalTransferFunction."""
        P, mse = fit_two_step(
            self.MEAN,
            self.fe_grids,
            self.fi_grids,
            self.w,
            self.params,
            alpha,
            verbose=verbose,
        )
        self._fit = self._make_result(alpha, P, mse)
        if save:
            self._fit.save(output_dir)
        return self._fit

    def _scan_alpha(self, alphas):
        out = [
            fit_two_step(
                self.MEAN,
                self.fe_grids,
                self.fi_grids,
                self.w,
                self.params,
                a,
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
            self.MEAN,
            self.fe_grids,
            self.fi_grids,
            self.w,
            self.params,
            alpha_fine[b],
            Ps[b],
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
        """3D surface of the current fit on top of the numerical dots (1-exc/1-inh only)."""
        if self._fit is None:
            raise RuntimeError("Run fit() or optimize_alpha() first.")
        return self._fit.plot_tf(numerical_tf=self.numerical_tf, save=save, output_dir=output_dir)

    def plot_fit(self, save=True, output_dir="."):
        """
        Sections at fixed (primary) inhibition: numerical (diamonds, ±SD) vs
        analytical (lines). Any extra excitatory/inhibitory channels are held at
        their mid-range value; the x axis is the first excitatory channel and the
        colour the first inhibitory channel.
        """
        if self._fit is None:
            raise RuntimeError("Run fit() or optimize_alpha() first.")
        fe_ax = self.Fe_axes[0]
        fi_ax = self.Fi_axes[0]
        fe = np.linspace(fe_ax[0], fe_ax[-1], 200)
        cmap = plt.cm.viridis
        n = len(fi_ax)

        # hold the non-plotted channels at their middle grid value
        extra_exc = [ax[len(ax) // 2] for ax in self.Fe_axes[1:]]
        extra_inh = [ax[len(ax) // 2] for ax in self.Fi_axes[1:]]
        mid_exc = [len(ax) // 2 for ax in self.Fe_axes[1:]]
        mid_inh = [len(ax) // 2 for ax in self.Fi_axes[1:]]

        fig, ax = plt.subplots(figsize=(6, 4.5))
        for j, fi0 in enumerate(fi_ax):
            col = cmap(j / max(n - 1, 1))
            sl = (slice(None), *mid_exc, j, *mid_inh)
            ax.errorbar(
                fe_ax,
                self.MEAN[sl],
                yerr=self.SD[sl],
                color=col,
                marker="D",
                ms=4,
                capsize=2,
                elinewidth=1,
                lw=0,
            )
            ax.plot(fe, self._fit.evaluate(fe, *extra_exc, fi0, *extra_inh), color=col)
        m = self._fit.metrics(self.MEAN)
        ax.set_title(
            f"{self.numerical_tf.cell_name}   α = {self._fit.alpha:.2f}   RMSE = {m['RMSE']:.1f} Hz"
        )
        ax.set_xlabel(f"{self.exc_labels[0]} [Hz]")
        ax.set_ylabel(f"{self.numerical_tf.cell_name.replace('_', ' ')} [Hz]")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(fi_ax[0], fi_ax[-1]))
        fig.colorbar(sm, ax=ax, label=f"{self.inh_labels[0]} [Hz]")
        fig.tight_layout()

        if save:
            os.makedirs(output_dir, exist_ok=True)
            fig.savefig(
                os.path.join(output_dir, f"analyticalTF_fit_{self.numerical_tf.cell_name}.png"),
                dpi=150,
            )
        return fig
