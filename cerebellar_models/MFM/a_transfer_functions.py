import json
import nest
import numpy as np
import yaml
import matplotlib.pyplot as plt
import scipy.special as sp_spec
from scipy.optimize import minimize
from bsb import from_storage
from cerebellar_models.analysis import StructureReport
from cerebellar_models.MFM.transfer_functions import TransferFunction, TransferFunctionSimulator
from matplotlib.colors import LinearSegmentedColormap


class TransferFunctionFitting:
    """
    Computes the transfer function of a single-compartment NEST neuron
    driven by Poisson inputs, given a BSB scaffold and circuit configuration.
    Supports multiprocessing: set ``n_workers`` to use a parallel pool.
    """

    def __init__(
        self,
        scaffold_path,
        yaml_path,
        tf_dict_path,
        cell_name,
        num_tf,
        n_points=20,
        params = {}
    ):

        scaffold = from_storage(scaffold_path)
        report = StructureReport(scaffold)
        conn_table = report.plots['connectivity_table']
        conn_table.update()

        self._convergences = {
            ps.tag: (np.mean(conv_arr), np.std(conv_arr))
            for ps, conv_arr in zip(
                scaffold.get_connectivity_sets(),
                conn_table.get_convergences().values(),
            )
        }
        self._synapses_per_pair = {
            ps.tag: np.mean(syn_arr)
            for ps, syn_arr in zip(
                scaffold.get_connectivity_sets(),
                conn_table.get_nb_synapse_per_pair().values(),
            )
        }
        self.scaffold_path = scaffold_path
        self.cell_name = cell_name
        self.sim_info = self._build_sim_info(yaml_path, tf_dict_path, cell_name)
        self.num_tf = self._load_num_tf(cell_name)
        self.n_points = n_points
        self.params = params

    def _load_num_tf(self, cell_name, output_dir="output_TF"):
        """Load the saved numerical TF matrix from ``output_TF/TF_{cell_name}.npz``."""
        tf = TransferFunction.load(cell_name, output_dir=output_dir)
        return tf.tf_mean

    @staticmethod
    def _extract_reversal_map(cell_params):
        """Build receptor_type -> E_rev mapping from cell constants."""
        if isinstance(cell_params.get("Erev"), dict):
            return {int(k): float(v) for k, v in cell_params["Erev"].items()}

        rev_map = {}
        for key, value in cell_params.items():
            if key.startswith("E_rev"):
                try:
                    receptor_id = int(key.split("E_rev", 1)[1])
                except ValueError:
                    continue
                rev_map[receptor_id] = float(value)

        if not rev_map and "E_e" in cell_params and "E_i" in cell_params:
            rev_map[1] = float(cell_params["E_e"])
            rev_map[2] = float(cell_params["E_i"])

        if not rev_map:
            raise KeyError(
                "No reversal potentials found in cell parameters; expected E_rev1/E_rev2/E_rev3 or E_e/E_i."
            )
        return rev_map

    _CELL_COLORS = {
        "granule_cell":   [0.7, 0.15, 0.15, 1.0],
        "golgi_cell":     [0, 0.45, 0.7, 1.0],
        "purkinje_cell":  [0.275, 0.800, 0.275, 1.0],
        "basket_cell":    [1, 0.647, 0, 1.0],
        "stellate_cell":  [1, 0.84, 0, 1.0],
        "dcn_p":          [0.3, 0.3, 0.3, 1.0],
        "dcn_i":          [0.635, 0, 0.145, 1.0],
        "io":             [0.46, 0.376, 0.54, 1.0],
    }

    def _cell_cmap(self):
        color = self._CELL_COLORS.get(self.cell_name, [0, 0, 0, 1.0])
        return LinearSegmentedColormap.from_list(self.cell_name, [(1, 1, 1, 1), color])


    def _build_sim_info(self, yaml_path, tf_dict_path, cell_name):
        simulator = TransferFunctionSimulator(
            self.scaffold_path,
            yaml_path,
            tf_dict_path,
            cell_name,
            n_points=1,
            n_reps=1,
            duration=1.0,
            transient=0.0,
            base_seed=1234,
            n_workers=1,
            show_rasters=False,
            raster_max_points=0,
            stochastic_convergence=False,
        )
        sim_info = simulator.sim_info.copy()
        sim_info['cell_params'] = simulator.sim_info['cell_params'].copy()

        if 'Erev' not in sim_info['cell_params']:
            sim_info['cell_params']['Erev'] = self._extract_reversal_map(
                sim_info['cell_params']
            )

        return sim_info


    def _setup_grid(self):
        label_to_range = {}
        for info in self.sim_info['connections'].values():
            if info['label'] not in label_to_range:
                label_to_range[info['label']] = info['rate_range']
        unique_labels = list(label_to_range.keys())
        label_to_idx = {lbl: i for i, lbl in enumerate(unique_labels)}
        freq_axes = [
            np.linspace(label_to_range[l][0], label_to_range[l][1], self.n_points)
            for l in unique_labels
        ]
        return unique_labels, label_to_idx, freq_axes


    def _rates_dict(self, idx, label_to_idx, freq_axes):
        return {
            conn_tag: freq_axes[label_to_idx[info['label']]][idx[label_to_idx[info['label']]]]
            for conn_tag, info in self.sim_info['connections'].items()
        }
    

    ## WANNA WORK WITH MATRICES
    def _rates_grid(self):
        """Build broadcasted rate arrays for every connection tag over the full TF grid."""
        unique_labels, label_to_idx, freq_axes = self._setup_grid()
        shape = tuple(len(ax) for ax in freq_axes)
        grid_axes = np.meshgrid(*freq_axes, indexing='ij')
        rates_grid = {
            conn_tag: grid_axes[label_to_idx[info['label']]]
            for conn_tag, info in self.sim_info['connections'].items()
        }
        return rates_grid, unique_labels, label_to_idx, freq_axes
    

    def _compute_muV(self, rates_grid, XX=0.0):
        """Compute muV for a full grid.

        rates_grid: dict mapping connection tag -> ndarray of firing rates (broadcasted
                    to the full TF grid shape).

        Returns: (muV, muGe_tot, muGi_tot, muG)
        """
        
        cell_params = self.sim_info['cell_params']
        #print(self.sim_info)
        #exit
        Gl = cell_params['C_m'] / cell_params['tau_m']
        El = cell_params['E_L']
        E_rev = cell_params['Erev']

        muV_terms = []
        muG_terms = []

        for tag, info in self.sim_info['connections'].items():
            ## CASE OF PF THIS MUST BE:
            ## (Q_ascending_axons*K_ascending_axons*T_ascending_axons + Q_parallel_fibers*K_parallel_fibers*T_parallel_fibers)*f_GrC +
            ## (Q_basket*K_basket*T_basket + Q_stellate*K_stellate*T_stellate)*f_MLI
            f = np.array(rates_grid[tag], dtype=float)*10**-3
            #print(np.shape(f), f)
            
            K = float(info.get('convergence', None))
            Q = float(info.get('weight', None))
            T = float(info.get('delay', None))
            receptor_type = int(info.get('receptor_type', None))

            #MuG single pop
            term_muG = Q * K * T * f
            muG_terms.append(term_muG)

            ## muG*Erev for the correct receptor_type
            term_muv = term_muG * E_rev.get(receptor_type, E_rev.get(1, 0.0))
            muV_terms.append(term_muv)

        muV_tot = np.sum(muV_terms, axis=0)
        muG = Gl + np.sum(muG_terms, axis=0)

        # compute muV (no adaptation term)
        muV = (np.e * (muV_tot + Gl * El) - XX) / (muG + 1e-20)
        # print(muV)

        return muV, muG
    

    def _compute_sigmaV(self, rates_grid, muV, muG):
        """Compute sigmaV for full grid.

        Returns sV (ndarray).
        """

        cell_params = self.sim_info['cell_params']
        E_rev = cell_params['Erev']

        Tm = cell_params['C_m'] / (muG + 1e-20)

        sV_terms = []

        for tag, info in self.sim_info['connections'].items():
            f = np.array(rates_grid[tag], dtype=float)*10**-3
            K = float(info.get('convergence', 0.0))
            Q = float(info.get('weight', 0.0))
            T = float(info.get('delay', 1.0))
            receptor_type = int(info.get('receptor_type', 1))

            U = Q / (muG + 1e-20) * (E_rev.get(receptor_type, E_rev.get(1, 0.0)) - muV)
            s = (2 * Tm + T) * ((np.e * U * T) / (2 * (T + Tm + 1e-20))) ** 2 * K * f

            sV_terms.append(s)
        
        sV_sq = np.sum(sV_terms, axis=0)
        sV = np.sqrt(sV_sq) + 1e-20

        return sV
    

    def _compute_tauV(self, rates_grid, muV, muG, sV):
        """Compute Tv and normalized TvN for full grid.

        Returns: TvN, Tv
        """
        cell_params = self.sim_info['cell_params']
        Cm = cell_params['C_m']
        Gl = Cm / cell_params['tau_m']
        E_rev = cell_params['Erev']

        Tv_num_terms = []

        for tag, info in self.sim_info['connections'].items():
            f = np.array(rates_grid[tag], dtype=float)*10**-3
            #print(np.shape(f), np.max(f))

            K = float(info.get('convergence', 0.0))
            Q = float(info.get('weight', 0.0))
            T = float(info.get('delay', 1.0))
            receptor_type = int(info.get('receptor_type', 1))

            U = Q / (muG + 1e-20) * (E_rev.get(receptor_type, E_rev.get(1, 0.0)) - muV)
            tv = K * f * U ** 2 * T ** 2 * np.e ** 2

            print('EREV', E_rev.get(receptor_type, E_rev.get(1, 0.0)))

            Tv_num_terms.append(tv)

        Tv_num = np.sum(Tv_num_terms, axis=0)
        Tv = 0.5 * Tv_num / ((sV + 1e-20) ** 2)
        TvN = Tv * Gl / Cm

        return TvN, Tv
    
    def _erfc_func(self, muV, sV, TvN, Vthre, Gl, Cm, alpha):
        """
        Error function-based analytical expression for the TF (eq. 19 in Zerlaut et al. 2016).
        Compute Fout_th from muV, sV, TvN, Vthre using the erfc function. 
        The expression is derived from the first-passage time problem for a leaky integrate-and-fire neuron 
        with Gaussian voltage fluctuations, where the firing rate is proportional to the probability of the 
        voltage crossing the threshold Vthre, which is given by the complementary error function (erfc) of 
        the normalized distance between muV and Vthre, scaled by the voltage noise sV.
        """
        return .5 / TvN * Gl / Cm * (sp_spec.erfc((Vthre - muV) / np.sqrt(2) / sV)) * alpha

    def _effective_Vthre(self, Fout, muV, sV, TvN, Gl, Cm, alpha):
        """
        Effective threshold computed by inverting the erf-based analytical expression for the TF (eq. 19 in Zerlaut et al. 2016).
        Compute effective threshold Vthre_eff from Fout using the inverse of the erfc function.
        """
        return muV + np.sqrt(2) * sV * sp_spec.erfcinv((1 / alpha) * (Fout * 2 * TvN * Cm / Gl))

    def _threshold_func(self, muV, sV, TvN, muGn, P0, P1, P2, P3, P4):
        """
        Polynomial exrpession of the phenomenological threshold function (eq. 25 in Zerlaut et al. 2016)
        """
        muV0, DmuV0 = -60e-3, 10e-3
        sV0, DsV0 = 4e-3, 6e-3
        TvN0, DTvN0 = 0.5, 1.0
        return P0 + P1 * (muV - muV0) / DmuV0 + \
               P2 * (sV - sV0) / DsV0 + P3 * (TvN - TvN0) / DTvN0 + P4 * np.log(muGn)

    def _TF_my_templateup_eglif(self, rates_grid, XX, alpha, P0, P1, P2, P3, P4):
        """
        Defining the template for the anlaytical TF shapes.
        Paarams P iteratively mdoifed during 2nd step fitting procedure to fit the TF shape to the data.
        """
        muV, muG = self._compute_muV(rates_grid, XX=XX)
        sV = self._compute_sigmaV(rates_grid, muV, muG)
        TvN, _ = self._compute_tauV(rates_grid, muV, muG, sV)
        sV = np.maximum(sV, 1e-4)
        cell_params = self.sim_info['cell_params']
        Gl = cell_params['C_m'] / cell_params['tau_m']
        Cm = cell_params['C_m']
        Vthre = self._threshold_func(muV, sV, TvN, muG / np.maximum(Gl, 1e-20), P0, P1, P2, P3, P4)
        Fout_th = self._erfc_func(muV, sV, TvN, Vthre, Gl, Cm, alpha)
        return np.maximum(Fout_th, 1e-8)
    

    def fit_tf(self, alpha, XX=0.0, maxiter=50000, xtol=1e-5):
        """
        FITTING OF THE ANALYTICAL TRANSFER FUNCTION
        Fit TF coefficients using the same two-step procedure as developed in Zerlaut et al 2016.
        Returns the fitted parameters P0...P4 of the phenomenological threshold function
        which is used to compute the effective threshold Vthre_eff and the final TF shape.
        """

        rates_grid, unique_labels, label_to_idx, freq_axes = self._rates_grid()
        muV, muG = self._compute_muV(rates_grid, XX=XX)
        sV = self._compute_sigmaV(rates_grid, muV, muG)
        TvN, _ = self._compute_tauV(rates_grid, muV, muG, sV)

        cell_params = self.sim_info['cell_params']
        Gl = cell_params['C_m'] / cell_params['tau_m']
        Cm = cell_params['C_m']
        
        foutlim = np.nanmax(self.num_tf)
        print(np.shape(self.num_tf))
        i_non_zeros = np.where((self.num_tf > 0.0) & (self.num_tf < foutlim))
        Vthre_eff = self._effective_Vthre(self.num_tf[i_non_zeros], muV[i_non_zeros], sV[i_non_zeros], TvN[i_non_zeros], Gl, Cm, alpha)
        P0 = np.array([-50, 0.0, 0.0, 0.0, 0.0])

        def res_vthre(p):
            """
            1-Residual function for the first step of the fitting procedure, 
            which fits the effective threshold Vthre_eff with the phenomenological threshold function (eq. 25 in Zerlaut et al. 2016) 
            to extract the parameters P0...P4.
            """
            return np.mean((Vthre_eff - self._threshold_func(muV[i_non_zeros], sV[i_non_zeros], TvN[i_non_zeros], muG[i_non_zeros] / Gl, *p)) ** 2)

        plsq = minimize(res_vthre, P0, options={'disp': True})
        P = plsq.x

        def res_fout(p):
            """
            2-Residual function for the second step of the fitting procedure, 
            which fits the TF shape by minimizing the residual between Fout (numerical template) and 
            the TF template with the fitted parameters P.
            """
            return np.mean((self.num_tf - self._TF_my_templateup_eglif(rates_grid, XX, alpha, *p)) ** 2)

        plsq = minimize(res_fout, P, method='nelder-mead', options={'xtol': xtol, 'disp': True, 'maxiter': maxiter})
        P = plsq.x
        print(P)
        self.params['P'] = P
        return P
