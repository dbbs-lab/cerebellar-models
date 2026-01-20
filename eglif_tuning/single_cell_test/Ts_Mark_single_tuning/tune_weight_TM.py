#!/usr/bin/env python3
from __future__ import annotations

import copy
import argparse
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml
import nest


# ------------------ YAML I/O ------------------

def load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        return yaml.safe_load(f)

def save_yaml(obj: Dict[str, Any], path: str) -> None:
    with open(path, "w") as f:
        yaml.safe_dump(obj, f, sort_keys=False)


# ------------------ Receptor conventions ------------------

RECEPTOR_KIND = {1: "AMPA", 2: "NMDA", 3: "GABA", 4: "AMPA2", 5: "AMPA3"}

RECEPTOR_TO_IVAR = {
    1: "I_syn_ampa",
    2: "I_syn_nmda",
    3: "I_syn_gaba",
    4: "I_syn_ampa2",
    5: "I_syn_ampa3",
}

RECEPTOR_TO_EREVKEY = {
    1: "AMPA_E_rev",
    2: "NMDA_E_rev",
    3: "GABA_E_rev",
    4: "AMPA2_E_rev",
    5: "AMPA3_E_rev",
}

MG_BLOCK_VAR = "Mg_block"

def is_nmda(rt: int) -> bool:
    return int(rt) == 2


# ------------------ Helpers ------------------

def infer_post_cell_key(conn_name: str, cell_models: Dict[str, Any]) -> Optional[str]:
    if "_to_" not in conn_name:
        return None
    post_token = conn_name.split("_to_", 1)[1]
    if post_token in cell_models:
        return post_token
    cand = f"{post_token}_cell"
    if cand in cell_models:
        return cand
    for k in cell_models.keys():
        if k.startswith(post_token):
            return k
    return None

def get_post_params(static_yaml: Dict[str, Any], sim_name: str, post_key: str) -> Optional[Dict[str, Any]]:
    cm = static_yaml["simulations"][sim_name]["cell_models"][post_key]
    if cm.get("model") != "eglif_multirec_opt":
        return None
    return dict(cm.get("constants", {}))

def index_synapses_by_receptor(syn_list: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for s in syn_list:
        if "receptor_type" in s:
            out[int(s["receptor_type"])] = s
    return out

def get_recordables(model_name: str) -> List[str]:
    return list(nest.GetDefaults(model_name).get("recordables", []))

def sanitize_eglif_params(params_raw: Dict[str, Any], model_name: str) -> Dict[str, Any]:
    defaults = nest.GetDefaults(model_name)
    allowed = set(defaults.keys())
    params: Dict[str, Any] = {}
    for k, v in params_raw.items():
        if isinstance(v, dict):
            continue
        if k in allowed:
            params[k] = v
    # Force Vm init = EL (requested)
    E_L = params_raw.get("E_L", defaults.get("E_L", None))
    if E_L is not None and "V_m" in allowed:
        params["V_m"] = float(E_L)
    return params

def compute_g_unblocked(I: np.ndarray, V: np.ndarray, Erev: float, mg: Optional[np.ndarray]) -> np.ndarray:
    # g = I/(Erev - V); if NMDA: g_unblocked = g / Mg_block
    denom = (Erev - V)
    denom = np.where(np.abs(denom) < 1e-9, np.sign(denom) * 1e-9, denom)
    g = I / denom
    if mg is not None:
        mb = np.where(np.abs(mg) < 1e-6, np.sign(mg) * 1e-6, mg)
        g = g / mb
    return g

def peak_abs_in_window(times: np.ndarray, x: np.ndarray, t0: float, t1: float) -> float:
    m = (times >= t0) & (times <= t1)
    if not np.any(m):
        return 0.0
    return float(np.max(np.abs(x[m])))


# ------------------ Simulation ------------------

@dataclass(frozen=True)
class TuneConfig:
    dt_ms: float = 0.1
    spike_times_ms: Tuple[float, ...] = (10.0, 30.0, 50.0, 70.0, 90.0)
    sim_time_ms: float = 400.0
    fast_window_ms: float = 8.0          # AMPA/GABA windows
    nmda_window_ms: float = 250.0        # after LAST spike


def simulate_isolated_receptor(
        params_post_raw: Dict[str, Any],
        syn: Dict[str, Any],
        cfg: TuneConfig,
        silent: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """
    One receptor at a time:
      spike_generator(5 spikes) -> parrot -> post
    Post forced to V_m=E_L always; optionally I_e=0 to avoid intrinsic firing.
    """
    nest.ResetKernel()
    nest.Install("cerebmodule")
    nest.SetKernelStatus({"resolution": float(cfg.dt_ms)})

    # Build post params
    post_params = dict(params_post_raw)
    if silent:
        post_params["I_e"] = 0.0
    post_params = sanitize_eglif_params(post_params, "eglif_multirec_opt")

    pre = nest.Create("parrot_neuron", 1)
    post = nest.Create("eglif_multirec_opt", 1, params=post_params)

    # hard-set V_m to E_L after creation too
    if "E_L" in params_post_raw:
        try:
            nest.SetStatus(post, {"V_m": float(params_post_raw["E_L"])})
        except Exception:
            pass

    sp_gen = nest.Create("spike_generator", params={"spike_times": [30.]}) #list(map(float, cfg.spike_times_ms))})
    nest.Connect(sp_gen, pre)

    model = syn.get("model")
    syn_spec: Dict[str, Any] = {"synapse_model": model}

    for k in ("weight", "delay"):
        if k in syn:
            syn_spec[k] = float(syn[k])

    if "receptor_type" in syn:
        syn_spec["receptor_type"] = int(syn["receptor_type"])

    if model == "tsodyks2_synapse":
        for k in ("U", "x", "tau_rec", "tau_fac"):
            if k in syn:
                syn_spec[k] = float(syn[k])

    nest.Connect(pre, post, syn_spec=syn_spec)

    recordables = set(get_recordables("eglif_multirec_opt"))
    record_from = ["V_m"]

    rt = int(syn.get("receptor_type", -1))
    ivar = RECEPTOR_TO_IVAR.get(rt)
    if ivar and ivar in recordables:
        record_from.append(ivar)

    # Mg_block recorded always if available, but used only for NMDA
    if MG_BLOCK_VAR in recordables:
        record_from.append(MG_BLOCK_VAR)

    mult = nest.Create("multimeter", params={
        "interval": float(cfg.dt_ms),
        "record_from": record_from,
        "record_to": "memory",
    })
    nest.Connect(mult, post)

    nest.Simulate(float(cfg.sim_time_ms))

    ev = mult.get()["events"]
    times = np.asarray(ev["times"])
    V = np.asarray(ev["V_m"])
    traces: Dict[str, np.ndarray] = {k: np.asarray(ev[k]) for k in record_from if k in ev}
    return times, V, traces


def tune_tm_weights_align(
        static_yaml: Dict[str, Any],
        tm_yaml: Dict[str, Any],
        sim_name: str,
        cfg: TuneConfig,
        silent_post: bool = True,
        verbose: bool = True,
) -> Dict[str, Any]:
    new_tm = copy.deepcopy(tm_yaml)

    static_cm = static_yaml["simulations"][sim_name]["connection_models"]
    tm_cm = new_tm["simulations"][sim_name]["connection_models"]
    cell_models = static_yaml["simulations"][sim_name]["cell_models"]

    common_conns = sorted(set(static_cm.keys()) & set(tm_cm.keys()))
    if verbose:
        print(f"Common connections: {len(common_conns)}")

    for conn in common_conns:
        synS_all = static_cm[conn].get("synapses", [])
        synT_all = tm_cm[conn].get("synapses", [])
        if not synS_all or not synT_all:
            continue

        idxS = index_synapses_by_receptor(synS_all)
        idxT = index_synapses_by_receptor(synT_all)
        rts = sorted(set(idxS.keys()) & set(idxT.keys()))
        if not rts:
            continue

        post_key = infer_post_cell_key(conn, cell_models)
        if post_key is None:
            if verbose:
                print(f"[SKIP] {conn}: post non inferibile")
            continue

        params_post = get_post_params(static_yaml, sim_name, post_key)
        if params_post is None:
            if verbose:
                print(f"[SKIP] {conn}: post '{post_key}' non eglif_multirec_opt")
            continue

        for rt in rts:
            rt = int(rt)
            kind = RECEPTOR_KIND.get(rt, f"rt={rt}")

            # isolate THIS receptor only
            synS = dict(idxS[rt])
            synT = dict(idxT[rt])

            # sanity: receptor must exist in post params
            erev_key = RECEPTOR_TO_EREVKEY.get(rt)
            ivar = RECEPTOR_TO_IVAR.get(rt)
            if erev_key is None or ivar is None or erev_key not in params_post:
                if verbose:
                    print(f"[SKIP] {conn} {kind}: manca {erev_key} o mapping")
                continue

            # simulate static + tm (isolated)
            tS, VS, trS = simulate_isolated_receptor(params_post, synS, cfg, silent=silent_post)
            tT, VT, trT = simulate_isolated_receptor(params_post, synT, cfg, silent=silent_post)

            if ivar not in trS or ivar not in trT:
                if verbose:
                    print(f"[SKIP] {conn} {kind}: '{ivar}' non registrato")
                continue

            mgS = trS.get(MG_BLOCK_VAR, None) if is_nmda(rt) else None
            mgT = trT.get(MG_BLOCK_VAR, None) if is_nmda(rt) else None

            if is_nmda(rt):
                if (MG_BLOCK_VAR not in trS) or (MG_BLOCK_VAR not in trT):
                    if verbose:
                        print(f"[SKIP] {conn} NMDA: Mg_block non disponibile")
                    continue

            Erev = float(params_post[erev_key])

            gS = compute_g_unblocked(trS[ivar], VS, Erev, mgS)
            gT = compute_g_unblocked(trT[ivar], VT, Erev, mgT)

            # define target windows:
            first_spike = float(cfg.spike_times_ms[0])
            last_spike = float(cfg.spike_times_ms[-1])

            if is_nmda(rt):
                # NMDA: match peak after LAST spike (stable accumulation proxy)
                peakS = peak_abs_in_window(tS, gS, last_spike, last_spike + cfg.nmda_window_ms)
                peakT = peak_abs_in_window(tT, gT, last_spike, last_spike + cfg.nmda_window_ms)
            else:
                # fast receptors: match peak after FIRST spike
                peakS = peak_abs_in_window(tS, gS, first_spike, first_spike + cfg.fast_window_ms)
                peakT = peak_abs_in_window(tT, gT, first_spike, first_spike + cfg.fast_window_ms)

            if not np.isfinite(peakS) or not np.isfinite(peakT) or peakS <= 0 or peakT <= 0:
                if verbose:
                    print(f"[WARN] {conn} {kind}: peakS={peakS:.3e}, peakT={peakT:.3e} -> skip")
                continue

            old_w = float(idxT[rt].get("weight", 1.0))
            new_w = old_w * (peakS / peakT)
            idxT[rt]["weight"] = float(new_w)

            if verbose:
                extra = ""
                if synT.get("model") == "tsodyks2_synapse" and "U" in synT:
                    extra = f" (U={float(synT['U']):.3g}, 1/U~{1/float(synT['U']):.3g})"
                print(f"[OK] {conn} {kind} rt={rt}: w {old_w:.6g} -> {new_w:.6g} | "
                      f"peakS={peakS:.3e}, peakT={peakT:.3e}{extra}")

    return new_tm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--static", default="../../../configurations/mouse/nest_examples/basal_vitro_rec.yaml")
    ap.add_argument("--tm", default="../../../configurations/mouse/nest_examples/basal_vitro_rec_tm.yaml")
    ap.add_argument("--out", default="../../../configurations/mouse/nest_examples/basal_vitro_rec_tm_tuned.yaml")
    ap.add_argument("--sim", default="nest_basal_activity")
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--spikes", type=str, default="10,30,50,70,90")
    ap.add_argument("--sim_time", type=float, default=400.0)
    ap.add_argument("--fast_win", type=float, default=8.0)
    ap.add_argument("--nmda_win", type=float, default=250.0)
    ap.add_argument("--keep_Ie", action="store_true")
    args = ap.parse_args()

    spike_times = tuple(float(x) for x in args.spikes.split(",") if x.strip())
    cfg = TuneConfig(
        dt_ms=args.dt,
        spike_times_ms=spike_times,
        sim_time_ms=args.sim_time,
        fast_window_ms=args.fast_win,
        nmda_window_ms=args.nmda_win,
    )

    static_yaml = load_yaml(args.static)
    tm_yaml = load_yaml(args.tm)

    tuned = tune_tm_weights_align(
        static_yaml=static_yaml,
        tm_yaml=tm_yaml,
        sim_name=args.sim,
        cfg=cfg,
        silent_post=(not args.keep_Ie),
        verbose=True,
    )

    save_yaml(tuned, args.out)
    print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
