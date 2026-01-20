#!/usr/bin/env python3
from __future__ import annotations

import os
import argparse
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml
import nest
import matplotlib.pyplot as plt


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


def load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        return yaml.safe_load(f)


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
    E_L = params_raw.get("E_L", defaults.get("E_L", None))
    if E_L is not None and "V_m" in allowed:
        params["V_m"] = float(E_L)
    return params


def compute_g_unblocked(I: np.ndarray, V: np.ndarray, Erev: float, mg: Optional[np.ndarray]) -> np.ndarray:
    denom = (Erev - V)
    denom = np.where(np.abs(denom) < 1e-9, np.sign(denom) * 1e-9, denom)
    g = I / denom
    if mg is not None:
        mb = np.where(np.abs(mg) < 1e-6, np.sign(mg) * 1e-6, mg)
        g = g / mb
    return g


@dataclass(frozen=True)
class PlotCfg:
    dt_ms: float = 0.1
    spike_times_ms: Tuple[float, ...] = (10.0, 30.0, 50.0, 70.0, 90.0)
    sim_time_ms: float = 800.0


def simulate_connection(params_post_raw: Dict[str, Any], syn_list: List[Dict[str, Any]], cfg: PlotCfg, silent: bool) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    nest.ResetKernel()
    nest.Install("cerebmodule")
    nest.SetKernelStatus({"resolution": float(cfg.dt_ms)})

    post_params = dict(params_post_raw)
    if silent:
        post_params["I_e"] = 0.0
    post_params = sanitize_eglif_params(post_params, "eglif_multirec_opt")

    pre = nest.Create("parrot_neuron", 1)
    post = nest.Create("eglif_multirec_opt", 1, params=post_params)

    if "E_L" in params_post_raw:
        try:
            nest.SetStatus(post, {"V_m": float(params_post_raw["E_L"])})
        except Exception:
            pass

    sp_gen = nest.Create("spike_generator", params={"spike_times": [30]})#list(map(float, cfg.spike_times_ms))})
    nest.Connect(sp_gen, pre)

    for s in syn_list:
        model = s.get("model")
        if model not in ("static_synapse", "tsodyks2_synapse"):
            continue
        syn_spec: Dict[str, Any] = {"synapse_model": model}
        if "weight" in s:
            syn_spec["weight"] = float(s["weight"])
        if "delay" in s:
            syn_spec["delay"] = float(s["delay"])
        if "receptor_type" in s:
            syn_spec["receptor_type"] = int(s["receptor_type"])
        if model == "tsodyks2_synapse":
            for k in ("U", "x", "tau_rec", "tau_fac"):
                if k in s:
                    syn_spec[k] = float(s[k])
        nest.Connect(pre, post, syn_spec=syn_spec)

    recordables = set(get_recordables("eglif_multirec_opt"))
    record_from = ["V_m"]

    # record all receptor currents present in syn_list
    for s in syn_list:
        rt = int(s.get("receptor_type", -1))
        ivar = RECEPTOR_TO_IVAR.get(rt)
        if ivar and ivar in recordables:
            record_from.append(ivar)

    if MG_BLOCK_VAR in recordables:
        record_from.append(MG_BLOCK_VAR)

    record_from = list(dict.fromkeys(record_from))

    mult = nest.Create("multimeter", params={
        "interval": float(cfg.dt_ms),
        "record_from": record_from,
        "record_to": "memory",
    })
    nest.Connect(mult, post)

    nest.Simulate(float(cfg.sim_time_ms))

    ev = mult.get()["events"]
    t = np.asarray(ev["times"])
    V = np.asarray(ev["V_m"])
    traces: Dict[str, np.ndarray] = {k: np.asarray(ev[k]) for k in record_from if k in ev}
    return t, V, traces


def plot_conn(conn: str, post_key: str, params_post: Dict[str, Any],
              syn_static: List[Dict[str, Any]], syn_tm: List[Dict[str, Any]],
              cfg: PlotCfg, outdir: str, silent: bool) -> None:
    os.makedirs(outdir, exist_ok=True)

    # receptor types present in both
    rts_s = {int(s["receptor_type"]) for s in syn_static if "receptor_type" in s}
    rts_t = {int(s["receptor_type"]) for s in syn_tm if "receptor_type" in s}
    rts = sorted(rts_s & rts_t)
    if not rts:
        return

    tS, VS, trS = simulate_connection(params_post, syn_static, cfg, silent=silent)
    tT, VT, trT = simulate_connection(params_post, syn_tm, cfg, silent=silent)

    nrows = len(rts) + 1
    fig, axs = plt.subplots(nrows, 2, figsize=(18, 3.0 * nrows), sharex=True)

    axs[0, 0].set_title(f"{conn} | STATIC")
    axs[0, 1].set_title(f"{conn} | TSODYKS tuned")

    mgS = trS.get(MG_BLOCK_VAR, None)
    mgT = trT.get(MG_BLOCK_VAR, None)

    # g panels
    for i, rt in enumerate(rts):
        kind = RECEPTOR_KIND.get(rt, f"rt={rt}")
        ivar = RECEPTOR_TO_IVAR.get(rt)
        erev_key = RECEPTOR_TO_EREVKEY.get(rt)
        if ivar is None or erev_key is None or erev_key not in params_post:
            continue
        if ivar not in trS or ivar not in trT:
            continue

        Erev = float(params_post[erev_key])

        gS = compute_g_unblocked(trS[ivar], VS, Erev, mgS if is_nmda(rt) else None)
        gT = compute_g_unblocked(trT[ivar], VT, Erev, mgT if is_nmda(rt) else None)

        axs[i, 0].plot(tS, np.abs(gS))
        axs[i, 1].plot(tT, np.abs(gT))
        axs[i, 0].set_ylabel(f"|g| {kind}")
        axs[i, 1].set_ylabel(f"|g| {kind}")

        for st in cfg.spike_times_ms:
            axs[i, 0].axvline(st, linewidth=0.8, alpha=0.25)
            axs[i, 1].axvline(st, linewidth=0.8, alpha=0.25)

        if is_nmda(rt):
            axs[i, 0].set_ylabel("|g| NMDA (unblocked)")
            axs[i, 1].set_ylabel("|g| NMDA (unblocked)")

    # Vm bottom
    axs[-1, 0].plot(tS, VS)
    axs[-1, 1].plot(tT, VT)
    axs[-1, 0].set_ylabel("V_m [mV]")
    axs[-1, 1].set_ylabel("V_m [mV]")
    axs[-1, 0].set_xlabel("time [ms]")
    axs[-1, 1].set_xlabel("time [ms]")

    for st in cfg.spike_times_ms:
        axs[-1, 0].axvline(st, linewidth=0.8, alpha=0.25)
        axs[-1, 1].axvline(st, linewidth=0.8, alpha=0.25)

    fig.suptitle(f"Post: {post_key} | V_m init = E_L | NMDA always unblocked by Mg_block", y=0.995, fontsize=11)
    fig.tight_layout()

    safe = conn.replace("/", "_")
    fig.savefig(os.path.join(outdir, f"verify_{safe}.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--static", default="../../../configurations/mouse/nest_examples/basal_vitro_rec.yaml")
    ap.add_argument("--tm_tuned", default="../../../configurations/mouse/nest_examples/basal_vitro_rec_tm_tuned.yaml")
    ap.add_argument("--sim", default="nest_basal_activity")
    ap.add_argument("--outdir", default="verify_imgs")
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--sim_time", type=float, default=800.0)
    ap.add_argument("--spikes", type=str, default="10,30,50,70,90")
    ap.add_argument("--keep_Ie", action="store_true")
    args = ap.parse_args()

    spike_times = tuple(float(x) for x in args.spikes.split(",") if x.strip())
    cfg = PlotCfg(dt_ms=args.dt, spike_times_ms=spike_times, sim_time_ms=args.sim_time)

    static_yaml = load_yaml(args.static)
    tm_yaml = load_yaml(args.tm_tuned)

    static_cm = static_yaml["simulations"][args.sim]["connection_models"]
    tm_cm = tm_yaml["simulations"][args.sim]["connection_models"]
    cell_models = static_yaml["simulations"][args.sim]["cell_models"]

    common_conns = sorted(set(static_cm.keys()) & set(tm_cm.keys()))
    print(f"Connessioni comuni: {len(common_conns)}")

    for conn in common_conns:
        synS = static_cm[conn].get("synapses", [])
        synT = tm_cm[conn].get("synapses", [])
        if not synS or not synT:
            continue

        post_key = infer_post_cell_key(conn, cell_models)
        if post_key is None:
            print(f"[SKIP] {conn}: post non inferibile")
            continue

        params_post = get_post_params(static_yaml, args.sim, post_key)
        if params_post is None:
            print(f"[SKIP] {conn}: post '{post_key}' non eglif_multirec_opt")
            continue

        print(f"[PLOT] {conn} -> post={post_key}")
        plot_conn(conn, post_key, params_post, synS, synT, cfg, args.outdir, silent=(not args.keep_Ie))

    print(f"Imgs in: {args.outdir}/")


if __name__ == "__main__":
    main()

