#!/usr/bin/env python3
"""Admission ledger for a phytochemical score into a tip-demonstration ODE.

The object is a rule. Thesis #16 claims are transcribed and not re-ranked.
Thesis #7 Fisher ranks are not recomputed. Kinetic Theta and the forcing
schedule are separate hashed objects. A refused call assigns neither.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
CANDIDATE_PATH = ROOT / "candidates.json"
# SHA-256 of sim/candidates.json as stored. A byte edit raises.
CANDIDATE_PIN = "6252fcb50627caba3e0f8eb99043ca52535c59e62350f0ed27dc3fd6653f34c4"

FORCING_SYMBOL = "u_phyto"
PROTOCOL_AMPLITUDE = "1"
HORIZON = "40"
THETA_NAMES = ("d0", "d_u", "g", "k_glyc", "s", "u_in")

# Declared kinetic vector. Decimal strings are the hashed object.
THETA = {
    "d0": "0.15",
    "d_u": "0.25",
    "g": "0.20",
    "k_glyc": "0.60",
    "s": "0.40",
    "u_in": "0.50",
}

G0, R0, A0 = 1.0, 0.20, 0.80
T_END = 40.0
N_SAMPLES = 401

# Affine map used only as a refused contrast: u = ALPHA0 + ALPHA1 * (-score).
AFFINE_ALPHA0 = 0.20
AFFINE_ALPHA1 = 0.50
# Normal-normal soft prior used only as a refused contrast.
PRIOR_SHIFT = 0.02
PRIOR_SD = 0.05
OBS_SD = 0.05


def canon_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def digest(obj) -> str:
    return hashlib.sha256(canon_bytes(obj)).hexdigest()


def load_candidates() -> dict:
    raw = CANDIDATE_PATH.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != CANDIDATE_PIN:
        raise SystemExit(f"candidate pin mismatch: {got}")
    return json.loads(raw.decode("utf-8"))


def theta_float() -> dict:
    return {k: float(THETA[k]) for k in THETA_NAMES}


def null_forcing() -> dict:
    return {
        "amplitude_rule": "null",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": "0"}],
        "provenance": None,
        "symbol": FORCING_SYMBOL,
    }


def legal_forcing(row: dict, library_sha: str) -> dict:
    return {
        "amplitude_rule": "protocol_constant",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": PROTOCOL_AMPLITUDE}],
        "provenance": {
            "claim": int(row["claim"]),
            "library_sha256": library_sha,
            "pains": int(row["pains"]),
            "row_id": row["id"],
            "score_copied": False,
        },
        "symbol": FORCING_SYMBOL,
    }


def row_by_id(doc: dict) -> dict:
    return {row["id"]: row for row in doc["rows"]}


def failed_checks(doc: dict, call: dict, slot_full: bool) -> list[str]:
    """Return every failed predicate. Empty means the call may be admitted."""
    reasons: list[str] = []
    kind = call["kind"]
    destination = call["destination"]
    rule = call["amplitude_rule"]

    if kind == "anecdote":
        reasons.append("not_a_score_record")
    if kind == "utility":
        reasons.append("utility_not_an_input")
    if rule == "soft_prior" or kind == "soft_prior":
        reasons.append("soft_prior")
    if rule == "soft_weight" or kind == "soft_weight":
        reasons.append("soft_weight")
    if destination in THETA_NAMES:
        reasons.append("theta_destination")
    elif destination != FORCING_SYMBOL:
        reasons.append("symbol")

    if kind == "anecdote":
        # The anecdote has no row-level record to audit further.
        return reasons

    row = row_by_id(doc)[call["row_id"]]
    library = call.get("library_sha256", doc["library_sha256"])
    engine = call.get("engine", doc["engine"])
    unit = call.get("unit", doc["unit"])
    may = call.get("may_enter_theta", doc["may_enter_theta"])
    claim = call.get("claim", row["claim"])
    pains = call.get("pains", row["pains"])

    if library != doc["library_sha256"]:
        reasons.append("library_sha")
    if engine is not None:
        reasons.append("engine")
    if unit != "arbitrary_surrogate":
        reasons.append("unit")
    if may is not False:
        reasons.append("may_enter_theta")
    if int(claim) != 3:
        reasons.append("claim")
    if int(pains) != 0:
        reasons.append("pains")
    if rule != "protocol_constant":
        reasons.append("amplitude_rule")

    if not reasons and slot_full:
        reasons.append("slot_occupied")
    return reasons


def apply_call(doc: dict, ledger: dict, call: dict) -> dict:
    reasons = failed_checks(doc, call, slot_full=ledger["forcing"]["amplitude_rule"] != "null")
    record = {
        "id": call["id"],
        "kind": call["kind"],
        "row_id": call.get("row_id"),
        "destination": call["destination"],
        "amplitude_rule": call["amplitude_rule"],
        "reasons": reasons,
        "status": "refused" if reasons else "admitted",
        "theta_sha256": digest(ledger["theta"]),
    }
    if record["status"] == "admitted":
        row = row_by_id(doc)[call["row_id"]]
        ledger["forcing"] = legal_forcing(row, doc["library_sha256"])
        record["admitted_row"] = row["id"]
    record["theta_sha256_after"] = digest(ledger["theta"])
    record["forcing_sha256_after"] = digest(ledger["forcing"])
    record["theta_unchanged"] = record["theta_sha256"] == record["theta_sha256_after"]
    ledger["calls"].append(record)
    return record


def equilibrium(theta: dict, u: float) -> dict:
    g_star = theta["u_in"] / theta["k_glyc"]
    r_star = theta["g"] * theta["u_in"] / (theta["d0"] + theta["d_u"] * u)
    a_star = theta["u_in"] / (theta["s"] * (1.0 + r_star))
    return {"G": g_star, "R": r_star, "A": a_star}


def closed_form(t: np.ndarray, theta: dict, u: float) -> tuple[np.ndarray, np.ndarray]:
    a = theta["k_glyc"]
    g_inf = theta["u_in"] / a
    g = g_inf + (G0 - g_inf) * np.exp(-a * t)
    b = theta["g"] * a
    c = theta["d0"] + theta["d_u"] * u
    if abs(c - a) < 1e-12:
        raise SystemExit("closed form assumes clearance != k_glyc")
    coef = b * (G0 - g_inf) / (c - a)
    particular = b * g_inf / c
    k_hom = R0 - coef - particular
    r = coef * np.exp(-a * t) + particular + k_hom * np.exp(-c * t)
    return g, r


def integrate(theta: dict, u: float) -> dict:
    t = np.linspace(0.0, T_END, N_SAMPLES)

    def field(_t, y):
        g, r, a_state = y
        d_g = theta["u_in"] - theta["k_glyc"] * g
        d_r = theta["g"] * theta["k_glyc"] * g - (theta["d0"] + theta["d_u"] * u) * r
        d_a = theta["k_glyc"] * g / (1.0 + r) - theta["s"] * a_state
        return (d_g, d_r, d_a)

    sol = solve_ivp(
        field,
        (0.0, T_END),
        (G0, R0, A0),
        t_eval=t,
        method="LSODA",
        rtol=1e-8,
        atol=1e-11,
    )
    if not sol.success:
        raise SystemExit(f"integration failed: {sol.message}")
    g_hat, r_hat = closed_form(sol.t, theta, u)
    g_err = float(np.max(np.abs(sol.y[0] - g_hat)))
    r_err = float(np.max(np.abs(sol.y[1] - r_hat)))
    if g_err > 1e-7 or r_err > 1e-7:
        raise SystemExit(f"closed form gap G {g_err} R {r_err}")
    if float(np.min(sol.y[1])) <= 0.0:
        raise SystemExit("stress coordinate left the positive half-line")
    packed = [
        [round(float(sol.t[i]), 10), round(float(sol.y[0, i]), 12), round(float(sol.y[1, i]), 12), round(float(sol.y[2, i]), 12)]
        for i in range(sol.t.size)
    ]
    eq = equilibrium(theta, u)
    terminal = {"G": float(sol.y[0, -1]), "R": float(sol.y[1, -1]), "A": float(sol.y[2, -1])}
    eq_gap = {k: abs(terminal[k] - eq[k]) for k in ("G", "R", "A")}
    return {
        "t": sol.t,
        "y": sol.y,
        "state_sha256": digest(packed),
        "terminal": terminal,
        "equilibrium": eq,
        "equilibrium_abs_gap": eq_gap,
        "closed_form_max_abs": {"G": g_err, "R": r_err},
        "u": u,
    }


def soft_prior_contrast(score: float, d0: float) -> dict:
    """Normal-normal posterior mean. Computed in a side buffer. Not written."""
    mu = d0 + PRIOR_SHIFT * (-score)
    # Equal variances => posterior mean is the average of prior mean and observation.
    y_obs = d0
    post = (mu / PRIOR_SD**2 + y_obs / OBS_SD**2) / (1.0 / PRIOR_SD**2 + 1.0 / OBS_SD**2)
    leaked = dict(THETA)
    leaked["d0"] = f"{post:.16f}"
    return {
        "prior_mean": mu,
        "observation": y_obs,
        "posterior_mean_d0": post,
        "theta_sha256_if_written": digest(leaked),
        "written_to_ledger_theta": False,
    }


def affine_amplitude(score: float) -> float:
    return AFFINE_ALPHA0 + AFFINE_ALPHA1 * (-score)


def eligibility(doc: dict) -> list[dict]:
    """Predicate on each transcribed row, with the legal rule and an empty slot."""
    out = []
    for row in doc["rows"]:
        call = {
            "id": f"predicate-{row['id']}",
            "kind": "score",
            "row_id": row["id"],
            "destination": FORCING_SYMBOL,
            "amplitude_rule": "protocol_constant",
        }
        reasons = failed_checks(doc, call, slot_full=False)
        out.append({"id": row["id"], "label": row["label"], "claim": row["claim"], "pains": row["pains"], "eligible": not reasons, "reasons": reasons})
    return out


def build_calls(doc: dict) -> list[dict]:
    pin = doc["library_sha256"]
    broken = pin[:-1] + ("0" if pin[-1] != "0" else "1")
    return [
        {"id": "R01", "kind": "score", "row_id": "L05", "destination": "d0", "amplitude_rule": "copy_score",
         "note": "PAINS-flagged AMPK score asked to replace d0"},
        {"id": "R02", "kind": "score", "row_id": "L12", "destination": "d_u", "amplitude_rule": "copy_score",
         "note": "C3 membrane score asked to replace d_u"},
        {"id": "R03", "kind": "soft_prior", "row_id": "L05", "destination": "d0", "amplitude_rule": "soft_prior",
         "note": "Gaussian prior mean of d0 shifted by the L05 AMPK score"},
        {"id": "R04", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "copy_score",
         "note": "C3 row, numeric score copied into the forcing amplitude"},
        {"id": "R05", "kind": "score", "row_id": "L05", "destination": FORCING_SYMBOL, "amplitude_rule": "affine_score",
         "note": "affine map of the L05 AMPK score proposed as u_phyto"},
        {"id": "R06", "kind": "anecdote", "row_id": None, "destination": FORCING_SYMBOL, "amplitude_rule": "copy_score",
         "note": "roadmap G_tip sentence proposed as a forcing"},
        {"id": "R07", "kind": "score", "row_id": "L05", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "note": "PAINS row asked to authorize the protocol schedule"},
        {"id": "R08", "kind": "score", "row_id": "L08", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "note": "C2 row, clear PAINS flag, asked to authorize the protocol schedule"},
        {"id": "R09", "kind": "utility", "row_id": "L10", "destination": "g", "amplitude_rule": "copy_score",
         "note": "utility of L10 asked to replace g"},
        {"id": "R10", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "library_sha256": broken, "note": "legal-looking call with a broken library pin"},
        {"id": "R11", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "engine": "vina", "note": "engine name set, no engine was run"},
        {"id": "R12", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "may_enter_theta": True, "note": "may_enter_theta flipped"},
        {"id": "R13", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "claim": 4, "note": "claim overwritten to 4"},
        {"id": "R14", "kind": "soft_weight", "row_id": "L10", "destination": "s", "amplitude_rule": "soft_weight",
         "note": "kinetic objective asked to carry a score weight"},
        {"id": "R15", "kind": "score", "row_id": "L10", "destination": "N_ROS", "amplitude_rule": "protocol_constant",
         "note": "protocol schedule aimed at a symbol this ledger does not declare"},
        {"id": "A01", "kind": "score", "row_id": "L10", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "note": "single legal admission"},
        {"id": "R16", "kind": "score", "row_id": "L01", "destination": FORCING_SYMBOL, "amplitude_rule": "protocol_constant",
         "note": "second C3 row after the slot is full"},
        {"id": "R17", "kind": "score", "row_id": "L12", "destination": FORCING_SYMBOL, "amplitude_rule": "copy_score",
         "note": "score copy after admission; must not replace the schedule"},
    ]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def draw_figures(elig: list[dict], calls: list[dict], curves: dict, eq_curve: dict) -> None:
    fig_dir = ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "figure.facecolor": "white",
    })

    # Eligibility. Not a rank plot.
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = [f"{e['id']} {e['label']}" for e in elig]
    colors = ["#1f4e79" if e["eligible"] else "#b9b9b9" for e in elig]
    y = np.arange(len(elig))
    ax.barh(y, [1 if e["eligible"] else 0 for e in elig], color=colors, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Admission predicate (1 = eligible, 0 = refused)")
    ax.set_xlim(0, 1.15)
    ax.invert_yaxis()
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / "eligibility.png", dpi=140)
    plt.close(fig)

    # Trajectories.
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3), sharex=True)
    names = [("G", 0, "glucose-like"), ("R", 1, "stress-like"), ("A", 2, "ATP-like")]
    series = [
        ("null", curves["null"], "#555555", "-"),
        ("admitted forcing", curves["legal"], "#1f4e79", "-"),
        ("score-copy contrast", curves["contrast"], "#a33c3c", ":"),
    ]
    for ax, (key, idx, title) in zip(axes, names):
        for name, curve, color, ls in series:
            ax.plot(curve["t"], curve["y"][idx], color=color, ls=ls, lw=1.6, label=name)
        ax.set_title(title)
        ax.set_xlabel("toy time")
        ax.set_ylabel(key)
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "trajectories.png", dpi=140)
    plt.close(fig)

    # Digest stability across calls.
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    xs = np.arange(len(calls))
    ax.scatter(xs, np.ones(len(calls)), s=28, color="#1f4e79", label="Θ digest unchanged", zorder=3)
    refused_x = [i for i, c in enumerate(calls) if c["status"] != "admitted"]
    admitted_x = [i for i, c in enumerate(calls) if c["status"] == "admitted"]
    ax.scatter(refused_x, np.zeros(len(refused_x)), s=22, color="#9a9a9a", label="forcing not written by this call", zorder=2)
    if admitted_x:
        ax.scatter(admitted_x, np.zeros(len(admitted_x)), s=42, color="#a33c3c", label="forcing schedule written", zorder=4)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["forcing slot", "Θ digest"])
    ax.set_xticks(xs)
    ax.set_xticklabels([c["id"] for c in calls], rotation=70, fontsize=7)
    ax.set_ylim(-0.25, 1.35)
    ax.set_xlabel("ledger call")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=7, loc="center right")
    fig.tight_layout()
    fig.savefig(fig_dir / "digest_stability.png", dpi=140)
    plt.close(fig)

    # Equilibrium versus a constant input. Protocol mark versus a refused amplitude.
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(eq_curve["u"], eq_curve["R"], color="#1f4e79", lw=1.6, label="R equilibrium")
    ax.plot(eq_curve["u"], eq_curve["A"], color="#3d6b4f", lw=1.6, label="A equilibrium")
    ax.axvline(1.0, color="#1f4e79", ls="--", lw=0.9, label="protocol amplitude 1")
    ax.scatter([eq_curve["copy_u"]], [eq_curve["copy_R"]], color="#a33c3c", s=28, zorder=3, label="score-copy amplitude, not admitted")
    ax.set_xlabel("constant u_phyto")
    ax.set_ylabel("equilibrium")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "equilibria.png", dpi=140)
    plt.close(fig)


def snapshot_times(curve: dict) -> dict:
    want = (0.0, 5.0, 10.0, 20.0, 40.0)
    out = {}
    for w in want:
        i = int(np.argmin(np.abs(curve["t"] - w)))
        out[f"{w:.0f}"] = {
            "G": round(float(curve["y"][0, i]), 6),
            "R": round(float(curve["y"][1, i]), 6),
            "A": round(float(curve["y"][2, i]), 6),
        }
    return out


def main() -> None:
    doc = load_candidates()
    if list(THETA) != list(THETA_NAMES):
        raise SystemExit("theta key order drifted")
    elig = eligibility(doc)
    eligible_ids = [e["id"] for e in elig if e["eligible"]]
    if set(eligible_ids) != {"L01", "L10", "L12"}:
        raise SystemExit(f"eligible set {eligible_ids}")

    ledger = {"theta": dict(THETA), "forcing": null_forcing(), "calls": []}
    theta_baseline = digest(ledger["theta"])
    forcing_null = digest(ledger["forcing"])
    theta_bytes = canon_bytes(ledger["theta"])

    for call in build_calls(doc):
        apply_call(doc, ledger, call)

    if canon_bytes(ledger["theta"]) != theta_bytes:
        raise SystemExit("theta bytes changed")
    if digest(ledger["theta"]) != theta_baseline:
        raise SystemExit("theta digest changed")
    admitted = [c for c in ledger["calls"] if c["status"] == "admitted"]
    if len(admitted) != 1 or admitted[0]["id"] != "A01":
        raise SystemExit(f"admission count {len(admitted)}")
    refused = [c for c in ledger["calls"] if c["status"] == "refused"]
    if any(not c["theta_unchanged"] for c in ledger["calls"]):
        raise SystemExit("a call rewrote theta")
    forcing_legal = digest(ledger["forcing"])
    if forcing_legal == forcing_null:
        raise SystemExit("legal forcing digest did not move")
    if ledger["forcing"]["provenance"]["score_copied"] is not False:
        raise SystemExit("provenance says the score was copied")
    # The admitted payload must not contain a transcribed score numeral.
    payload = canon_bytes(ledger["forcing"]).decode("utf-8")
    for row in doc["rows"]:
        for value in row["scores"].values():
            if value in payload:
                raise SystemExit(f"score {value} leaked into the forcing payload")
        if row["utility"] in payload:
            raise SystemExit("utility leaked into the forcing payload")

    # Refusals before admission leave the null forcing in place.
    pre = [c for c in ledger["calls"] if c["id"].startswith("R") and int(c["id"][1:]) < 16]
    if any(c["forcing_sha256_after"] != forcing_null for c in pre):
        raise SystemExit("a refusal wrote the forcing")
    post = [c for c in ledger["calls"] if c["id"] in {"R16", "R17"}]
    if any(c["forcing_sha256_after"] != forcing_legal for c in post):
        raise SystemExit("a later refusal replaced the admitted schedule")

    theta = theta_float()
    curve_null = integrate(theta, 0.0)
    curve_legal = integrate(theta, float(PROTOCOL_AMPLITUDE))
    l10 = row_by_id(doc)["L10"]
    copy_u = -float(l10["scores"]["ampk"])
    curve_copy = integrate(theta, copy_u)
    if curve_null["state_sha256"] == curve_legal["state_sha256"]:
        raise SystemExit("state digest did not move under the admitted forcing")
    if digest(THETA) != theta_baseline:
        raise SystemExit("integration rewrote theta")

    l05 = row_by_id(doc)["L05"]
    prior = soft_prior_contrast(float(l05["scores"]["ampk"]), float(THETA["d0"]))
    if prior["theta_sha256_if_written"] == theta_baseline:
        raise SystemExit("soft prior contrast would not have changed the digest")
    if prior["written_to_ledger_theta"]:
        raise SystemExit("contrast was written")
    if digest(ledger["theta"]) != theta_baseline:
        raise SystemExit("contrast touched the ledger")

    affine_u = affine_amplitude(float(l05["scores"]["ampk"]))
    # Side integration of the affine input. Not a ledger write.
    curve_affine = integrate(theta, affine_u)

    # Equilibrium identities.
    eq0 = equilibrium(theta, 0.0)
    eq1 = equilibrium(theta, 1.0)
    if abs(eq0["G"] - eq1["G"]) > 1e-15:
        raise SystemExit("forcing moved G equilibrium")
    if abs(eq0["R"] - (0.20 * 0.50) / 0.15) > 1e-12:
        raise SystemExit("null R equilibrium")
    if abs(eq1["R"] - 0.25) > 1e-12:
        raise SystemExit("legal R equilibrium")
    if abs(eq1["A"] - 1.0) > 1e-12:
        raise SystemExit("legal A equilibrium")
    if abs(eq0["A"] - 0.75) > 1e-12:
        raise SystemExit("null A equilibrium")

    u_grid = np.linspace(0.0, 2.0, 81)
    eq_curve = {
        "u": u_grid,
        "R": np.array([equilibrium(theta, float(u))["R"] for u in u_grid]),
        "A": np.array([equilibrium(theta, float(u))["A"] for u in u_grid]),
        "copy_u": copy_u,
        "copy_R": equilibrium(theta, copy_u)["R"],
    }
    draw_figures(elig, ledger["calls"], {"null": curve_null, "legal": curve_legal, "contrast": curve_copy}, eq_curve)

    def pack_curve(curve: dict) -> dict:
        return {
            "u": curve["u"],
            "state_sha256": curve["state_sha256"],
            "terminal": {k: round(v, 6) for k, v in curve["terminal"].items()},
            "equilibrium": {k: round(v, 6) for k, v in curve["equilibrium"].items()},
            "equilibrium_abs_gap": {k: float(f"{v:.3e}") for k, v in curve["equilibrium_abs_gap"].items()},
            "closed_form_max_abs": {k: float(f"{v:.3e}") for k, v in curve["closed_form_max_abs"].items()},
            "samples": snapshot_times(curve),
        }

    max_state_gap = float(np.max(np.abs(curve_legal["y"] - curve_null["y"])))
    max_g_gap = float(np.max(np.abs(curve_legal["y"][0] - curve_null["y"][0])))
    results = {
        "rng_draws": 0,
        "candidate_sha256": CANDIDATE_PIN,
        "library_sha256_transcribed": doc["library_sha256"],
        "t16_theta_sha256_cited_not_used": doc["t16_theta_sha256_cited_not_used"],
        "theta": THETA,
        "theta_sha256": theta_baseline,
        "theta_sha256_after_all_calls": digest(ledger["theta"]),
        "theta_unchanged": True,
        "forcing_null_sha256": forcing_null,
        "forcing_admitted_sha256": forcing_legal,
        "admitted_forcing": ledger["forcing"],
        "n_calls": len(ledger["calls"]),
        "n_refused": len(refused),
        "n_admitted": len(admitted),
        "eligible_ids": eligible_ids,
        "eligibility": elig,
        "calls": [
            {
                "id": c["id"],
                "kind": c["kind"],
                "row_id": c["row_id"],
                "destination": c["destination"],
                "amplitude_rule": c["amplitude_rule"],
                "status": c["status"],
                "reasons": c["reasons"],
                "theta_unchanged": c["theta_unchanged"],
                "theta_sha256_after": c["theta_sha256_after"],
                "forcing_sha256_after": c["forcing_sha256_after"],
            }
            for c in ledger["calls"]
        ],
        "curves": {
            "null": pack_curve(curve_null),
            "admitted": pack_curve(curve_legal),
            "score_copy_contrast_not_a_result": pack_curve(curve_copy),
            "affine_contrast_not_a_result": pack_curve(curve_affine),
        },
        "max_abs_state_gap_admitted_vs_null": round(max_state_gap, 6),
        "max_abs_G_gap_admitted_vs_null": float(f"{max_g_gap:.3e}"),
        "contrasts_not_results": {
            "soft_prior": {
                "row_id": "L05",
                "column": "ampk",
                "score": l05["scores"]["ampk"],
                "prior_mean_d0": round(prior["prior_mean"], 6),
                "posterior_mean_d0": round(prior["posterior_mean_d0"], 6),
                "theta_sha256_if_written": prior["theta_sha256_if_written"],
                "written_to_ledger_theta": False,
                "ledger_theta_sha256_after_contrast": digest(ledger["theta"]),
            },
            "score_copy_amplitude": {
                "row_id": "L10",
                "column": "ampk",
                "score": l10["scores"]["ampk"],
                "proposed_u": round(copy_u, 6),
                "written_to_ledger_forcing": False,
            },
            "affine_amplitude": {
                "row_id": "L05",
                "column": "ampk",
                "score": l05["scores"]["ampk"],
                "map": "0.20 + 0.50 * (-score)",
                "proposed_u": round(affine_u, 6),
                "written_to_ledger_forcing": False,
            },
            "anecdote_text": doc["anecdote"]["text"],
            "anecdote_recomputed_here": False,
        },
        "checks": {
            "eligible_set_is_L01_L10_L12": True,
            "single_admission_A01": True,
            "theta_bytes_unchanged": True,
            "refusals_before_admission_left_null_forcing": True,
            "refusals_after_admission_left_admitted_forcing": True,
            "score_numerals_absent_from_forcing_payload": True,
            "closed_form_G_R": True,
            "G_equilibrium_independent_of_u": True,
            "state_digest_moves_when_forcing_is_admitted": True,
            "soft_prior_digest_differs_and_was_not_written": True,
        },
    }
    out = ROOT / "results.json"
    out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"theta_sha256 {theta_baseline}")
    print(f"forcing_null {forcing_null}")
    print(f"forcing_admitted {forcing_legal}")
    print(f"refused {len(refused)} admitted {len(admitted)}")
    print(f"eligible {eligible_ids}")
    print(f"state null {curve_null['state_sha256']}")
    print(f"state legal {curve_legal['state_sha256']}")
    print(f"max gap {max_state_gap:.6f} G gap {max_g_gap:.3e}")
    print(f"prior posterior d0 {prior['posterior_mean_d0']:.6f} sha {prior['theta_sha256_if_written']}")
    print(f"copy u {copy_u:.6f} affine u {affine_u:.6f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
