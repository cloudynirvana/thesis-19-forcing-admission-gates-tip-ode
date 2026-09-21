# Forcing admission under evidence gates: which phytochemical screen scores may enter a tip ODE as known forcings?

**Thesis #19.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

**Depends on:** [Thesis #1](https://github.com/cloudynirvana/thesis-01-confluence-onco) (evidence gates: knowledge ≠ Θ), [Thesis #7](https://github.com/cloudynirvana/thesis-07-tnbc-tipping-identifiability) (TNBC tip ODE under known forcings), [Thesis #16](https://github.com/cloudynirvana/thesis-16-mitochondrial-dpsim-phytochemical-screen) (gated phytochemical / ΔΨm screen; scores ≠ Θ).

When a phytochemical screen is barred from writing scores into Θ, which of those scores—if any—may be admitted as declared known forcings into a tip ODE without bypassing evidence gates or leaking through soft priors?

Twelve claim records are transcribed from Thesis #16 and are not re-ranked. Eighteen ledger calls are issued. Seventeen return refused. One call, A01, admits row L10 as the provenance of a protocol-constant forcing `u_phyto = 1` on a bounded three-state demonstration field. The SHA-256 of kinetic Θ is `1ace2bb9eed19f85046f25ab690076d0d5920cf19219ddb362c775be9f41d0e8` before the refusals and the same string after the admission. The forcing digest changes. The glucose-like coordinate does not see the forcing. The stress-like equilibrium moves from 0.666667 to 0.250000. A soft-prior contrast that would have set `d0` to 0.186979 is computed beside the ledger and is not written.

The field is not the Thesis #7 right-hand side, and no Fisher table is recomputed. The scores are not a new screen. An admitted schedule is not a dose and not an efficacy.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/admission.py` | Admission ledger and bounded tip-demonstration field |
| `sim/candidates.json` | Transcribed Thesis #16 claim records (not recomputed) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Eligibility, trajectories, digests, equilibria |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/admission.py
python3 build_pdf.py
```

NumPy, SciPy, and Matplotlib are required for the ledger. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`. The candidate-file SHA-256 is pinned inside `sim/admission.py`. A byte edit that does not update the pin raises.

## Cite

Ogbonna KE. Forcing admission under evidence gates: which phytochemical screen scores may enter a tip ODE as known forcings? [Internet]. Thesis #19 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-19-forcing-admission-gates-tip-ode

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
