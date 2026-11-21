# Berea CO₂–water drainage and imbibition: relative permeability, hysteresis and trapping (12 Nov 2026)

462 quasi-static runs of the unmodified upstream PNflow executable on a Berea network re-extracted from the Imperial College micro-CT image. All 462 exited with code 0.

## Why the network was replaced

The first set of runs (kept in `old_network_DongBlunt/`) used the published 2007 Dong & Blunt network. Its kr curves were poor for two reasons that belong to that network, not to the solver or the settings:

| Symptom | Old network (6,298 pores) | New network (4,670 pores) |
|---|---|---|
| Gas breakthrough in drainage | Sw ≈ 0.59 | Sw ≈ 0.80 |
| Imbibition krw / drainage krw at equal Sw | 2–3 | about 1 |
| Sgr from Sgi 0.75, θa 45–65° | 0.51 | 0.34 |
| Throat share of pore volume | 11% | 54% |

In the old network almost all volume sits in pore bodies and throats are thin connectors. Gas must fill 40% of the pore volume before it spans the sample, and during imbibition water refills the throats (large gain in conductance) while gas stays in the bodies (small change in saturation), which produces the large apparent water hysteresis. The same test on the Bentheimer network already in the project gave breakthrough at Sw 0.80 and a ratio of 1.2–1.6, which confirmed the cause.

The new network was extracted with Imperial College's current `pnextract` tool from the same Berea image (400³ voxels, 5.345 µm, porosity 19.65%). Files and extraction log are in `network_pnextract/`.

## Set-up

| Item | Value |
|---|---|
| Network | 4,670 pores, 7,933 throats, 2.138 mm cube, porosity 20.1%, 1,516 mD |
| Fluids | Water and supercritical CO₂ at 50 °C, 10 MPa; IFT 37.37 mN/m (base), 25 and 50 mN/m |
| Primary drainage | Receding angle 0–10°; stopped at a prescribed maximum Pc, so Sgi is an output |
| Pc ladder | 21 levels, 4.4–30 kPa (equivalent CO₂ column 0.74–5.0 m, Δρ = 608 kg/m³) |
| Imbibition | Advancing angle 20–40°, 45–65°, 70–90° or 95–115° (through water); run to Pc = −200 kPa |
| Boundaries | Inject left, produce both faces; averages over x = 0.1–0.9 |
| Repeats | 5 random seeds per wettability and Pc level; IFT cases use one seed |

## Results

**Relative permeability (fig1, fits.json).** Base case θa 45–65°, drainage to 5.8 kPa: Sgi 0.75 → Sgr 0.35.

| Quantity | Value |
|---|---:|
| Critical gas saturation in drainage | 0.19 |
| krg at Swi | 0.55 |
| krw at Sgr | 0.23 |
| Corey exponent, water (drainage / imbibition) | 1.9 / 1.8 |
| Corey exponent, CO₂ (drainage / imbibition) | 1.5 / 1.2 |

- CO₂ shows clear hysteresis: imbibition krg falls below drainage krg and reaches zero at Sgr.
- Water shows almost none in the water-wet cases: imbibition krw follows drainage krw. This is the expected behaviour for the wetting phase.

**Trapping curve (fig3, trapping_table.csv).** Land model Sgr = Sgi / (1 + C·Sgi), fitted for Sgi ≤ 0.8:

| Advancing angle | Land C | Fit RMSE | Sgr at Sgi 0.75 |
|---|---:|---:|---:|
| 20–40° | 1.92 ± 0.08 | 0.019 | 0.32 |
| 45–65° | 1.71 ± 0.07 | 0.017 | 0.34 |
| 70–90° | 1.35 ± 0.04 | 0.014 | 0.38 |
| 95–115° | 0.96 ± 0.03 | 0.015 | 0.44 |

- Berea CO₂–brine corefloods are usually summarised by Land C ≈ 1 (Krevor et al. 2012; Niu et al. 2015 use 1.1), i.e. Sgr ≈ 0.43 at Sgi 0.75. The water-wet cases trap somewhat less than that.
- Above Sgi ≈ 0.8, Sgr falls again in every case. Those states need more than ~6.5 kPa, where the network has almost no water left; they are shown as open symbols and excluded from the fits.

**Wettability (fig4).** On this network trapping rises as the advancing angle increases. On the old network it fell. The direction of the contact-angle effect therefore depends on the pore-space representation and should not be presented as a general result. What is consistent in both: Sgi is the dominant control, and the wettability effect within the water-wet range is modest (Sgr 0.32–0.38 at Sgi 0.75 here). Less water-wet cases have much lower imbibition krw.

**Pc hysteresis (fig2).** Imbibition scanning curves from five reversal points lie well below primary drainage.

**IFT (fig5).** Pc scales in proportion to IFT. At an applied Pc of 5.8 kPa, Sgi is 0.93, 0.75 and 0.05 for 25, 37.4 and 50 mN/m. The Sgr-versus-Sgi relationship is the same for all three.

## Audit

- 462 of 462 runs exit code 0. Sgi, Sgr and kr all within [0, 1]; no case with Sgr > Sgi.
- Breakthrough Sw between 0.795 and 0.803 in all runs that reach it.
- Seed standard deviation of Sgr is typically 0.01–0.02 (up to 0.04 for 70–90°).

## Limitations

- **One 2.1 mm sample.** Roughly 17 pores across. The drainage krg still has a shelf between Sw 0.70 and 0.80, a finite-size feature that the Corey fit smooths over. Seeds vary only the contact-angle assignment.
- **Network extraction matters.** Two extractions of the same image give Sgr 0.34 and 0.51 for the same case, and opposite wettability trends. Absolute permeability is 1,516 mD here against 1,111 mD for the old network.
- **Connectivity.** 18% of pores in the new network are isolated and 17% are dead ends; mean coordination number 3.3.
- **No sub-resolution porosity.** Sw falls below 0.05 by about 10 kPa. The Pc curve is flat, so Sgi is very sensitive to Pc between 4.7 and 6 kPa.
- **Low-Sgi points.** Runs with Sgi below 0.2 reversed before gas breakthrough; their Sgr is valid but they have no krg.
- **Quasi-static.** No viscous or rate effects, dissolution or reaction.
- **Contact angles are assumed ranges**, not measurements.

## Files

- `run_matrix.py`, `analyse.py`: reproduce runs and figures.
- `summary_all_runs.csv`: one row per run. `trapping_table.csv`: seed means. `fits.json`: Land and Corey fits.
- `fig1`–`fig5` PNG. `raw_runs_462.zip`: PNflow decks, logs and curve tables for every run.
- `network_pnextract/`: network files, image header and extraction log.
- `old_network_DongBlunt/`: the superseded first set (306 runs).

## Sources

- Krevor, Pini, Zuo & Benson (2012), WRR, https://doi.org/10.1029/2011WR010859
- Niu, Al-Menhali & Krevor (2015), WRR, https://doi.org/10.1002/2014WR016441
- PNflow and pnextract: https://github.com/ImperialCollegeLondon/pnflow
- Berea image and original network: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/berea-sandstone/
