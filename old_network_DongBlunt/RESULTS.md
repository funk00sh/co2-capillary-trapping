# Berea CO₂–water drainage and imbibition: trapping and hysteresis (12 Nov 2026)

306 quasi-static runs of the unmodified upstream PNflow executable on the Dong & Blunt Berea network. All 306 exited with code 0.

## Set-up

| Item | Value |
|---|---|
| Network | Imperial College Berea (Dong & Blunt): 6,298 pores, 12,545 throats, 2.138 mm cube, porosity 19.6%, 1,111 mD |
| Fluids | Water and supercritical CO₂ at 50 °C, 10 MPa; IFT 37.37 mN/m (base), 25 and 50 mN/m (sensitivity) |
| Primary drainage | Receding angle 0–10°; stopped at a prescribed maximum capillary pressure, so Sgi is an output |
| Pc ladder | 18 levels, 5.8–50 kPa (equivalent CO₂ column 0.97–8.4 m, Δρ = 608 kg/m³) |
| Imbibition | Advancing angle 20–40°, 45–65° or 70–90° (through water); run to Pc = −200 kPa |
| Boundaries | Inject left, produce both faces (upstream sample convention); averages over x = 0.1–0.9 |
| Repeats | 5 random seeds per wettability and Pc level; IFT cases use one seed |

Differences from the 10 November pilot: drainage no longer forced to Sgi ≈ 1; separate receding and advancing angles; production from both faces; end regions excluded from averages; a full Sgi ladder instead of one reversal.

## Results

**Trapping curve (fig3, trapping_table.csv, land_fits.json).** Land model Sgr = Sgi / (1 + C·Sgi):

| Advancing angle | Land C | Fit RMSE | Largest Sgr |
|---|---:|---:|---:|
| 20–40° | 0.54 ± 0.02 | 0.018 | 0.65 |
| 45–65° | 0.73 ± 0.02 | 0.018 | 0.56 |
| 70–90° | 1.27 ± 0.07 | 0.041 | 0.47 (at Sgi 0.80) |

- Trapping falls as the rock becomes less water-wet. Below Sgi ≈ 0.25 the three cases coincide.
- The 70–90° case is non-monotonic: Sgr peaks at 0.47 near Sgi 0.80 and falls to 0.35 at Sgi 0.99. Land's model does not describe this, which is why its RMSE is larger.
- Berea CO₂–brine corefloods are usually summarised by Land C ≈ 1 (Krevor et al. 2012; Niu et al. 2015 use 1.1). The model brackets this: C = 1 lies between the 45–65° and 70–90° cases. The strongly water-wet cases trap more than the experiments.

**kr hysteresis (fig1).** Base case 45–65°, Pc,max 8 kPa: Sgi 0.76 → Sgr 0.51. Imbibition krg lies below drainage krg at the same Sw and reaches zero at Sgr; imbibition krw lies above drainage krw. Seed-to-seed scatter is small.

**Pc hysteresis (fig2).** Imbibition scanning curves from five reversal points all sit well below primary drainage and end at successively larger trapped saturations.

**Wettability (fig4).** At the same Sgi, the 70–90° case keeps CO₂ mobile to higher Sw and ends with a much lower krw. Trapping efficiency Sgr/Sgi is 0.65–0.75 for the strongly water-wet case and falls to 0.35 for 70–90° at high Sgi.

**IFT (fig5).** Pc scales in proportion to IFT. At a fixed applied Pc this changes Sgi strongly: at 8 kPa, Sgi is 0.92, 0.76 and 0.10 for 25, 37.4 and 50 mN/m. The Sgr-versus-Sgi relationship is the same for all three IFT values. In a capillary-dominated model, IFT therefore controls how much CO₂ enters, not what fraction of it is trapped.

## Limitations

- **Network size.** A larger Berea network could not be obtained: the Imperial College download is this same network, the 12,349-pore Øren–Bakke network was not found in a public repository, and this PC cannot extract one from a larger image. Gas breakthrough occurs at Sw ≈ 0.59 in every run and the drainage krg is stepped; both are finite-size features. Seeds change only the contact-angle assignment, not the network.
- **No sub-resolution porosity.** Sw falls below 0.05 by about 15 kPa, so the Sgi points above ~0.85 are not reachable in a real Berea core. The Pc curve is flat, so Sgi is very sensitive to the applied Pc between 6 and 9 kPa.
- **Low-Sgi points.** Runs with Sgi below ~0.4 reversed before gas breakthrough. Their Sgr values are valid; their krg is zero by construction.
- **Quasi-static.** No viscous or rate effects, dissolution or reaction. IFT cannot affect trapping through capillary number here.
- **Contact angles are assumed ranges**, not measurements, and are not correlated with pore size.
- **Unexplained diagnostic.** PNflow's snap-off share of throat fillings is lowest for the 20–40° case (column `snap_frac` in trapping_table.csv), opposite to the simple expectation. It was not investigated and is not used in any figure.

## Files

- `run_matrix.py`, `analyse.py`: reproduce runs and figures.
- `summary_all_runs.csv`: one row per run. `trapping_table.csv`: seed means. `land_fits.json`: fits.
- `fig1`–`fig5` PNG. `raw_runs_306.zip`: PNflow decks, logs and curve tables for every run.

## Sources

- Krevor, Pini, Zuo & Benson (2012), WRR, https://doi.org/10.1029/2011WR010859
- Niu, Al-Menhali & Krevor (2015), WRR, https://doi.org/10.1002/2014WR016441
- Valvatne & Blunt (2004), WRR, https://doi.org/10.1029/2003WR002627
- PNflow: https://github.com/ImperialCollegeLondon/pnflow
- Berea network: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/berea-sandstone/
