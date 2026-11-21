# Why the angle trend in Sgr differs between network A and network B

Same rock, same PNflow deck, Pc-controlled drainage, seeds 11/37/50 (24 snapshot runs, Sgi ≈ 0.75: A at 5.8 kPa, B at 8.0 kPa).
All trapped CO2 verified disconnected from both faces. Sgr reproduces the 780-run matrix (A WW 0.315, B WW 0.543).

## What was ruled out
Throat snap-off counts per 1,000 throats are the same on both networks (about 50, 92, 83, ~0 for 20-40, 45-65, 70-90, 95-115 deg).
The snap-off pattern therefore cannot explain opposite Sgr trends.

## Cause
Two local effects with opposite sign, present in BOTH networks:

| Location of the CO2 | Trapped fraction as angle rises 20-40 -> 95-115 deg |
|---|---|
| Throats | A 0.30 -> 0.70, B 0.45 -> 0.60 (rises) |
| Pores with <= 3 throats | A 0.18 -> 0.73, B 0.28 -> 0.72 (rises) |
| Pores with >= 8 throats | A 0.72 -> 0.36, B 0.87 -> 0.19 (falls) |

Cooperative pore-body filling (I4+, a pore filled with >= 4 water-filled neighbours) rises from ~10 events at 20-40 deg to
228 (A) and 336 (B) at 95-115 deg, while snap-off disappears. Highly connected pores are refilled cooperatively and stop trapping;
poorly connected pores and throats are bypassed and trap more.

The networks differ in WHERE the drained CO2 sits (identical at every angle, because drainage uses the same receding angles):

| Share of initial CO2 | A | B |
|---|---|---|
| Throats | 53 % | 9 % |
| Pores with 4-7 throats or fewer | 22 % | 37 % |
| Pores with >= 8 throats | 26 % | 55 % |

A holds most of its CO2 where trapping rises with angle (throats), so Sgr rises. B holds most of it in highly connected pores
where trapping collapses with angle, so Sgr falls.

## Reweighting test (exact: quasi-static fill sequence does not depend on element volume)
Giving A the pore/throat volume split of B: Sgr = 0.41, 0.38, 0.39, 0.39 (A's rise disappears).
Giving B the pore/throat volume split of A: Sgr = 0.40, 0.42, 0.42, 0.31 (B's fall is reduced to a drop at 95-115 deg).
Throat volume share explains A's rise; B's fall comes from the highly connected pores.

## Limits
- Pore-body filling counts come from the PNflow log; the link to trapping efficiency is a measured association, not a proof of mechanism.
- One Sgi (~0.75), three seeds per case; B's fit region has few points near Sgi 0.8.
- Throat-volume share depends on how the extraction assigns volume between pores and throats; that choice is part of the network, not of the rock.

Scripts: run_structure_snapshots.py, analyse_trapping_structure.py, analyse_trapping_by_size.py, analyse_cause_tests.py, plot_cause_A_vs_B.py, parse_filling_events.py.
Tables: trapping_structure_mean.csv, cause_test_reweighting.csv, cause_test_pore_efficiency_by_coordination.csv, filling_events_all_runs.csv. Figure: fig11_cause_trapping_A_vs_B.png.

## Geometry sensitivity (numerical experiment; 432 runs, Sgi ~0.60, 3 seeds, one knob at a time)
Metric: delta = Sgr/Sgi(95-115 deg) - Sgr/Sgi(20-40 deg); delta > 0 means less water-wet traps more.

| Knob | Network A (delta at baseline +0.17) | Network B (delta at baseline -0.25) |
|---|---|---|
| Coordination (delete 10-40% of throats, z 3.2 -> 1.9 / 3.8 -> 2.3) | stays positive (0.17 -> 0.11) | stays negative, steeper (-0.25 -> -0.44) |
| Pore/throat radius ratio (scale throat radii) | FLIPS at ratio ~5.0 (baseline 2.95); -0.08 at ratio 8.4 | approaches zero (-0.04 at ratio 2.0), no flip |
| Throat shape factor G (scale G) | stays positive (0.09-0.22) | approaches zero (-0.03 at G = 0.063), no flip |

Conclusions
- The trend is controlled mainly by the pore-to-throat radius ratio: narrow throats and wide pores (ratio > ~5) give the classical order (more water-wet traps more); wide throats give the reverse.
- Throat angularity matters in network B (angular throats strongly favour water-wet trapping); rounder throats remove B's trend.
- Coordination number alone does not flip either network (an earlier suspicion that highly connected pores decide it is not confirmed as the controlling lever).
- Network B cannot be pushed to A's behaviour with these knobs; its flip also needs a larger throat share of pore volume (~78 % on its own fill sequence, see fig12). The two properties (radius ratio, volume split) act together.
- Edited networks are not real rocks. Throat volumes were held fixed when radii/shape were scaled to isolate the fill sequence; one network realisation per deletion level and seed.

Files: run_geometry_sweep.py, make_variant_network.py, summarise_geometry_sweep.py, geometry_sweep_table.csv, geometry_sweep_runs.jsonl, fig13_geometry_sweep_flip.png.

## Combined test: can A act "normally" at a realistic permeability? (84 runs, Sgi ~0.60, K held at 1,140 mD)
Throat length scaled to hold K at B's value. delta = Sgr/Sgi(95-115) - Sgr/Sgi(20-40); B's own delta at Sgi 0.6 is -0.25.

| Variant of A | aspect | G | throat vol share | K (mD) | delta |
|---|---|---|---|---|---|
| K-matched control | 2.96 | 0.044 | 54 % | 1147 | +0.16 |
| volume share 11 % | 2.96 | 0.044 | 11 % | 1147 | -0.01 |
| radius x0.78 + share 11 % | 3.78 | 0.044 | 11 % | 1139 | -0.11 |
| radius x0.78 + G x0.69 + share 11 % (B-like) | 3.78 | 0.030 | 11 % | 1121 | -0.22 |
| radius x0.65 + share 11 % | 4.54 | 0.044 | 11 % | 1137 | -0.23 |
| radius x0.65 + G x0.69 + share 11 % | 4.54 | 0.030 | 11 % | 1138 | -0.31 |
| radius x0.55 + share 11 % | 5.37 | 0.044 | 11 % | 1135 | -0.31 |

A flips to the classical trend at realistic K with two moderate changes (aspect 3.0 -> 3.8, throat volume share 54 % -> 11 %); adding B's throat angularity reproduces B's trend (-0.22 vs -0.25).
Files: run_combo_test.py, combo_test_table.csv, combo_test_runs.json.
