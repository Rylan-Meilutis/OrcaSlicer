# Material-specific arc cooling

Enable **Automatic arc cooling estimates** under Filament → Cooling → Overhang
and bridge cooling to derive limits from the material. The UI shows an
**Auto — estimated** reference for a 0.45 mm round strand. Slicing uses the arc's
equivalent round diameter from deposited volume per millimeter. Covering layers
use a reference diameter of 1.125 times the selected nozzle diameter.

Two advanced overrides remain alongside **Arc overhang fan speed**:

- **Minimum arc path time** slows each native arc path, not each fitted segment.
  The larger of the resolved filament and process minimum times is used. Positive
  filament values override Auto; zero uses Auto when enabled, otherwise the
  process setting. Automatic estimates respect the filament's existing minimum
  print speed (`slow_down_min_speed`), except paths shorter than four equivalent
  strand diameters gradually relax that floor. Automatic arc slowdown cannot
  more than double each path's otherwise-planned time. This is a bounded heuristic,
  not a physical cooling guarantee. Explicitly configured minimum times retain
  the 0.1 mm/s floor. Neither can increase a speed already limited by flow or
  another process setting. This does not
  insert a cooling dwell or guarantee a minimum time between nearby deposits.
- **Maximum arc covering speed** caps deposition on subsequent layers for the
  process's configured arc-overhang/arc-bridge speed-transition duration. Unlike
  the process speed ramp, a manual ceiling stays constant throughout that duration.
  It applies to the whole object, including its walls, using the filament being
  deposited. Positive values override Auto; zero uses Auto when enabled, otherwise
  leaves the process behavior unchanged. A zero transition-layer
  count disables the corresponding covering window.

Automatic covering-speed estimates also respect the filament minimum print speed.
They recover with the configured layer-speed ramp and no longer cap the final
recovery layer; explicitly set covering-speed limits remain fixed.
The estimated cooling time is a soft target, not a dwell required for every small
arc: material properties and an assumed convection coefficient cannot establish
a calibrated deposition speed. Reaching the speed floor may leave the target
unmet. Manual covering-speed overrides can deliberately request lower speeds.

These settings are saved in filament presets and exported with the job. Old
profiles and projects inherit Auto off and zero overrides, preserving their existing output.
Changes require G-code regeneration but not regeneration of arc geometry.

The existing filament **Thermal transition smoothing** controls can separately
smooth abrupt layer-time reductions. Arc fan speed and increased-cooling duration
remain independent of the covering-speed limit; use material-appropriate cooling.

These are thermal-process controls, not dimensional compensation. They do not
alter wall coordinates, arc flow, overlap, or temperature. Bulk filament shrinkage
does not predict local contraction of an arc-supported surface. Validate settings
on a physical coupon with both narrow and wide arc regions and measure dimensions
before applying them to tolerance-critical parts. No setting guarantees flatness
for every geometry or material.

## Bundled reference table, version 5

Only exact types `PLA`, `PETG`, `ABS`, `ASA`, `PC`, `HIPS`, `PMMA`, `PSU`
and `PES`/`PESU` select reference data, not an exact spool grade. `PP`, `PPS`
and `PPS-CF` additionally use the semi-crystalline scheduling approximation below.
Except for the explicitly identified PPS-CF matrix proxy, filled/foamed types and
unknown materials are not mapped to the base polymer.
Positive profile density overrides the table; zero uses the reference density.
Nothing is fetched during slicing and estimates never overwrite manual values.

| Reference | Density (g/cm³) | Specific heat (J/kg K) | Conductivity (W/m K) | Tg (°C) |
|---|---:|---:|---:|---:|
| PLA modeling reference | 1.24 | 1800 | 0.20 | 62.8 |
| Eastar 6763 reference for PETG | 1.27 | 1760/1880/1970/2050 at 100/150/200/250°C | 0.21 | 80 |
| ABS generic reference | 1.05 | 1470 | 0.175 | 100 |
| ASA generic reference | 1.055 | 1350 | 0.18 | 100 |
| PC generic reference | 1.22 | 1335 | 0.20 | 145 |
| HIPS combined reference | 1.04 | 1351 | 0.128 | 101.1 |
| PMMA generic, atactic reference | 1.17 | 1460 | 0.205 | 105 |
| PSU generic reference | 1.245 | 1370 | 0.15 | 187 |
| PES/PESU polyethersulfone reference | 1.37 | 1370 | 0.18 | 228 |

The ABS, ASA and PC entries use midpoints of published ranges, not measured melt
properties. ABS and ASA use the rigid matrix transition, not the rubber phase.
They are starting estimates requiring physical tuning, especially for low-fan
or heated-chamber printing. Auto never raises the configured fan speed.

Other catalog materials still use process/manual limits: flexible materials
(TPU/FLEX/SBS), other semicrystalline polymers (such as PA/PET/PEEK), filled/foamed
grades other than the explicit PPS-CF proxy, and unspecified blends must not inherit an unrelated base-polymer
estimate. Density alone cannot determine specific heat or crystallization.

Sources checked 2026-09-21:

- NETZSCH generic ranges for [ABS](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/commodity-thermoplastics/abs-acrylonitrile-butadiene-styrene-copolymer),
  [ASA](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/commodity-thermoplastics/asa-acrylonitrile-styrene-acrylate-copolymer),
  and [PC](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/engineering-thermoplastics/pc-polycarbonate).

- PLA density, heat capacity and conductivity: Table 1 of
  [Mathematical Model of the Layer-by-Layer FFF/FGF Polymer Extrusion Process](https://pmc.ncbi.nlm.nih.gov/articles/PMC10707973/).
  These are modeling constants, not measurements of every PLA formulation.
- PLA transition proxy: Table 2 (unfilled 4032D) of
  [Blanco et al., Entropy 2022, 24, 654](https://doi.org/10.3390/e24050654).
  Combining this grade-specific Tg with generic constants is an approximation.
- PETG reference: [Eastman Eastar 6763 technical data sheet](https://productcatalog.eastman.com/tds/ProdDatasheet.aspx?pn=Eastar+676&product=71040786).
  Heat capacity is interpolated and clamped at the tabulated temperature endpoints,
  not claimed to represent molten properties outside the measured range.

HIPS sources checked 2026-10-17:

- [UL FSRI HIPS measurements](https://materials.fsri.org/materialdetail/high-impact-polystyrene-hips):
  specific heat at 40°C and conductivity at 45°C. These fixed solid-state
  references are not molten-property measurements or temperature extrapolations.
- [TotalEnergies Impact 3450 data sheet](https://polymers.totalenergies.com/sites/g/files/wompnd5016/files/site_collection_documents/Technical%20Datasheets/3450-EU-EN.pdf):
  reference density, overridden by positive profile density.
- [Sieradzka et al., Materials 2021, 14, 7008](https://pmc.ncbi.nlm.nih.gov/articles/PMC8623337/),
  section 3.4: 101.1°C matrix transition for unfilled HIPS filament made from
  Impact 3450, not its graphene-filled variants. Combining this with FSRI data
  from another sample is an approximation requiring physical tuning; it does
  not establish a universal HIPS cooling curve.

## Model assumptions

### Semi-crystalline scheduling references

| Type | Density (g/cm³) | cp (J/kg K) | k (W/m K) | Cooling target (°C) | Latent heat allowance (J/kg) | Release temperature (°C) |
|---|---:|---:|---:|---:|---:|---:|
| PP | 0.905 | 1800 | 0.21 | 100 | 208000 | 120 |
| PPS | 1.35 | 1830 | 0.30 | 90 | 150400 | 205 |
| PPS-CF matrix proxy | 1.35 | 1830 | 0.30 | 90 | 150400 | 205 |

PP properties and full-crystal fusion enthalpy use
[NETZSCH isotactic PP](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/commodity-thermoplastics/pp-isotactic-polypropylene).
The 100°C target and 120°C heat-release point are scheduling assumptions, not
Tg or predictions of load-bearing strength. PPS uses the unfilled column of
[Celanese's Fortron design manual, Table 1](https://www.celanese.com/-/media/Engineered-Materials/Files/Product-Technical-Guides/PPS-050_Fortron_DesignManualTG_AM_0913.pdf):
melt heat capacity and conductivity at 300°C. Its 205°C release point is the
midpoint of the reported 160–250°C grade-dependent recrystallization range.
The 150.4 J/g full-crystal allowance is the high reference discussed in
[Spruiell and Janke, ORNL](https://info.ornl.gov/sites/publications/Files/Pub57451.pdf).
These deliberately include the full latent-heat allowance rather than claiming
to know the crystallinity reached during fast printing.

PPS-CF is explicitly a **matrix proxy**, not measured composite data: positive
profile density wins, but fiber fraction, orientation and anisotropic heat
transport are unknown. Its UI warns to calibrate the existing manual overrides.
No arbitrary suffix matching maps other composites to this proxy.

The additional heat-release time is `rho*L*d / (4*h*(Trelease-Tambient))`,
only when the nozzle temperature exceeds the release temperature. The same
30-second validity ceiling, minimum-speed constraints and 2× automatic arc
slowdown limit apply. This is not a crystallization-kinetics simulation.

Bulk XY/Z shrinkage compensation is **not** inferred from this cooling table.
Crystallization, constraints and fiber orientation can change shrinkage and
warpage independently of heat capacity. Existing measured compensation values
remain authoritative; injection-molding shrinkage is not an FFF correction.

Additional references checked 2026-10-17: NETZSCH
[PMMA](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/engineering-thermoplastics/pmma-polymethylmethacrylate),
[PSU](https://analyzing-testing.netzsch.com/en-AU/polymers-netzsch-com/high-temperature-resistant-thermoplastics/psu-polysulfone),
and [PESU](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/high-temperature-resistant-thermoplastics/pesu-polyethersulfone).
PMMA uses midpoints of the density, heat-capacity and conductivity ranges and
the atactic transition. PSU uses midpoint density and the measured second-heating
transition; PESU uses its measured transition. These fixed reference properties
are not temperature-dependent melt measurements. `PES` is the filament catalog's
name for polyethersulfone, not a catch-all for polyester.

Generated calibration projects explicitly disable arc overhangs on their model
objects and any overriding volumes or height ranges. This preserves conventional
bridges in temperature towers and other calibration models without changing the
saved process preset. Existing calibration projects must be regenerated to receive
the override.

Arc fill uses the round strand width as its nominal pitch, without the extra gap
used for straight bridges. Straight-bridge density does not alter arc pitch;
the separate arc-overlap control reduces it to establish bead contact. Native
arc generation reserves a minimum tessellation/fitting allowance (at most 5% of
strand width) when the requested overlap is smaller, so polygonized neighbouring
courses remain in contact. Otherwise,
the support checks reject separated concentric courses and recursive filling
replaces the long courses with many tiny restart families. The extrusion's round
cross-section and volume are unchanged by this pitch correction.

These estimates do not equalize the interval between successive wall deposits.
Arc-bearing regions now print the local supporting wall sections first, then
arcs, then the remaining retained wall sections and ordinary fill in the configured
order: walls before infill unless infill-first is enabled. The anchor
neighborhood covers the arc endpoints plus the configured wall-stack width
(at least 1 mm for a lead-in). Complementary clipping preserves wall geometry,
flow and height without double extrusion. Contact walls still print early;
putting them last could leave an arc unsupported. This conservative neighborhood
can retain more wall than strictly necessary, especially on narrow features.
Other regions retain their ordinary ordering. This reduces the extra cooling
interval for deferred wall sections, but is not a guarantee of dimensional accuracy.

The scheduling heuristic adds a lumped convection time and radial diffusion scale:

`t = rho*cp*d/(4*h) * ln((Tnozzle-Tambient)/(Tg-Tambient)) + rho*cp*d²/(16*k)`

Units are SI. Ambient is the larger of 25°C and the configured chamber temperature.
The **uncalibrated airflow assumption** is `h = 20 + 80*fan_percent/100 W/(m² K)`
using the arc fan setting. PWM is not a measured heat-transfer coefficient.
Bed heating, shielding, contact conduction, radiation, crystallization kinetics,
shrinkage and actual firmware temperatures are not modeled. PP/PPS add the lumped
heat allowance above and use their listed cooling targets in place of Tg.
Tg is a scheduling proxy, not
proof a strand can carry load.

Time has a one-second lower bound. Cover speed is `15 mm / t`, bounded to
0.5–15 mm/s. The 15 mm reference traverse is a tuning heuristic, not a sourced
material property. Existing process and volumetric limits may slow it further.

Unknown types, invalid inputs, zero fan, ambient within 5°C of Tg, density outside
0.5–2 g/cm³, diameter above 2 mm, or times exceeding 30 s produce **no estimate**.
The process/manual settings remain in effect. Excessive estimates are not clipped
and called adequate cooling. Auto is also skipped during initial fan-off layers.
The UI reports unavailable estimates for its reference
strand. Other materials require documented data and physical validation.
