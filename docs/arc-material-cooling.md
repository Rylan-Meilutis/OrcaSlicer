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

The existing filament **Localized shrinkage mitigation** controls can separately
smooth abrupt layer-time reductions. Arc fan speed and increased-cooling duration
remain independent of the covering-speed limit; use material-appropriate cooling.

These are thermal-process controls, not dimensional compensation. They do not
alter wall coordinates, arc flow, overlap, or temperature. Bulk filament shrinkage
does not predict local contraction of an arc-supported surface. Validate settings
on a physical coupon with both narrow and wide arc regions and measure dimensions
before applying them to tolerance-critical parts. No setting guarantees flatness
for every geometry or material.

## Bundled reference table, version 2

Only exact types `PLA`, `PETG`, `ABS`, `ASA` and `PC` select reference data, not an exact spool grade.
Filled/foamed types and unknown materials are not mapped to the base polymer.
Positive profile density overrides the table; zero uses the reference density.
Nothing is fetched during slicing and estimates never overwrite manual values.

| Reference | Density (g/cm³) | Specific heat (J/kg K) | Conductivity (W/m K) | Tg (°C) |
|---|---:|---:|---:|---:|
| PLA modeling reference | 1.24 | 1800 | 0.20 | 62.8 |
| Eastar 6763 reference for PETG | 1.27 | 1760/1880/1970/2050 at 100/150/200/250°C | 0.21 | 80 |
| ABS generic reference | 1.05 | 1470 | 0.175 | 100 |
| ASA generic reference | 1.055 | 1350 | 0.18 | 100 |
| PC generic reference | 1.22 | 1335 | 0.20 | 145 |

The new generic entries use midpoints of published ranges, not measured melt
properties. ABS and ASA use the rigid matrix transition, not the rubber phase.
They are starting estimates requiring physical tuning, especially for low-fan
or heated-chamber printing. Auto never raises the configured fan speed.

Other catalog materials still use process/manual limits: flexible materials
(TPU/FLEX/SBS), semicrystalline polymers (such as PA/PP/PET/PEEK), filled/foamed
grades, and unspecified blends must not inherit an unrelated base-polymer
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

## Model assumptions

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
Bed heating, shielding, contact conduction, radiation, crystallization, shrinkage
and actual firmware temperatures are not modeled. Tg is a scheduling proxy, not
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
