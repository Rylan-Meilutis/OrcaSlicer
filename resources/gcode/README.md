# INDX / RME custom G-code templates

These templates preserve the custom Cygnas X1 INDX startup and tool-change settings.
They require the corresponding RME firmware commands (including `M976` and `M870`)
and a calibrated INDX cleaner. They are not generic Prusa printer defaults.

Startup requires RME firmware supporting the explicit physical-tool list in
`G427 T... R2 P3`. The list comes from `is_extruder_used`: a job using tools 0
and 7 emits `G427 T0,7 R2 P3`. This template assumes Orca tool IDs correspond
directly to physical INDX tool IDs; do not use that list unchanged with remapping.
Older firmware without this parameter may still calibrate all tools.

- `indx_rme_start.gcode`: printer preset **Machine start G-code**.
- `indx_rme_toolchange.gcode`: printer preset **Change filament G-code**.

These files do not overwrite saved or cloud printer presets. Copy each template
into its corresponding preset field, save the preset, and re-slice the project.

The bundled **Prusa CORE One INDX RME 0.4 nozzle** preset is based on the
Cygnas X1 INDX profile, with eight independent tools and the 248 × 205 mm bed.
It deliberately contains no OctoPrint address or API key. Select the model in
the printer setup wizard and configure its connection separately.

Its adaptive-mesh fallback uses Orca's first-layer bounds (including supports,
brims and wipe towers), clamped to the bed with a 1 mm margin. Set **Bed mesh
min** to **0,0**, **Bed mesh max** to **248,205**, and **Adaptive bed mesh margin**
to **1 mm** when copying the standalone start template into another profile.
The start sequence restores these bounds with M555 after calibration and
immediately before ordinary G29 P1. Native probe spacing and safety margins
remain in effect; small corner prints do not request a full-bed rectangle.

RME Compatibility b118+ and firmware advertising mesh_area=1 can supply
file-derived bounds through their capability-gated MESH service before P1.
The preset intentionally does not embed serial-only @RME commands or explicit
G29 XY/size arguments, so that host override remains available and file printing
can use the M555 fallback. Verify the probe area on hardware after re-slicing;
the slicer cannot guarantee a particular probe count.

Bucket priming replaces a bed purge line. After purging, each sequence retracts,
waits for motion completion, wipes with `G12 S1`, waits again, ejects the new
pellet with `G12 S30`, and exits with `G12 S91`. Exit alone is not a wipe.
The initial leftover-pellet ejection also follows a quick wipe.

Startup homes XY once, before vent control and tool pickup, then homes Z after
tool pickup. Do not add another XY or all-axis home later in the template.

Auto-PA (`M976`) runs after homing and before tool-offset calibration (`G427`).
The bed heats during PA; the template waits for its requested temperature
before offset calibration and restores the initial print tool only after both
calibration stages. No startup command turns off the bed. Optional chamber
warm-up still raises the bed target and then restores the printing target;
its `M104 S0` turns off the **nozzle**, not the bed, and is retained.
Heater shutdown at the end of the print must remain in the end G-code.

Do not insert `[filament_start_gcode]` into the change-filament template: Orca
already evaluates and emits the selected filament's start G-code after the tool
change. Raw substitution duplicates it and can leave expressions in the output.

The `[INDX]` tests exercise startup with single, nonconsecutive, and all tools and tool changes with the
wipe tower enabled/disabled. Physical pellet release and nozzle priming still
need verification on the printer. Existing exported G-code must be re-sliced.
